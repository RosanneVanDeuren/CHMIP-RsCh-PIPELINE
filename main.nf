#!/usr/bin/env nextflow

/*
 * Author: Rosanne C. van Deuren
 *
 *   This is the Nextflow main to run 'CHMIP-RsCh-PIPELINE'.
 *
 *   CHMIP-RsCh-PIPELINE is a pipeline designed to filter candidate clonal hematopoiesis
 *   variant calls as identified by the commercial software JSI Medical Systems Sequence Pilot
 *   using Molecular Inversion Probe (MIP) research data for the purpose of conducting cohort based analyses.
 *
 *   For more information or questions please contact: rc.vandeuren@outlook.com
 */

// Message.
println """
Running::

      ===  ==  == ==    ==  ==  ====        ====             ====             ====    ==  ====    ======  ==      ==  ==   ==  ======
   //      ||  || ||\\\\//||  ||  ||  \\\\      ||  \\\\    **   //       *         ||  \\\\  ||  ||  \\\\  ||      ||      ||  |\\  ||  ||
  ||       ||__|| ||    ||  ||  ||  //  --  ||  //   *    ||        *     --  ||  //  ||  ||  //  ||__    ||      ||  ||\\ ||  ||__
  ||       ||  || ||    ||  ||  ||==    --  ||==      **  ||        ***   --  ||==    ||  ||==    ||      ||      ||  || \\||  ||
   \\\\      ||  || ||    ||  ||  ||          ||  \\\\      *  \\\\       *  *      ||      ||  ||      ||      ||      ||  ||  \\|  ||
      ===  ==  == ==    ==  ==  ==          ==   ==   **     ====   *  *      ==      ==  ==      ======  ======  ==  ==   ==  ======

@author: Rosanne C. van Deuren
"""

// Input sampleIDs for which candidate JSI variants have been merged.
input_sampleIDs_ch = Channel
  .fromPath( params.input_sampleIDs_path + params.sampleIDs_file )
  .splitCsv()
  .map{ line ->
    def sampleID = "${line[0]}";
  }

/*
 * >I.bedtoolsCoverage<
 *
 * Process bedtoolsCoverage computes the bedtools coverage depth of each 
 * of the PCR-replicate bamfiles containing CHMIP sequencing data, 
 * based on the targetCallingFile specified in the nextflow.config.
 *
 */

// Process bedtoolsCoverage.
process bedtoolsCoverage {

  tag "${sampleID}"

  input:
  val(sampleID) from input_sampleIDs_ch

  output:
  set val(sampleID), file("${sampleID}-PCR1.coverages.txt"), file("${sampleID}-PCR2.coverages.txt") into rawCoverages_ch

  """
  #!/bin/bash

  bedtools coverage -a ${params.reference_path}${params.reference_targetCallingFile} -b ${params.input_bams_path}${sampleID}-PCR1${params.bams_suffix} -d > "${sampleID}-PCR1.coverages.txt"
  bedtools coverage -a ${params.reference_path}${params.reference_targetCallingFile} -b ${params.input_bams_path}${sampleID}-PCR2${params.bams_suffix} -d > "${sampleID}-PCR2.coverages.txt"

  """
}

/*
 * >II.parseCoverage<
 *
 * Process parseCoverage parses and sorts the bedtools coverage depth
 * by running PythonScript.CHMIP-RsCh-PIPELINE.parseCoverages.py. 
 *
 */
 
// Process parseCoverage.
process parseCoverages {

  publishDir params.output_path + params.output_coverages, mode: "copy", overwrite: true
  tag "${sampleID}"

  input:
  set val(sampleID), file(rawCoverage_PCR1), file(rawCoverage_PCR2) from rawCoverages_ch

  output:
  set val(sampleID), file("${rawCoverage_PCR1.baseName}.parsed.txt"), file("${rawCoverage_PCR2.baseName}.parsed.txt") into parsedCoverages_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_parseCoverages} ${rawCoverage_PCR1} "${rawCoverage_PCR1.baseName}.tmpparsed.txt" "${rawCoverage_PCR1.baseName}.parsed.txt"
  python ${params.script_path}${params.script_parseCoverages} ${rawCoverage_PCR2} "${rawCoverage_PCR2.baseName}.tmpparsed.txt" "${rawCoverage_PCR2.baseName}.parsed.txt"

  """
}

/*
 * >III.averageCoverages<
 *
 * Process averageCoverages computes the average coverage depth over
 * the entire CHMIP-panel for PCR-replicates separately.
 * The average coverage is published in a summary file, and at the same
 * time used as input for the next process coverageFilter.
 *
 */
 
// Process averageCoverages.
process averageCoverages {

  tag "${sampleID}"

  input:
  set val(sampleID), file(parsedCoverage_PCR1), file(parsedCoverage_PCR2) from parsedCoverages_ch

  output:
  tuple file("${parsedCoverage_PCR1.baseName}.average.txt"), file("${parsedCoverage_PCR2.baseName}.average.txt") into averageCoverages_ch

  """
  #!/bin/bash

  echo "${parsedCoverage_PCR1.baseName}"
  python ${params.script_path}${params.script_averageCoverages} ${parsedCoverage_PCR1} "${parsedCoverage_PCR1.baseName}" "${parsedCoverage_PCR1.baseName}.average.txt"

  echo "${parsedCoverage_PCR2.baseName}"
  python ${params.script_path}${params.script_averageCoverages} ${parsedCoverage_PCR2} "${parsedCoverage_PCR2.baseName}" "${parsedCoverage_PCR2.baseName}.average.txt"

  """
}

// Combine average coverages into one file for output, but also duplicate for next process coverageFilter.
averageCoverages_ch
  .flatten()
  .collectFile(name: params.cohort + "_CHMIP-RsCh-PIPELINE_averageCoverages.txt", newLine: false, keepHeader: true, storeDir: params.output_path + params.output_coverages)
  .set { averageCoverages_coverageFilter_ch }

/*
 * >IV.coverageFilter<
 *
 * Process coverageFilter determines which samples pass the coverage-threshold
 * as specified by the coverage_threshold parameter in the nextflow.config.
 * By running PythonScript.CHMIP-RsCh-PIPELINE.coverageFilter.py it produces
 * a file containing only the sampleIDs that should be included in filtering
 * JSI candidate variant calls. The output channel is duplicated, as this file
 * is also used to determine the threshold in process excludeRunArtifacts.
 *
 */
 
// Process coverageFilter.
process coverageFilter {

  publishDir params.output_path + params.output_coverages, mode: "copy", overwrite: true

  input:
  file(averageCoverages) from averageCoverages_coverageFilter_ch

  output:
  file("${params.cohort}_samples_included.txt") into coveragePass_ch
  
  """
  #!/bin/bash

  python ${params.script_path}${params.script_coverageFilter} ${averageCoverages} ${params.coverage_threshold} "${params.cohort}_samples_included.txt"

  """
}

// Duplicate included samples channel for checkDuplicateCalls and excludeRunArtifacts.
coveragePass_ch.into{ coveragePass_checkDuplicates_ch; coveragePass_runArtifacts_ch }

// Input JSImerged.
input_JSImerged_ch = Channel
  .fromPath( params.input_JSI_path + params.JSImerged_file )
  .map{ file -> tuple(file.simpleName, file) }

/*
 * >V.checkDuplicateCalls<
 *
 * Process checkDuplicateCalls processes all candidate JSI variant calls, by running
 * PythonScript.CHMIP-RsCh-PIPELINE.extractDuplicateCalls.py.
 * 1) the entire file is recoded to parseable format with candidate 
 *    variant calls of both PCRs collapsed into one line.
 * 2) candidate variant calls of sampleIDs that do not pass the coverage
 *    threshold are excluded.
 * 3) all PTPN11 candidate variant calls are excluded based on likely homology.
 * 4) all candidate variant calls with an N-allele in the reference are excluded.
 * 5) candidate variant calls identified only in one of the two technical
 *    replicates (PCRs) are excluded.
 * 6) the output channel comprises candidate variant calls identified in both 
 *    technical replicates.
 *
 */
 
// Process checkDuplicateCalls.
process checkDuplicateCalls {

  publishDir params.output_path + params.output_JSI_excl, mode: "copy", overwrite: true, pattern: "*_EXCL_*"
  publishDir params.output_path + params.output_JSI_intermediate, mode: "copy", overwrite: true, pattern: "*_duplicates.txt"
  tag "${cohort}"

  input:
  tuple val(cohort), file(JSImerged) from input_JSImerged_ch
  file(inclSamples) from coveragePass_checkDuplicates_ch

  output:
  set val(cohort), file("${cohort}_recoded.txt"), file("${cohort}_recoded_EXCL_coverage.txt"), file("${cohort}_recoded_EXCL_PTPN11.txt"), file("${cohort}_recoded_EXCL_Nalt.txt"), file("${cohort}_recoded_EXCL_singletons.txt") into output_JSImerged_duplicates_ch
  tuple val(cohort), file("${cohort}_recoded_duplicates.txt") into JSImerged_duplicates_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_checkDuplicateCalls} ${JSImerged} "${cohort}_recoded.txt" "${cohort}_recoded_EXCL_coverage.txt" "${cohort}_recoded_EXCL_PTPN11.txt" "${cohort}_recoded_EXCL_Nalt.txt" "${cohort}_recoded_EXCL_singletons.txt" "${cohort}_recoded_duplicates.txt" ${inclSamples}
  """
}

/*
 * >VI.selectCodingNonsynSomatic<
 *
 * Process selectCodingNonsynSomatic processes duplicate candidate JSI variant calls,
 * by running PythonScript.CHMIP-RsCh-PIPELINE.selectCodingNonsynSomatic.py.
 * 1) non-coding candidate variant calls based on p.HGVS are excluded.
 * 2) synonymous candidate variant calls based on p.HGVS are excluded.
 * 3) germline candidate variant calls based on JSI VAF are excluded.
 * 4) the output channel comprises duplicate, coding, non-synonymous, somatic candidate variant calls.
 *
 */

// Process selectCodingNonsynSomatic.
process selectCodingNonsynSomatic {

  publishDir params.output_path + params.output_JSI_excl, mode: "copy", overwrite: true, pattern: "*_EXCL_*"
  publishDir params.output_path + params.output_JSI_intermediate, mode: "copy", overwrite: true, pattern: "*_codingNonsynonymousSomatic.txt"
  tag "${cohort}"

  input:
  tuple val(cohort), file(JSImerged_recoded_duplicates) from JSImerged_duplicates_ch

  output:
  set val(cohort), file("${cohort}_recoded_duplicates_EXCL_noncoding.txt"), file("${cohort}_recoded_duplicates_coding_EXCL_synonymous.txt"), file("${cohort}_recoded_duplicates_codingNonsynonymous_EXCL_germline.txt") into output_JSImerged_duplicates_codingNonsynSomatic_ch
  tuple val(cohort), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic.txt") into JSImerged_duplicates_codingNonsynSomatic_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_selectCodingNonsynSomatic} ${JSImerged_recoded_duplicates} "${cohort}_recoded_duplicates_EXCL_noncoding.txt" "${cohort}_recoded_duplicates_coding_EXCL_synonymous.txt" "${cohort}_recoded_duplicates_codingNonsynonymous_EXCL_germline.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic.txt"
  """

}

/*
 * >VII.excludeRunArtifacts<
 *
 * Process excludeRunArtifacts processes duplicate, coding, non-synonymous, somatic 
 * candidate JSI variant calls.
 * Depending on the nextflow.config parameter "artifacts" either all run-artifacts (common and
 * additional based on this particular run) are excluded (artifacts = "all", runs 
 * PythonScript.CHMIP-RsCh-PIPELINE.excludeRunArtifacts.py.), or only common run-artifacts are
 * excluded (artifacts = "common", runs PythonScript.CHMIP-RsCh-PIPELINE.excludeRunArtifacts_onlyCommon.py.
 * Note: Run-artifacts are overruled by known drivers as specified by the reference_knownDrivers parameter
 * in the nextflow.config.
 *
 */
 
// Process excludeRunArtifacts.
process excludeRunArtifacts {

  publishDir params.output_path + params.output_JSI_excl, mode: "copy", overwrite: true, pattern: "*_EXCL_*"
  publishDir params.output_path + params.output_JSI_intermediate, mode: "copy", overwrite: true, pattern: "*_clean.txt"
  tag "${cohort}"

  input:
  tuple val(cohort), file(JSImerged_recoded_duplicates_codingNonsynSomatic) from JSImerged_duplicates_codingNonsynSomatic_ch
  file(inclSamples) from coveragePass_runArtifacts_ch

  output:
  set val(cohort), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_EXCL_commonArtifacts.txt"), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_EXCL_addArtifacts.txt") into output_JSImerged_duplicates_codingNonsynSomatic_clean_ch
  tuple val(cohort), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean.txt") into JSImerged_duplicates_codingNonsynSomatic_clean_ch

  script:
  if( params.artifacts == "all" )
  """
  #!/bin/bash

  python ${params.script_path}${params.script_allRunArtifacts} ${JSImerged_recoded_duplicates_codingNonsynSomatic} ${params.reference_path}${params.reference_commonArtifacts} ${params.reference_path}${params.reference_knownDrivers} "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_excl_commonArtifacts.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_excl_addArtifacts.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean.txt" ${inclSamples}
  """
  else if( params.artifacts == "common" )
  """
  #!/bin/bash

  python ${params.script_path}${params.script_commonRunArtifacts} ${JSImerged_recoded_duplicates_codingNonsynSomatic} ${params.reference_path}${params.reference_commonArtifacts} ${params.reference_path}${params.reference_knownDrivers} "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_EXCL_commonArtifacts.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_EXCL_addArtifacts.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean.txt" ${inclSamples}
  """
  else
    error "Invalid argument for parameter artifacts in nextflow.config: ${params.artifacts}"
}

/*
 * >VIII.finetuneCandidateCalls<
 *
 * Process finetuneCandidateCalls processes duplicate, coding, non-synonymous, somatic 
 * candidate JSI variant calls for which run-artifacts have been excluded, by running
 * PythonScript.CHMIP-RsCh-PIPELINE.finetuneCandidateDalls.py. Candidate variant calls
 * are further finetuned on three aspects:
 * 1) inspection-flags: manual inspection provided green and red flags for candidate variant
 *    calls that should pass and not pass respectively.
 * 2) variant-statistics: for candidate variant calls identified in 4 or more individuals, an
 *    allele-count of 16 or higher is required for at least one of the candidate calls; for 
 *    candidate variant calls identified in less than 4 individuals, an allele-count of 24
 *    or higher is required for at least one of the candidate calls.
 * 3) complex-statistics: for candidate variant calls identified in the same individual, no
 *    more than three calls are allowed within 54bp.
 * Note: Aspect 1) overrules aspect 2) and 3).
 * The output consists of final candidate variant calls and the corresponding positions to
 * be submitted for mpileup in the next process.
 *
 */
 
// Process finetuneCandidateCalls.
process finetuneCandidateCalls {

  publishDir params.output_path + params.output_JSI_excl, mode: "copy", overwrite: true, pattern: "*_EXCL_*"
  publishDir params.output_path + params.output_JSI_final, mode: "copy", overwrite: true, pattern: "*_finetuned_inspectIGV.txt"
  publishDir params.output_path + params.output_JSI_intermediate, mode: "copy", overwrite: true, pattern: "*_finetuned.txt"
  publishDir params.output_path + params.output_mpileups, mode: "copy", overwrite: true, pattern: "*_mpileuppositions.txt"
  tag "${cohort}"

  input:
  tuple val(cohort), file(JSImerged_recoded_duplicates_codingNonsynSomatic_clean) from JSImerged_duplicates_codingNonsynSomatic_clean_ch

  output:
  set val(cohort), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_inspectionFlag.txt"), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_varstats.txt"), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_complexstats.txt"), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned_inspectIGV.txt") into output_JSImerged_duplicates_codingNonsynSomatic_clean_final_ch
  tuple val(cohort), file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned.txt") into JSImerged_final_forMpileup_ch
  file("${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned_mpileuppositions.txt") into mpileup_positions_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_finetuneCandidateCalls} ${JSImerged_recoded_duplicates_codingNonsynSomatic_clean} ${params.reference_path}${params.reference_inspectionFlags} "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_inspectionFlag.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_varstats.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_EXCL_complexstats.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned_mpileuppositions.txt" "${cohort}_recoded_duplicates_codingNonsynonymousSomatic_clean_finetuned_inspectIGV.txt"
  """

}

/*
 * Mpileup-positions stored in output-file are splitted so that samtools mpileup can be
 * run in parallel in the next process.
 *
 */
 
// Split mpileup positions into multiple streams.
mpileup_positions_ch
  .splitCsv(sep: "\t", header: true)
  .map{ line ->
    def sampleID = "${line.sampleID}";
    def chromosome = "${line.chromosome}";
    def bpStart = "${line.bpStart}";
    def bpEnd = "${line.bpEnd}";
    def position = "${line.chromosome}:${line.bpStart}-${line.bpEnd}";
    def ref = "${line.ref}";
    def alt = "${line.alt}";
    def gene = "${line.gene}";

    [sampleID, chromosome, bpStart, bpEnd, position, ref, alt, gene]
  }
  .set { input_mpileup_positions_ch }

/*
 * >IX.samtoolsMpileup<
 *
 * Process samtoolsMpileup computes an mpileup on both technical replicates (PCRs) for
 * each of the final candidate variant calls.
 *
 */
 
// Process samtoolsMpileup.
process samtoolsMpileup {

  tag "${sampleID}"

  input:
  set val(sampleID), val(chromosome), val(bpStart), val(bpEnd), val(position), val(ref), val(alt), val(gene) from input_mpileup_positions_ch

  output:
  set val(sampleID), file("${sampleID}-PCR1.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}.mpileup.txt"), file("${sampleID}-PCR2.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}.mpileup.txt") into output_mpileups_ch

  """
  #!/bin/bash

  echo "${sampleID}-PCR1.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}"
  samtools mpileup -r ${position} -d 10000 -q 30 -Q 20 -B -f ${params.reference_hg19} ${params.input_bams_path}${sampleID}-PCR1${params.bams_suffix} > "${sampleID}-PCR1.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}.mpileup.txt"

  echo "${sampleID}-PCR2.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}"
  samtools mpileup -r ${position} -d 10000 -q 30 -Q 20 -B -f ${params.reference_hg19} ${params.input_bams_path}${sampleID}-PCR2${params.bams_suffix} > "${sampleID}-PCR2.${chromosome}.${bpStart}.${bpEnd}.${ref}.${alt}.${gene}.mpileup.txt"

  """

}

/*
 * >X.preprocessMpileups<
 *
 * Process preprocessMpileups parses mpileup-sequences and saves information on
 * read-marker presence.
 *
 */
 
// Process preprocessMpileups.
process preprocessMpileups {

  publishDir params.output_path + params.output_mpileups, mode: "copy", overwrite: true, pattern: "*.prepped.txt"
  tag "${sampleID}"

  input:
  set val(sampleID), file(mpileupPCR1), file(mpileupPCR2) from output_mpileups_ch

  output:
  set val(sampleID), file("${mpileupPCR1.baseName}.prepped.txt"), file("${mpileupPCR2.baseName}.prepped.txt") into output_preppedMpileups_ch
  set val(sampleID), file("${mpileupPCR1.baseName}.readmarkers.txt"), file("${mpileupPCR2.baseName}.readmarkers.txt") into output_readmarkerMpileups_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_preprocessMpileups} ${mpileupPCR1} "${mpileupPCR1.baseName}.prepped.txt" "${mpileupPCR1.baseName}.readmarkers.txt"
  python ${params.script_path}${params.script_preprocessMpileups} ${mpileupPCR2} "${mpileupPCR2.baseName}.prepped.txt" "${mpileupPCR2.baseName}.readmarkers.txt"

  """
}

/*
 * >XI.parseMpileups<
 *
 * Process parseMpileups runs PerlScript.CHMIP-RsCh-PIPELINE.pileup2baseindel.no.strand.pl
 * which is a script from library pileup2base (https://github.com/riverlee/pileup2base.git)
 * that parses our preprocessed mpileup-sequences so that variant allele frequencies may be
 * calculated.
 *
 */
 
// Process parseMpileups.
process parseMpileups {

  tag "${sampleID}"

  input:
  set val(sampleID), file(preppedMpileupPCR1), file(preppedMpileupPCR2) from output_preppedMpileups_ch

  output:
  set val(sampleID), file("${preppedMpileupPCR1.baseName}.parsed.txt"), file("${preppedMpileupPCR2.baseName}.parsed.txt") into output_parsedMpileups_ch

  """
  #!/bin/bash

  perl ${params.script_path}${params.script_parseMpileups} ${preppedMpileupPCR1} 1 "${preppedMpileupPCR1.baseName}.parsed.txt"
  perl ${params.script_path}${params.script_parseMpileups} ${preppedMpileupPCR2} 1 "${preppedMpileupPCR2.baseName}.parsed.txt"

  """
}

// Combine parsed mpileups with readmarkers.
output_parsedMpileups_ch
  .join(output_readmarkerMpileups_ch)
  .map{ it ->
    def sampleID = file(it[0]).getName();
    def parsedMpileupPCR1 = file(it[1]);
    def parsedMpileupPCR2 = file(it[2]);
    def readmarkerInfoPCR1 = file(it[3]);
    def readmarkerInfoPCR2 = file(it[4]);

    [sampleID, parsedMpileupPCR1, readmarkerInfoPCR1, parsedMpileupPCR2, readmarkerInfoPCR2]
  }
  .set { parsedMpileups_readmarkerInfo_ch }

/*
 * >XII.concludeMpileups<
 *
 * Process concludeMpileups runs PythonScript.CHMIP-RsCh-PIPELINE.concludeMpileups.py to
 * calculate variant allele frequencies based on bamfiles. 
 *
 */
 
// Process concludeMpileups.
process concludeMpileups {

  tag "${sampleID}"

  input:
  set val(sampleID), file(parsedMpileupPCR1), file(readmarkerInfoPCR1), file(parsedMpileupPCR2), file(readmarkerInfoPCR2) from parsedMpileups_readmarkerInfo_ch

  output:
  tuple file("${parsedMpileupPCR1.baseName}.concluded.txt"), file("${parsedMpileupPCR2.baseName}.concluded.txt") into concludedMpileups_ch

  """
  #!/bin/bash

  echo "${parsedMpileupPCR1.baseName}"
  python ${params.script_path}${params.script_concludeMpileupVAFs} ${parsedMpileupPCR1} ${readmarkerInfoPCR1} "${parsedMpileupPCR1.baseName}" "${parsedMpileupPCR1.baseName}.concluded.txt"

  echo "${parsedMpileupPCR2.baseName}"
  python ${params.script_path}${params.script_concludeMpileupVAFs} ${parsedMpileupPCR2} ${readmarkerInfoPCR2} "${parsedMpileupPCR2.baseName}" "${parsedMpileupPCR2.baseName}.concluded.txt"

  """
}

// Save concluded mpileups to file and create new channel of that file for final process.
concludedMpileups_ch
  .flatten()
  .collectFile(name: params.cohort + "CHMIP-RsCh-PIPELINE_concludedMpileups.txt", newLine: false, keepHeader: true, storeDir: params.output_path + params.output_mpileups )
  .set { concludedMpileups_merged_ch }

/*
 * >XIII.integrateMpileupVAFs<
 *
 * Process integrateMpileupVAFs runs PythonScript.CHMIP-RsCh-PIPELINE.integrateMpileupVAFfinalJSIvariants.py
 * which combines the concluded mpileup information with final JSI variant calls.
 * Variants for which mpileup did not provide VAFs are recommended for inspection in IGV.
 *
 */
 
// Process integrateMpileupVAFs.
process integrateMpileupVAFs {

  publishDir params.output_path + params.output_JSI_final, mode: "copy", overwrite: true
  tag "${params.cohort}"

  input:
  file(concludedMpileups_merged) from concludedMpileups_merged_ch
  tuple val(cohort), file(finalJSIvariants_merged) from JSImerged_final_forMpileup_ch

  output:
  set file("${params.cohort}_CHMIP-RsCh-PIPELINE_final.txt"), file("${finalJSIvariants_merged.baseName}_noMpileupVAF_inspectIGV.txt") into finalOutput_ch

  """
  #!/bin/bash

  python ${params.script_path}${params.script_integrateMpileupVAFs} ${concludedMpileups_merged} ${finalJSIvariants_merged} "${params.cohort}_CHMIP-RsCh-PIPELINE_final.txt" "${finalJSIvariants_merged.baseName}_noMpileupVAF_inspectIGV.txt"

  """
}