## Python script to finetune candidate calls based on inspection-flags, variant- and complex statistics - Nextflow executable.

# Import modules.
import pandas as pd
import numpy as np
import datetime as dt
import sys

# Define parameters.
INPUT_JSIMERGED_CLEAN = str(sys.argv[1])
REFERENCE_INSPECTIONFLAGS = str(sys.argv[2])
OUTPUT_JSIMERGED_EXCL_FLAG = str(sys.argv[3])
OUTPUT_JSIMERGED_EXCL_VARSTATS = str(sys.argv[4])
OUTPUT_JSIMERGED_EXCL_COMPLEXSTATS = str(sys.argv[5])
OUTPUT_JSIMERGED_FINAL = str(sys.argv[6])
OUTPUT_MPILEUPPOSITIONS = str(sys.argv[7])
OUTPUT_IGV = str(sys.argv[8])

# Define function that extracts information from JSIduplicates.
def getJSIduplicatesInfo(row):
    sampleID = row[0]
    chromosome = row[1]
    bpStart = row[2]
    bpEnd = row[3]
    gene = row[4]
    cHGVS = row[5]
    pHGVS = row[6]
    BASEchange = row[7]
    AAchange = row[8]
    JSIaltVAFPCR1 = row[9]
    JSIaltCountPCR1 = row[10]
    JSIaltVAFPCR2 = row[11]
    JSIaltCountPCR2 = row[12]
    return(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2)
# Define function that extracts information for mpileups.
def getMpileupInformation(chromosome, bpStart, BASEchange):
    if ">" in BASEchange:
        ref, alt = BASEchange.split(">")
        if ref == "*": # insertion.
            CHROMOSOME = chromosome
            BPstart = bpStart
            BPend = bpStart
            REF = "NA"
            ALT = ("ins" + alt)
        elif alt == "*": # deletion.
            CHROMOSOME = chromosome
            BPstart = (bpStart-1)
            BPend = (bpStart-1)
            REF = ref
            ALT = ("del" + ref)
        elif len(ref) == 1 and len(alt) == 1: # substitution.
            CHROMOSOME = chromosome
            BPstart = bpStart
            BPend = bpStart
            REF = ref
            ALT = alt
        else:
            CHROMOSOME = "manual"
            BPstart = "manual"
            BPend = "manual"
            REF = "manual"
            ALT = "manual"
    else:
        CHROMOSOME = "manual"
        BPstart = "manual"
        BPend = "manual"
        REF = "manual"
        ALT = "manual"
    return(CHROMOSOME, BPstart, BPend, REF, ALT)

print(dt.datetime.now())
print("Start.")

# Import cleaned variants.
JSIcleaned_df = pd.read_csv(INPUT_JSIMERGED_CLEAN, header = 0, sep = "\t")

# Import referenceFile inspectionFlags.
inspectionflags_df = pd.read_csv(REFERENCE_INSPECTIONFLAGS, header = 0, sep = "\t")

# Iterate over inspectionflags_df and store in dictionary.
inspectionFlags_dict = {}
for index, row in inspectionflags_df.iterrows():
    variant_key = inspectionflags_df.loc[index, 'variant_key']
    flag = inspectionflags_df.loc[index, 'flag']
    inspectionFlags_dict[variant_key] = flag

# Iterate over cleaned variants, and count and keep information of variants for variant- and complex-statistics.
variantSampleID_dict = {}
sampleIDvariants_dict = {}
variantAC_dict = {}
for index, row in JSIcleaned_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2 = getJSIduplicatesInfo(row)
    variant_key = "{0}|{1}|{2}|{3}|{4}|{5}".format(chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS)
    variantSampleID_dict.setdefault(variant_key, []).append(sampleID)
    sampleIDvariants_dict.setdefault(sampleID, []).append(variant_key)
    AC_PCR1 = int(JSIaltCountPCR1[1:(len(JSIaltCountPCR1)-1)])
    AC_PCR2 = int(JSIaltCountPCR2[1:(len(JSIaltCountPCR2)-1)])
    variantAC_dict.setdefault(variant_key, []).append(AC_PCR1)
    variantAC_dict.setdefault(variant_key, []).append(AC_PCR2)

# Create variantStat_list.
variantStat_EXCL_list = []
for key, value in variantSampleID_dict.items():
    maxAC = max(variantAC_dict[key])
    if len(value) >= 4: # Ncalls >= 4 --> highest AC needs to be >= 16.
        if maxAC < 16:
            variantStat_EXCL_list.append(key) # NB. not sample-specific, all samples with this variant pass.
    else: # Ncalls < 4 --> highest AC needs to be >= 24.
        if maxAC < 24:
            variantStat_EXCL_list.append(key) # NB. not sample-specific, all samples with this variant pass.

# Create complexStat_list.
complexStat_EXCL_list = []
for key, value in sampleIDvariants_dict.items():
    if len(value) > 3:
        chromosomeBPstart_dict = {}
        chromosomeVariant_dict = {}
        for variant in value:
            chromosome, bpStart, bpStop, gene, cHGVS, pHGVS = variant.split('|')
            chromosomeBPstart_dict.setdefault(chromosome, []).append(bpStart)
            chromosomeVariant_dict.setdefault(chromosome, []).append(variant)
        # check proximity of variants on the same chromosome.
        for CHROM, BP in chromosomeBPstart_dict.items():
            if len(BP) > 3:
                sortedBP = BP
                sortedBP.sort()
                if len(BP) == 4:
                    if ((int(sortedBP[-1]) - int(sortedBP[0])) < 54): # all variants should receive "NO_PASS" flag.
                        variant_keys = chromosomeVariant_dict[CHROM]
                        for variant_key in variant_keys:
                            complexKey = "{0}:{1}".format(key, variant_key)
                            complexStat_EXCL_list.append(complexKey)
                else:
                    maxItem = len(sortedBP)-1
                    for ix in range(0, len(sortedBP)):
                        if (ix+3) <= maxItem:
                            if (int(sortedBP[ix+3]) - int(sortedBP[ix])) < 54:
                                ixVar1, ixVar2, ixVar3, ixVar4 = BP.index(sortedBP[ix]), BP.index(sortedBP[ix+1]), BP.index(sortedBP[ix+2]), BP.index(sortedBP[ix+3])
                                complexStat_EXCL_list.append(complexKey).append(key + ":" + chromosomeVariant_dict[CHROM][ixVar1])
                                complexStat_EXCL_list.append(complexKey).append(key + ":" + chromosomeVariant_dict[CHROM][ixVar2])
                                complexStat_EXCL_list.append(complexKey).append(key + ":" + chromosomeVariant_dict[CHROM][ixVar3])
                                complexStat_EXCL_list.append(complexKey).append(key + ":" + chromosomeVariant_dict[CHROM][ixVar4])

# Create new JSI-files to write to.
JSIcleaned_final_open = open(OUTPUT_JSIMERGED_FINAL, "w")
JSIcleaned_final_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIcleaned_exclude_redflag_open = open(OUTPUT_JSIMERGED_EXCL_FLAG, "w")
JSIcleaned_exclude_redflag_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIcleaned_exclude_varstat_open = open(OUTPUT_JSIMERGED_EXCL_VARSTATS, "w")
JSIcleaned_exclude_varstat_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIcleaned_exclude_complexstat_open = open(OUTPUT_JSIMERGED_EXCL_COMPLEXSTATS, "w")
JSIcleaned_exclude_complexstat_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
# Create mpileup-positions file to write to.
mpileup_positions_open = open(OUTPUT_MPILEUPPOSITIONS, "w")
mpileup_positions_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "ref", "alt", "gene"))
# Create file for which no mpileup possible, so inspection IGV needed.
noMpileupIGV_open = open(OUTPUT_IGV, "w")
noMpileupIGV_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))


# Iterate over cleaned variants and exclude or include variants based on inspection-flag, variant-stats or complex-stats.
# For included variants extract the positions for mpileup at the same time.
for index, row in JSIcleaned_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2 = getJSIduplicatesInfo(row)
    variant_key = "{0}|{1}|{2}|{3}|{4}|{5}".format(chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS)
    sampleIDvariant_key = "{0}:{1}|{2}|{3}|{4}|{5}|{6}".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS)
    if variant_key in inspectionFlags_dict.keys():
        if inspectionFlags_dict[variant_key] == "green":
            CHROMOSOME, BPstart, BPend, REF, ALT = getMpileupInformation(chromosome, bpStart, BASEchange)
            if CHROMOSOME == "manual":
                noMpileupIGV_open.write(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2)
            else:
                mpileup_positions_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(sampleID, CHROMOSOME, BPstart, BPend, REF, ALT, gene))
                JSIcleaned_final_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
        else:
            JSIcleaned_exclude_redflag_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    elif variant_key in variantStat_EXCL_list:
        JSIcleaned_exclude_varstat_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    elif sampleIDvariant_key in complexStat_EXCL_list:
        JSIcleaned_exclude_varstat_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    else:
        CHROMOSOME, BPstart, BPend, REF, ALT = getMpileupInformation(chromosome, bpStart, BASEchange)
        if CHROMOSOME == "manual":
            noMpileupIGV_open.write(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2)
        else:
            mpileup_positions_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(sampleID, CHROMOSOME, BPstart, BPend, REF, ALT, gene))
            JSIcleaned_final_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))

# Close files.
JSIcleaned_final_open.close()
JSIcleaned_exclude_redflag_open.close()
JSIcleaned_exclude_varstat_open.close()
JSIcleaned_exclude_complexstat_open.close()
mpileup_positions_open.close()
noMpileupIGV_open.close()

print("Done.")
print(dt.datetime.now())
