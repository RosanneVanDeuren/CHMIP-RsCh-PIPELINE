## Python script to conclude mpileup VAFs of final JSI variants - Nextflow executable.

# Import modules.
import pandas as pd
import numpy as np
import sys

# Define parameters.
INPUT_PILEUP_PARSED = str(sys.argv[1])
INPUT_PILEUP_READMARKERS = str(sys.argv[2])
INPUT_PILEUP_INFO = str(sys.argv[3])
OUTPUT_PILEUP_CONCLUDED = str(sys.argv[4])

# Define functions.
def getMpileupInfo(baseName):
    sampleID, rest = baseName.split("-")
    pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, baseNameRest1, baseNameRest2, baseNameRest3 = rest.split(".")
    return(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene)

# List with genes on the minus strand.
minus_strand_genes_list = ["GNB1", "NRAS", "DNMT3A", "SF3B1", "BRAF", "FGFR2", "HRAS", "KRAS", "IDH2", "TP53", "STAT3", "SRSF2", "U2AF1"]
# Alternative allele dictionary for genes on minus strand.
alternative_allele_dict = {'A':'T','C':'G','G':'C','T':'A'}

# Open output pileup concluded file.
concludedMpileup_open = open(OUTPUT_PILEUP_CONCLUDED, "w")
concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "pcrID", "chromosome", "bpStart", "bpEnd", "ref", "alt", "gene", "mpileupConcluded", "mpileupVAF", "mpileupALTcount", "mpileupTOTALcount", "mpileupNote"))

# Open parsed mpileup, compute total read count and calculate VAF using pileup info.
parsedMpileup_df = pd.read_csv(INPUT_PILEUP_PARSED, sep = "\t", header = 0)
parsedMpileup_df["total"] = parsedMpileup_df["A"] + parsedMpileup_df["C"] + parsedMpileup_df["G"] + parsedMpileup_df["T"]
sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene = getMpileupInfo(INPUT_PILEUP_INFO)
# Check if there is an mpileup.
if not parsedMpileup_df.size > 0:
    concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, "FP", "NA", "NA", "NA", "NO_MPILEUP_CHECK_IGV"))
else:
    # Check for total read count.
    if parsedMpileup_df.loc[0, "total"] < 100:
        concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, "FP", "NA", "NA", "NA", "NO_PASS_TOTAL_COUNT"))
    else:
        # Check for read marker.
        try:
            readMarker_df = pd.read_csv(INPUT_PILEUP_READMARKERS, sep = "\t", header = 0)
            mpileupNote = readMarker_df.loc[0, "mpileupSeqMarker"]
        except:
            mpileupNote = "NA"
        # Calculate VAF.
        if "ins" in alt:
            if pd.isnull(parsedMpileup_df.loc[0, "Insertion"]):
                ALTcount = 0
                TOTALcount = parsedMpileup_df.loc[0, "total"]
                VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
                conclusion = "FP"
            else:
                if gene in minus_strand_genes_list:
                    ALT = alt[3:]
                    mpileupALT = ""
                    for base in ALT:
                        baseREV = alternative_allele_dict[base]
                        mpileupALT = str(baseREV) + mpileupALT
                    INSopts = (parsedMpileup_df.loc[0, 'Insertion']).split('|')
                    for INSopt in INSopts:
                        INScount, INSallele = INSopt.split(':')
                        if INSallele == mpileupALT:
                            ALTcount = INScount
                            TOTALcount = parsedMpileup_df.loc[0, "total"]
                            VAF = round((float(INScount) / TOTALcount) * 100, 2)
                else:
                    mpileupALT = alt[3:]
                    INSopts = (parsedMpileup_df.loc[0, 'Insertion']).split('|')
                    for INSopt in INSopts:
                        INScount, INSallele = INSopt.split(':')
                        if INSallele == mpileupALT:
                            ALTcount = INScount
                            TOTALcount = parsedMpileup_df.loc[0, "total"]
                            VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
                if not VAF == 0.0:
                    conclusion = "TP"
                else:
                    conclusion = "FP"
            concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, conclusion, VAF, ALTcount, TOTALcount, mpileupNote))
        elif "del" in alt:
            if pd.isnull(parsedMpileup_df.loc[0, "Deletion"]):
                ALTcount = 0
                TOTALcount = parsedMpileup_df.loc[0, "total"]
                VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
                conclusion = "FP"
            else:
                if gene in minus_strand_genes_list:
                    ALT = alt[3:]
                    mpileupALT = ""
                    for base in ALT:
                        baseREV = alternative_allele_dict[base]
                        mpileupALT = str(baseREV) + mpileupALT
                    DELopts = (parsedMpileup_df.loc[0, 'Deletion']).split('|')
                    for DELopt in DELopts:
                        DELcount, DELallele = DELopt.split(':')
                        if DELallele == mpileupALT:
                            ALTcount = DELcount
                            TOTALcount = parsedMpileup_df.loc[0, "total"]
                            VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
                else:
                    mpileupALT = alt[3:]
                    INSopts = (parsedMpileup_df.loc[0, 'Insertion']).split('|')
                    for INSopt in INSopts:
                        INScount, INSallele = INSopt.split(':')
                        if INSallele == mpileupALT:
                            ALTcount = INScount
                            TOTALcount = parsedMpileup_df.loc[0, "total"]
                            VAF = round((float(INScount) / TOTALcount) * 100, 2)
                if not VAF == 0.0:
                    conclusion = "TP"
                else:
                    conclusion = "FP"
            concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, conclusion, VAF, ALTcount, TOTALcount, mpileupNote))
        else:
            if gene in minus_strand_genes_list:
                mpileupALT = alternative_allele_dict[alt]
                ALTcount = int(parsedMpileup_df.loc[0, mpileupALT])
                TOTALcount = int(parsedMpileup_df.loc[0, "total"])
                VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
            else:
                ALTcount = int(parsedMpileup_df.loc[0, alt])
                TOTALcount = int(parsedMpileup_df.loc[0, "total"])
                VAF = round((float(ALTcount) / TOTALcount) * 100, 2)
            if not VAF == 0.0:
                conclusion = "TP"
            else:
                conclusion = "FP"
            # Write conclusion mpileup
            concludedMpileup_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, conclusion, VAF, ALTcount, TOTALcount, mpileupNote))

# Close conclusion mpileup file.
concludedMpileup_open.close()
