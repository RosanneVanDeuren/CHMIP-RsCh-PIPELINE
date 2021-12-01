## Python script to extract coding, non-synonymous, somatic variants from duplicates - nextflow executable.

# Import modules.
import numpy as np
import pandas as pd
import datetime as dt
import sys

# Define parameters.
INPUT_JSIMERGED_DUPLICATES = str(sys.argv[1])
OUTPUT_JSIMERGED_EXCL_NONCODING = str(sys.argv[2])
OUTPUT_JSIMERGED_EXCL_SYNONYMOUS = str(sys.argv[3])
OUTPUT_JSIMERGED_GERMLINE = str(sys.argv[4])
OUTPUT_JSIMERGED_CODINGNONSYNSOMATIC = str(sys.argv[5])

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
# Define function that recodes VAF.
def recodeVAF(JSIaltVAF):
    recJSIaltVAF = float(JSIaltVAF[0:len(JSIaltVAF)-1])
    return(recJSIaltVAF)

print(dt.datetime.now())
print("Start.")

# Import variants.
JSIduplicates_df = pd.read_csv(INPUT_JSIMERGED_DUPLICATES, header = 0, sep = "\t")

# Open new noncoding, synonymous, germline, and codingNonsynSomatic files to write to.
JSIduplicates_noncoding_open = open(OUTPUT_JSIMERGED_EXCL_NONCODING, "w")
JSIduplicates_noncoding_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIduplicates_synonymous_open = open(OUTPUT_JSIMERGED_EXCL_SYNONYMOUS, "w")
JSIduplicates_synonymous_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIduplicates_germline_open = open(OUTPUT_JSIMERGED_GERMLINE, "w")
JSIduplicates_germline_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIduplicates_codingNonsynSomatic_open = open(OUTPUT_JSIMERGED_CODINGNONSYNSOMATIC, "w")
JSIduplicates_codingNonsynSomatic_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))

# Iterate over duplicate variants and store variants in appropriate files.
for index, row in JSIduplicates_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2 = getJSIduplicatesInfo(row)
    if pHGVS is np.nan:
        JSIduplicates_noncoding_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    elif "=" in pHGVS:
        JSIduplicates_synonymous_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    elif (recodeVAF(JSIaltVAFPCR1) >= 40) or (recodeVAF(JSIaltVAFPCR2) >= 40):
        JSIduplicates_germline_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    else:
        JSIduplicates_codingNonsynSomatic_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))

# Close opened files.
JSIduplicates_noncoding_open.close()
JSIduplicates_synonymous_open.close()
JSIduplicates_germline_open.close()
JSIduplicates_codingNonsynSomatic_open.close()

print("Done.")
print(dt.datetime.now())
