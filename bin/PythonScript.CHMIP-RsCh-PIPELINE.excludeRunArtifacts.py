## Python script to exclude all (common and current run) artifacts - Nextflow executable.

# Import modules.
import numpy as np
import pandas as pd
import datetime as dt
import sys

# Define parameters.
INPUT_JSIMERGED_CODINGNONSYNSOMATIC = str(sys.argv[1])
REFERENCE_COMMONARTIFACTS = str(sys.argv[2])
REFERENCE_KNOWNDRIVERS = str(sys.argv[3])
OUTPUT_JSIMERGED_EXCL_COMMONARTIFACTS = str(sys.argv[4])
OUTPUT_JSIMERGED_EXCL_ADDARTIFACTS = str(sys.argv[5])
OUTPUT_JSIMERGED_CLEAN = str(sys.argv[6])
INCL_SAMPLES = str(sys.argv[7])

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

print(dt.datetime.now())
print("Start.")

# Import variants.
JSIcodingNonsynSomatic_df = pd.read_csv(INPUT_JSIMERGED_CODINGNONSYNSOMATIC, header = 0, sep = "\t")

# Included sampleIDs list for threshold. NB. will be supplied on command line in process.
sampleIDs_included_list = set((np.loadtxt(INCL_SAMPLES, dtype='string')).tolist())
THRESHOLD = np.floor(0.05*(len(sampleIDs_included_list)))

# Import common run artifacts (overallTHRESHOLD-file based on SOSruns).
commonArtifacts_df = pd.read_csv(REFERENCE_COMMONARTIFACTS, header = None, sep = "\t")
commonArtifacts_list = list(commonArtifacts_df.iloc[:, 0])

# Import known drivers.
knownDrivers_df = pd.read_csv(REFERENCE_KNOWNDRIVERS, header = None, sep = "\t")
knownDrivers_list = list(knownDrivers_df.iloc[:, 0])

# Open new common artifacts, additional artifacts and cleaned files to write to.
JSIcodingNonsynSomatic_commonArtifacts_open = open(OUTPUT_JSIMERGED_EXCL_COMMONARTIFACTS, "w")
JSIcodingNonsynSomatic_commonArtifacts_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIcodingNonsynSomatic_addArtifacts_open = open(OUTPUT_JSIMERGED_EXCL_ADDARTIFACTS, "w")
JSIcodingNonsynSomatic_addArtifacts_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSIcodingNonsynSomatic_clean_open = open(OUTPUT_JSIMERGED_CLEAN, "w")
JSIcodingNonsynSomatic_clean_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))

# Iterate over table and count variants.
variantCount_dict = {}
for index, row in JSIcodingNonsynSomatic_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2 = getJSIduplicatesInfo(row)
    variant_key = "{0}|{1}|{2}|{3}|{4}|{5}".format(chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS)
    if variant_key not in variantCount_dict.keys():
        variantCount_dict[variant_key] = 1
    else:
        variantCount_dict[variant_key] +=1

# Iterate over variantCount_dict and store additonal artifacts that are not known drivers in list.
addArtifacts_list = []
for key, value in variantCount_dict.items():
    if value >= THRESHOLD:
        if key not in knownDrivers_list:
            addArtifacts_list.append(key)

# Iterate over table and exclude common artifacts, and then additional artifacts.
for index, row in JSIcodingNonsynSomatic_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2 = getJSIduplicatesInfo(row)
    variant_key = "{0}|{1}|{2}|{3}|{4}|{5}".format(chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS)
    if variant_key in commonArtifacts_list:
        JSIcodingNonsynSomatic_commonArtifacts_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    elif variant_key in addArtifacts_list:
        JSIcodingNonsynSomatic_addArtifacts_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))
    else:
        JSIcodingNonsynSomatic_clean_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAFPCR1, JSIaltCountPCR1, JSIaltVAFPCR2, JSIaltCountPCR2))

# Close opened files.
JSIcodingNonsynSomatic_commonArtifacts_open.close()
JSIcodingNonsynSomatic_addArtifacts_open.close()
JSIcodingNonsynSomatic_clean_open.close()

print("Done.")
print(dt.datetime.now())
