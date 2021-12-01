## Python script to extract duplicate variant calls from JSImerged - nextflow executable.

# Import modules.
import pandas as pd
import numpy as np
import datetime as dt
import sys

# Define parameters.
INPUT_JSIMERGED = str(sys.argv[1])
OUTPUT_JSIMERGED_RECODED = str(sys.argv[2])
OUTPUT_JSIMERGED_RECODED_EXCL_COVERAGE = str(sys.argv[3])
OUTPUT_JSIMERGED_RECODED_EXCL_PTPN11 = str(sys.argv[4])
OUTPUT_JSIMERGED_RECODED_EXCL_NALT = str(sys.argv[5])
OUTPUT_JSIMERGED_RECODED_EXCL_SINGLETONS = str(sys.argv[6])
OUTPUT_JSIMERGED_RECODED_DUPLICATES = str(sys.argv[7])
INCL_SAMPLES = str(sys.argv[8])

# Define function that recodes JSI information.
def recodeJSIinfo(row):
    sampleID = row[0].split("-")[0]
    pcrID = row[0].split("-")[1][0:4]
    posInfo = row[2].split("[")[1]
    if "/" in posInfo:
        chromosome = posInfo.split(":")[0]
        if "_" in posInfo.split(":")[1]:
            bpStart = posInfo.split(" / ")[1].split("_")[0].replace(" (hg19)]", "")
            bpEnd = posInfo.split(" / ")[1].split("_")[1].replace(" (hg19)]", "")
        else:
            bpStart = posInfo.split(" / ")[1].replace(" (hg19)]", "")
            bpEnd = posInfo.split(" / ")[1].replace(" (hg19)]", "")
    else:
        chromosome = posInfo.split(":")[0]
        if "_" in posInfo.split(":")[1]:
            bpStart = posInfo.split(":")[1].split("_")[0].replace("g.", "").replace(" (hg19)]", "")
            bpEnd = posInfo.split(":")[1].split("_")[1].replace("g.", "").replace(" (hg19)]", "")
        else:
            bpStart = posInfo.split(":")[1].replace("g.", "").replace(" (hg19)]", "")
            bpEnd = posInfo.split(":")[1].replace("g.", "").replace(" (hg19)]", "")
    gene = row[1]
    cHGVS = row[7]
    pHGVS = row[8]
    if row[5] is np.nan:
        BASEchange = row[5]
    elif "->" in row[5]:
        BASEref = row[5].split(" -> ")[0]
        BASEalt = row[5].split(" -> ")[1].split(" (")[0]
        BASEchange = "{0}>{1}".format(BASEref, BASEalt)
    elif row[3] == "D":
        BASEref = row[5].split(" (")[0]
        BASEalt = "*"
        BASEchange = "{0}>{1}".format(BASEref, BASEalt)
    elif row[3] == "I (Dup)":
        BASEref = "*"
        BASEalt = row[5].split(" (")[0]
        BASEchange = "{0}>{1}".format(BASEref, BASEalt)
    else:
        BASEchange = row[5]
    if row[6] is np.nan:
        AAchange = row[6]
    elif "->" in row[6]:
        AAref = row[6].split(" -> ")[0]
        AAalt = row[6].split(" -> ")[1].split(" (")[0]
        AAchange = "{0}>{1}".format(AAref, AAalt)
    elif "[STOP]" in row[6]:
        AAchange = row[6].split(" (")[0].replace(" ", "_")
    else:
        AAchange = row[6]
    JSIaltVAF = row[4].split(" ")[0]
    JSIaltCount = row[4].split(" ")[1]
    return(sampleID, pcrID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF, JSIaltCount)
# Define function that extracts information from PCR1 and PCR2 dictionaries to write to file.
def getVAFtoWrite(variant):
    if (variant in variant_PCR1_dict.keys()) and (variant in variant_PCR2_dict.keys()):
        JSIaltVAF_PCR1, JSIaltCount_PCR1 = variant_PCR1_dict[variant].split("||")
        JSIaltVAF_PCR2, JSIaltCount_PCR2 = variant_PCR2_dict[variant].split("||")
    elif variant in variant_PCR1_dict.keys():
        JSIaltVAF_PCR1, JSIaltCount_PCR1 = variant_PCR1_dict[variant].split("||")
        JSIaltVAF_PCR2, JSIaltCount_PCR2 = "NA", "NA"
    elif variant in variant_PCR2_dict.keys():
        JSIaltVAF_PCR1, JSIaltCount_PCR1 = "NA", "NA"
        JSIaltVAF_PCR2, JSIaltCount_PCR2 = variant_PCR2_dict[variant].split("||")
    return(JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2)

print(dt.datetime.now())
print("Start.")

# Included sampleIDs list. NB. will be supplied on command line in process.
sampleIDs_included_list = set((np.loadtxt(INCL_SAMPLES, dtype='string')).tolist())

# Import variants without redundant metadata columns.
JSImerged_df = pd.read_csv(INPUT_JSIMERGED, header = 0, sep = "\t", index_col = False, usecols = ["Sample", "Gene", "Pos.", "Type", "Coverage", "Nuc Change", "AA Change", "c.HGVS", "p.HGVS"])

# Dictionary to store variants in.
variant_PCR1_dict = {}
variant_PCR2_dict = {}
variant_list = []

# Recode information and store for filtering.
for index, row in JSImerged_df.iterrows():
    sampleID, pcrID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF, JSIaltCount = recodeJSIinfo(row)
    variant_key = "{0}:{1}|{2}|{3}|{4}|{5}|{6}|{7}|{8}".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange)
    variant_PCR_value = "{0}||{1}".format(JSIaltVAF, JSIaltCount)
    variant_list.append(variant_key)
    if pcrID == "PCR1":
        variant_PCR1_dict[variant_key] = variant_PCR_value
    elif pcrID == "PCR2":
        variant_PCR2_dict[variant_key] = variant_PCR_value

# Store duplicate variants (identified in both PCR1 and PCR2) in list.
duplicates_list = variant_PCR1_dict.viewkeys() & variant_PCR2_dict.viewkeys()

# Open new coverage, PTPN11, Nalt, singletons and duplicates files to write to.
JSImerged_recoded_open = open(OUTPUT_JSIMERGED_RECODED, "w")
JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSImerged_recoded_coverage_open = open(OUTPUT_JSIMERGED_RECODED_EXCL_COVERAGE, "w")
JSImerged_recoded_coverage_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSImerged_recoded_PTPN11_open = open(OUTPUT_JSIMERGED_RECODED_EXCL_PTPN11, "w")
JSImerged_recoded_PTPN11_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSImerged_recoded_Nalt_open = open(OUTPUT_JSIMERGED_RECODED_EXCL_NALT, "w")
JSImerged_recoded_Nalt_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSImerged_recoded_singletons_open = open(OUTPUT_JSIMERGED_RECODED_EXCL_SINGLETONS, "w")
JSImerged_recoded_singletons_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))
JSImerged_recoded_duplicates_open = open(OUTPUT_JSIMERGED_RECODED_DUPLICATES, "w")
JSImerged_recoded_duplicates_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format("sampleID", "chromosome", "bpStart", "bpEnd", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2"))

# Loop over variant_list and write each variant to appropriate file.
for variant in set(variant_list):
    sampleID, pcrID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = None, None, None, None, None, None, None, None, None, None, None, None, None, None
    sampleID, variantInfo = variant.split(":")
    chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange = variantInfo.split("|")
    if sampleID not in sampleIDs_included_list:
        JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = getVAFtoWrite(variant)
        JSImerged_recoded_coverage_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
        JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
    elif gene == "PTPN11":
        JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = getVAFtoWrite(variant)
        JSImerged_recoded_PTPN11_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
        JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
    elif (cHGVS is np.nan) or ("N" in cHGVS):
        JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = getVAFtoWrite(variant)
        JSImerged_recoded_Nalt_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
        JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
    elif variant in duplicates_list:
        JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = getVAFtoWrite(variant)
        JSImerged_recoded_duplicates_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
        JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
    else:
        JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2 = getVAFtoWrite(variant)
        JSImerged_recoded_singletons_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))
        JSImerged_recoded_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2))

# Close open files.
JSImerged_recoded_open.close()
JSImerged_recoded_PTPN11_open.close()
JSImerged_recoded_Nalt_open.close()
JSImerged_recoded_singletons_open.close()
JSImerged_recoded_duplicates_open.close()
JSImerged_recoded_coverage_open.close()

print("Done.")
print(dt.datetime.now())
