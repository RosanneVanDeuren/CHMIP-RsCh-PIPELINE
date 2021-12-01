## Python script to combine mpileup VAFs with final JSI variants - Nextflow executable.

# Import modules.
import pandas as pd
import numpy as np
import sys

# Define parameters.
INPUT_PILEUPS_CONCLUDED = str(sys.argv[1])
INPUT_JSI_FINAL = str(sys.argv[2])
OUTPUT_JSI_FINAL_PILEUP_VAF = str(sys.argv[3])
OUTPUT_JSI_FINAL_NO_PILEUP_INSPECT_IGV = str(sys.argv[4])

# Define function that extracts mpileups concluded information.
def getMpileupsConcludedInfo(row):
    sampleID = row[0]
    pcrID = row[1]
    chromosome = row[2]
    bpStart = row[3]
    bpEnd = row[4]
    ref = row[5]
    alt = row[6]
    gene = row[7]
    conclusion = row[8]
    VAF = row[9]
    ALTcount = row[10]
    TOTALcount = row[11]
    note = row[12]
    return(sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, conclusion, VAF, ALTcount, TOTALcount, note)
# Define function that extracts variant_key from final JSI variants.
def getJSIvariantKey(row):
    sampleID = row[0]
    chromosome = row[1]
    bpStart = row[2]
    bpEnd = row[3]
    gene = row[4]
    BASEchange = row[7]
    ref, alt = BASEchange.split(">")
    if ref == "*": # insertion.
        CHROMOSOME = chromosome
        BPstart = bpStart
        BPend = bpStart
        REF = "nan"
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
    JSIkey = "{0}:{1}|{2}|{3}|{4}|{5}|{6}".format(sampleID, CHROMOSOME, BPstart, BPend, REF, ALT, gene)
    return(JSIkey)
# Define function that extracts information from final JSI variants including mpileupVAFs.
def getJSIfinalInfo(row):
    sampleID = row[0]
    chromosome = row[1]
    bpStart = row[2]
    bpEnd = row[3]
    gene = row[4]
    cHGVS = row[5]
    pHGVS = row[6]
    BASEchange = row[7]
    AAchange = row[8]
    JSIaltVAF_PCR1 = row[9]
    JSIaltCount_PCR1 = row[10]
    JSIaltVAF_PCR2 = row[11]
    JSIaltCount_PCR2 = row[12]
    mpileupConclusion_PCR2 = row[13]
    mpileupVAF_PCR2 = row[14]
    mpileupALTcount_PCR2 = row[15]
    mpileupTOTALcount_PCR2 = row[16]
    mpileupNote_PCR2 = row[17]
    mpileupConclusion_PCR1 = row[18]
    mpileupVAF_PCR1 = row[19]
    mpileupALTcount_PCR1 = row[20]
    mpileupTOTALcount_PCR1 = row[21]
    mpileupNote_PCR1 = row[22]
    return(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2, mpileupConclusion_PCR1, mpileupVAF_PCR1, mpileupALTcount_PCR1, mpileupTOTALcount_PCR1, mpileupNote_PCR1, mpileupConclusion_PCR2, mpileupVAF_PCR2, mpileupALTcount_PCR2, mpileupTOTALcount_PCR2, mpileupNote_PCR2)

# Open mpileups concluded.
mpileups_concluded_df = pd.read_csv(INPUT_PILEUPS_CONCLUDED, header = 0, sep = "\t")

# Iterate over concluded mpileups, compute averages and store in dictionary.
mpileups_concluded_dict = {}
for index, row in mpileups_concluded_df.iterrows():
    sampleID, pcrID, chromosome, bpStart, bpEnd, ref, alt, gene, conclusion, VAF, ALTcount, TOTALcount, note = getMpileupsConcludedInfo(row)
    variant_key = "{0}:{1}|{2}|{3}|{4}|{5}|{6}".format(sampleID, chromosome, bpStart, bpEnd, ref, alt, gene)
    variant_value = "{0}:{1}|{2}|{3}|{4}|{5}".format(pcrID, conclusion, VAF, ALTcount, TOTALcount, note)
    mpileups_concluded_dict.setdefault(variant_key, []).append(variant_value)

# Open final JSI variants.
JSIfinal_df = pd.read_csv(INPUT_JSI_FINAL, header = 0, sep = "\t")

print(mpileups_concluded_dict)

# Iterate over final JSI variants and add mpileup information.
for index, row in JSIfinal_df.iterrows():
    JSIkey = getJSIvariantKey(row)
    for valuePCR in mpileups_concluded_dict[JSIkey]:
        pcrID, rest = valuePCR.split(":")
        conclusion, VAF, ALTcount, TOTALcount, note = rest.split("|")
        if pcrID == "PCR1":
            JSIfinal_df.loc[index, "mpileupConclusion_PCR1"] = conclusion
            JSIfinal_df.loc[index, "mpileupVAF_PCR1"] = VAF
            JSIfinal_df.loc[index, "mpileupALTcount_PCR1"] = ALTcount
            JSIfinal_df.loc[index, "mpileupTOTALcount_PCR1"] = TOTALcount
            JSIfinal_df.loc[index, "mpileupNote_PCR1"] = note
        elif pcrID == "PCR2":
            JSIfinal_df.loc[index, "mpileupConclusion_PCR2"] = conclusion
            JSIfinal_df.loc[index, "mpileupVAF_PCR2"] = VAF
            JSIfinal_df.loc[index, "mpileupALTcount_PCR2"] = ALTcount
            JSIfinal_df.loc[index, "mpileupTOTALcount_PCR2"] = TOTALcount
            JSIfinal_df.loc[index, "mpileupNote_PCR2"] = note

# Iterate once more over final JSI variants, compute averages and write output files.
JSIfinal_mpileupVAFs_open = open(OUTPUT_JSI_FINAL_PILEUP_VAF, "w")
JSIfinal_mpileupVAFs_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\n".format("sampleID", "chromosome", "BPstart", "BPend", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2", "mpileupConclusion_PCR1", "mpileupVAF_PCR1", "mpileupALTcount_PCR1", "mpileupTOTALcount_PCR1", "mpileupNote_PCR1", "mpileupConclusion_PCR2", "mpileupVAF_PCR2", "mpileupALTcount_PCR2", "mpileupTOTALcount_PCR2", "mpileupNote_PCR2", "finalALTcount", "finalTOTALcount", "finalVAF"))
JSIfinal_noMpileupVAF_inspectIGV_open = open(OUTPUT_JSI_FINAL_NO_PILEUP_INSPECT_IGV, "w")
JSIfinal_noMpileupVAF_inspectIGV_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\n".format("sampleID", "chromosome", "BPstart", "BPend", "gene", "cHGVS", "pHGVS", "BASEchange", "AAchange", "JSIaltVAF_PCR1", "JSIaltCount_PCR1", "JSIaltVAF_PCR2", "JSIaltCount_PCR2", "mpileupConclusion_PCR1", "mpileupVAF_PCR1", "mpileupALTcount_PCR1", "mpileupTOTALcount_PCR1", "mpileupNote_PCR1", "mpileupConclusion_PCR2", "mpileupVAF_PCR2", "mpileupALTcount_PCR2", "mpileupTOTALcount_PCR2", "mpileupNote_PCR2"))

for index, row in JSIfinal_df.iterrows():
    sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2, mpileupConclusion_PCR1, mpileupVAF_PCR1, mpileupALTcount_PCR1, mpileupTOTALcount_PCR1, mpileupNote_PCR1, mpileupConclusion_PCR2, mpileupVAF_PCR2, mpileupALTcount_PCR2, mpileupTOTALcount_PCR2, mpileupNote_PCR2 = getJSIfinalInfo(row)
    if (mpileupConclusion_PCR1 == "FP") and (mpileupConclusion_PCR2 == "FP"):
        JSIfinal_noMpileupVAF_inspectIGV_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2, mpileupConclusion_PCR1, mpileupVAF_PCR1, mpileupALTcount_PCR1, mpileupTOTALcount_PCR1, mpileupNote_PCR1, mpileupConclusion_PCR2, mpileupVAF_PCR2, mpileupALTcount_PCR2, mpileupTOTALcount_PCR2, mpileupNote_PCR2))
    else:
        finalALTcount = round((float(mpileupALTcount_PCR1)+int(mpileupALTcount_PCR2))/2, 0)
        finalTOTALcount = round((float(mpileupTOTALcount_PCR1)+int(mpileupTOTALcount_PCR2))/2, 0)
        finalVAF = round((float(finalALTcount)/finalTOTALcount)*100, 2)
        JSIfinal_mpileupVAFs_open.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\n".format(sampleID, chromosome, bpStart, bpEnd, gene, cHGVS, pHGVS, BASEchange, AAchange, JSIaltVAF_PCR1, JSIaltCount_PCR1, JSIaltVAF_PCR2, JSIaltCount_PCR2, mpileupConclusion_PCR1, mpileupVAF_PCR1, mpileupALTcount_PCR1, mpileupTOTALcount_PCR1, mpileupNote_PCR1, mpileupConclusion_PCR2, mpileupVAF_PCR2, mpileupALTcount_PCR2, mpileupTOTALcount_PCR2, mpileupNote_PCR2, finalALTcount, finalTOTALcount, finalVAF))

# Close files.
JSIfinal_mpileupVAFs_open.close()
JSIfinal_noMpileupVAF_inspectIGV_open.close()
