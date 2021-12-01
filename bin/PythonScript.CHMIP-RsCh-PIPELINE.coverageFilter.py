## Python script that checks whether coverage is above THRESHOLD for both PCRs - Nextflow executable.

# Import  modules.
import sys
import pandas as pd
import numpy as np

# Define parameters.
INPUT_AVERAGE_COVERAGES = str(sys.argv[1])
INPUT_THRESHOLD = int(sys.argv[2])
OUTPUT_COVERAGE_PASS = str(sys.argv[3])

# Define function that extracts information from average coverages.
def getAverageCoverageInfo(row):
    sampleID = row[0]
    pcrID = row[1]
    panelAverage = int(row[2])
    return(sampleID, pcrID, panelAverage)

# Open average coverages.
averageCoverages_df = pd.read_csv(INPUT_AVERAGE_COVERAGES, sep = "\t", header = 0)

# Iterate over file and store sampleIDs with coverage INPUT_THRESHOLD THRESHOLD in PCR1 and PCR2 lists.
PCR1_PASS_list = []
PCR2_PASS_list = []
for index, row in averageCoverages_df.iterrows():
    sampleID, pcrID, panelAverage = getAverageCoverageInfo(row)
    if int(panelAverage) > INPUT_THRESHOLD:
        if pcrID == "PCR1":
            PCR1_PASS_list.append(sampleID)
        elif pcrID == "PCR2":
            PCR2_PASS_list.append(sampleID)

# Intersection between two lists.
sampleIDs_PASS_list = list(set(PCR1_PASS_list) & set(PCR2_PASS_list))

# Write sampleIDs that pass to file.
coveragePass_open = open(OUTPUT_COVERAGE_PASS, "w")
for sampleID in sampleIDs_PASS_list:
    print(sampleID)
    coveragePass_open.write("{0}\n".format(sampleID))

# Close file.
coveragePass_open.close()
