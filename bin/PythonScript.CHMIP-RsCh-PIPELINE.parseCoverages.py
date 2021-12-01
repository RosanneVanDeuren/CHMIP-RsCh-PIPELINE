## Python script that parses and sorts raw bedtools coverage depth - Nextflow executable.

# Import  modules.
import sys
import pandas as pd
import numpy as np

# Define parameters.
INPUT_RAW_COVERAGE = str(sys.argv[1])
OUTPUT_PARSED_COVERAGE = str(sys.argv[2])
OUTPUT_PARSED_SORTED_COVERAGE = str(sys.argv[3])

# Define function that extracts raw coverage depth information.
def getRawCovInfo(row):
    chrom = row[0]
    startMIP = row[1]
    endMIP = row[2]
    nameMIP = row[3]
    strandMIP = row[4]
    ixMIP = row[5]
    covDepth = row[6]
    return(chrom, startMIP, endMIP, nameMIP, strandMIP, ixMIP, covDepth)

# Open raw bedtools coverage depth.
rawCovDepth_df = pd.read_csv(INPUT_RAW_COVERAGE, sep = "\t", header = None)

# Create dictionary to store coverage and MIP information.
coverageDict = {}

# Loop over raw coverage dataframe and store coverage information in dictionary.
for index, row in rawCovDepth_df.iterrows():
    chrom, startMIP, endMIP, nameMIP, strandMIP, ixMIP, covDepth = getRawCovInfo(row)
    bp = (startMIP + (ixMIP - 1))
    MIPinfo  = "{0}|({1})|[{2}-{3}]".format(nameMIP, strandMIP, startMIP, endMIP)
    covKey = "{0}:{1}".format(chrom, bp)
    covVal = "{0};{1}".format(covDepth, MIPinfo)
    coverageDict.setdefault(covKey, []).append(covVal)

# Create temporary unsorted output-file to write to & write header:
parsedCoverage_open = open(OUTPUT_PARSED_COVERAGE, "w")
parsedCoverage_open.write("{0}\t{1}\t{2}\t{3}\n".format("Chromosome", "Position", "CoverageDepth", "MIPs"))

# Loop over coverage dictionary and write information to parsed file.
for key in coverageDict:
    chromosome = key.split(':')[0]
    position  = key.split(':')[1]
    samePosMIPinfo = {}
    for item in coverageDict[key]:
        coverage = item.split(';')[0]
        newKey = "{0}:{1}:{2}".format(chromosome, position, coverage)
        posMIPinfo = item.split(';')[1]
        samePosMIPinfo.setdefault(newKey, []).append(posMIPinfo)
    linetowrite = '{0}\t{1}\t{2}\t{3}\n'.format(chromosome, position, coverage, ', '.join(samePosMIPinfo[newKey]))
    parsedCoverage_open.write(linetowrite)

# Close temporary unsorted output-file.
parsedCoverage_open.close()

# Import temporary unsorted coverage as pandas dataframe.
parsedCovDepth_df = pd.read_csv(OUTPUT_PARSED_COVERAGE, sep = "\t", header = 0)

# Compute noChromosome-column, replace X with 23, and convert to numeric for sorting.
parsedCovDepth_df["noChromosome"] = parsedCovDepth_df.Chromosome.str[3:]
parsedCovDepth_df.loc[parsedCovDepth_df["noChromosome"] == "X", "noChromosome"] = 23
parsedCovDepth_df["noChromosome"] = pd.to_numeric(parsedCovDepth_df["noChromosome"])

# Sort on 1) noChromosome, and 2) position.
parsedCovDepth_df = parsedCovDepth_df.sort_values(["noChromosome", "Position"], ascending=[True, True])

# Remove noChromosome-column.
parsedCovDepth_df = parsedCovDepth_df.drop(["noChromosome"], axis = 1)

# Write parsed sorted coverage to file.
parsedCovDepth_df.to_csv(OUTPUT_PARSED_SORTED_COVERAGE, sep = "\t", header = True, index = False)
