## Python script to compute the average coverage depth over the entire panel and exclude sampleIDs with coverage below THRESHOLD - Nextflow executable.

# Import modules.
import sys
import pandas as pd
import numpy as np

# Define parameters.
INPUT_PARSED_COVERAGE = str(sys.argv[1])
INPUT_FILENAME = str(sys.argv[2])
OUTPUT_AVERAGE_COVERAGE = str(sys.argv[3])

# Open parsedCoverage-file.
parsedCovDepth_df = pd.read_csv(INPUT_PARSED_COVERAGE, sep = "\t", header = 0)

# Define sampleID and pcrID from fileName, and THRESHOLD from parameter.
fileName = INPUT_FILENAME
sampleID, rest = fileName.split("-")
pcrID = rest.split(".")[0]

# Compute average over the entire panel.
panelAverage = round(parsedCovDepth_df["CoverageDepth"].mean(), 2)

# Store average coverage over the entire panel in average-file.
panelAverageFile_open = open(OUTPUT_AVERAGE_COVERAGE, "w")
panelAverageFile_open.write("{0}\t{1}\t{2}\n".format("sampleID", "pcrID", "panelAverage"))
panelAverageFile_open.write("{0}\t{1}\t{2}\n".format(sampleID, pcrID, panelAverage))
panelAverageFile_open.close()
