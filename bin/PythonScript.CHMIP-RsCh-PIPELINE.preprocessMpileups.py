## Python script to prefilter mpileups before parsing by perl - Nextflow executable.

# Import relevant modules.
import sys

# Define parameters.
INPUT_PILEUP = str(sys.argv[1])
OUTPUT_PILEUP_PREPPED = str(sys.argv[2])
OUTPUT_PILEUP_READMARKERS = str(sys.argv[3])

# Define function that extracts information from mpileup file.
def getinfompileup(line):
    colmns = line.strip('\n').split('\t')
    pileup_chr = colmns[0]
    pileup_bp = int(colmns[1])
    pileup_ref = colmns[2]
    pileup_rc = int(colmns[3])
    pileup_seq = (colmns[4]).upper().replace(',', '.')
    pileup_bq = colmns[5]
    return(pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq)
# Define function that writes to readmarker files.
def writeReadMarker(pileup_seq):
    if '^' in pileup_seq:
        readMarker_open.write('{0}\t{1}\t{2}\t{3}\n'.format(pileup_chr, pileup_bp, pileup_ref, "readstart"))
    elif '$' in pileup_seq:
        readMarker_open.write('{0}\t{1}\t{2}\t{3}\n'.format(pileup_chr, pileup_bp, pileup_ref, "readend"))

# Open mpileups.merged file.
pileup_open = open(INPUT_PILEUP, "r")
# Create mpileups.merged.prepped to write to file.
pileupPrepped_open = open(OUTPUT_PILEUP_PREPPED, "w")
# Create readmarker file to write to.
readMarker_open = open(OUTPUT_PILEUP_READMARKERS, "w")
readMarker_open.write("{0}\t{1}\t{2}\t{3}\n".format("chromosome", "bpposition", "ref", "mpileupSeqMarker"))

# Process mpileup.
for line in pileup_open:
    pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq = getinfompileup(line)
    if '^' in pileup_seq or '$' in pileup_seq:
        writeReadMarker(pileup_seq)
        pileupPrepped_open.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq))
    else:
        pileupPrepped_open.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq))

# Close files.
pileupPrepped_open.close()
pileup_open.close()
readMarker_open.close()
