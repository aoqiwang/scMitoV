import pysam
import argparse

# Function to check if a read is spliced
def is_spliced(read):
    for op, length in read.cigar:
        if op == 3:  # 3 corresponds to 'N' in the CIGAR string
            return True
    return False

# Parse command line arguments
parser = argparse.ArgumentParser(description='Filter BAM files.')
parser.add_argument('-m', '--mt_bam', required=True, help='Path to the mitochondrial BAM file')
parser.add_argument('-n', '--nuclear_bam', required=True, help='Path to the nuclear BAM file')
parser.add_argument('-o', '--output_bam', required=True, help='Path to the output BAM file')
args = parser.parse_args()

# Open the BAM files
nuclear_bam = pysam.AlignmentFile(args.nuclear_bam, "rb")
mt_bam = pysam.AlignmentFile(args.mt_bam, "rb")

# Create a dictionary to store the mapping quality and splicing status of each read in the nuclear genome
nuclear_reads = {}
for read in nuclear_bam:
    # Skip secondary alignments, indicating multi-mapped reads
    if read.is_secondary or read.is_supplementary:
        continue
    nuclear_reads[read.query_name] = {"mapping_quality": read.mapping_quality, "is_spliced": is_spliced(read)}

# Open a new BAM file to write the reads that are better mapped to the mitochondrial genome
filtered_mt_bam = pysam.AlignmentFile(args.output_bam, "wb", template=mt_bam)

# Iterate over the reads in the mitochondrial BAM file
for read in mt_bam:
    # If the read is multi-mapped on the mitochondrial genome, skip it
    if read.is_secondary or read.is_supplementary:
        continue
    # Get the nuclear read
    nuclear_read = nuclear_reads.get(read.query_name)
    if nuclear_read is not None:
        mt_spliced = is_spliced(read)
        nuclear_spliced = nuclear_read["is_spliced"]
        # If the read is spliced only in one genome, assign it to the other genome
        if mt_spliced and not nuclear_spliced:
            continue
        if nuclear_spliced and not mt_spliced:
            filtered_mt_bam.write(read)
            continue
        # If the read is spliced in both genomes or not spliced in both, assign it based on mapping quality
        if nuclear_read["mapping_quality"] > read.mapping_quality:
            continue
    # If the read is not in the nuclear reads, or if it is better mapped to the mitochondrial genome, write it to the new BAM file
    filtered_mt_bam.write(read)

# Close the BAM files
nuclear_bam.close()
mt_bam.close()
filtered_mt_bam.close()
