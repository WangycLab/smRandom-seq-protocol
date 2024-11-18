import gzip
from collections import defaultdict

# Input and output file names
barcode_file = "your_prefix_1.fq.gz"  # Contains barcodes and UMIs
cDNA_file = "your_prefix_2.fq.gz"     # Contains cDNA sequences
output_barcode_file = "your_prefix_1_extracted.fq.gz"
output_cDNA_file = "your_prefix_2_extracted.fq.gz"

# Parameters
barcode_length = 20
umi_length = 8

# Function to read FASTQ files
def read_fastq(file_path):
    with gzip.open(file_path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:  # End of file
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            yield header, sequence, plus, quality

# Extract barcodes and count UMIs
barcode_umi_counts = defaultdict(set)
for header, sequence, _, _ in read_fastq(barcode_file):
    barcode = sequence[:barcode_length]
    umi = sequence[barcode_length:barcode_length + umi_length]
    barcode_umi_counts[barcode].add(umi)

# Convert to a dictionary of barcode counts
barcode_counts = {barcode: len(umis) for barcode, umis in barcode_umi_counts.items()}

# Select barcodes above threshold (adjust threshold as needed)
threshold = 500  # Example threshold for barcodes
selected_barcodes = {barcode for barcode, count in barcode_counts.items() if count > threshold}

# Extract reads matching selected barcodes
def filter_fastq(input_file, output_file, selected_barcodes):
    with gzip.open(input_file, "rt") as infile, gzip.open(output_file, "wt") as outfile:
        while True:
            header = infile.readline().strip()
            if not header:  # End of file
                break
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            barcode = sequence[:barcode_length]
            if barcode in selected_barcodes:
                outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")

# Write filtered FASTQ files
filter_fastq(barcode_file, output_barcode_file, selected_barcodes)
filter_fastq(cDNA_file, output_cDNA_file, selected_barcodes)

print(f"Filtered barcodes saved to {output_barcode_file} and {output_cDNA_file}")
