# smRandom-seq Protocol Analysis Workflow

This guide provides a comprehensive, step-by-step workflow for analyzing single-bacteria total-RNA-seq data using smRandom-seq. It is applicable to various sample types, including laboratory-cultured microbial samples, human feces, and bovine rumen fluid.

## Contents
1. [Setup and Dependencies](#1-setup-and-dependencies)
2. [Sequencing Data Pre-processing](#2-sequencing-data-pre-processing)
3. [Count Matrix Generation](#3-count-matrix-generation)
4. [Objects Generation and Process](#4-objects-generation-and-process)
---

### 1. Setup and Dependencies
Ensure the following dependencies are installed:
#### 1.1 Required for all samples：
- **Python 3.8+**: For running custom scripts.
- **R 4.0+**: For statistical analysis and visualization.
- **Samtools 1.15**: (https://github.com/samtools/samtools) For sorting and indexing BAM files.
- **bedtools v2.29.2**: For genomic data manipulation.
- **STAR 2.7.10a**:(https://github.com/alexdobin/STAR) For genome indexing and read alignment.
- **FeatureCounts 2.0.3**: (https://sourceforge.net/projects/subread) For exon-based gene expression counting.
- **Umi_tools 1.1.2**: (https://github.com/CGATOxford/UMI-tools) For UMI-based read counting.
#### 1.2 Additional Dependencies for Complex Microbial Samples：
- **Kraken2 2.1.2**: (https://github.com/DerrickWood/kraken2) For K-mer-based taxonomic classification.
- **Bracken 2.6.0**: (https://github.com/jenniferlu717/Bracken) For calculating taxonomic abundance after Kraken2 classification.
Note: Ensure all software versions are compatible. Update tools or databases as needed for optimal performance.

### 2. Sequencing Data Pre-processing
Process raw FASTQ data with Raw_data_preprocessing.sh using the custom tool anchoradp.o to separate reads, extract barcodes and UMIs, and prepare data for downstream analyses.
---
#### 2.1 About `anchoradp.o` and `Raw_data_preprocessing.sh`
`anchoradp.o`: Processes raw FASTQ data, including separating reads, identifying PCR adapters, and extracting barcodes and UMIs.
`Raw_data_preprocessing.sh`: A wrapper script that uses anchoradp.o for streamlined data preprocessing.
The `Raw_data_preprocessing.sh` Script Content:
```bash
# Raw data preprocessing
anchoradp.o -1 <(zcat R1) -2 <(zcat R2) -p -o your_prefix
```
Parameters:
`-1 <(zcat R1)`: Unzips R1.fastq.gz and passes the decompressed file to anchoradp.o as the first input. 
`-2 <(zcat R2)`: Unzips R2.fastq.gz and passes the decompressed file to anchoradp.o as the second input. 
`-p`: Enables parallel processing to improve speed (if supported).
`-o` your_prefix: Specifies the prefix for output files.

#### 2.2 Running the Script:
Use the following command to preprocess the raw FASTQ data:
```bash
bash `Raw_data_preprocessing.sh` R1.fastq.gz R2.fastq.gz your_prefix
```
#### Intput:
R1.fastq.gz and R2.fastq.gz: Illumina PE150 paired-end sequencing raw data that have already undergone adapter trimming. 
`your_prefix`: The prefix for the output files.

Compress FASTQ Files:
After processing, compress the output FASTQ files to prepare for alignment. Ensure gzip is installed in your environment:
```bash
gzip your_prefix_1.fq
gzip your_prefix_2.fq
```

#### Output:
`your_prefix_1.fq.gz`: Contains compressed cell barcodes (20bp) and UMIs (8bp).
`your_prefix_2.fq.gz`: Contains compressed cDNA reads.
The sequence IDs in the two output files are aligned.

### 3. Count Matrix Generation
The workflow for genome index creation, read alignment, and count matrix generation varies by sample type.
---
### 3.1 Laboratory-Cultured Microbial Samples
Use the `Count_matrix.sh` script to perform all steps of genome index creation, read alignment, and gene expression counting. This script integrates all necessary commands for STAR indexing, alignment, feature counting, and barcode selection into a single workflow.
About `Count_matrix.sh` script details:：
Genome Index Creation: Uses the STAR aligner to generate genome index files from the input genome FASTA and GTF annotation files. `sjdbOverhang`: Set to 122 for 150-bp paired-end reads. `genomeDir`: Directory where index files are stored.
Read Alignment: Aligns cDNA sequencing reads `your_prefix_2.fq.gz` to the indexed genome. Outputs sorted BAM files for downstream analysis.
Gene Expression Counting: Uses `featureCounts` to count gene-level expression by mapping aligned reads to exon regions in the genome annotation. Produces a BAM file containing unique gene alignments.
Sorting and Indexing BAM Files: `samtools` is used to sort and index the BAM files for fast access.
UMI-Based Read Counting: Extracts unique molecular identifiers (UMIs) and calculates gene expression per cell using `umi_tools`. Generates a wide-format gene expression matrix.
Cell Barcode Selection: Uses the Python script `selectResult.py` to filter the expression matrix based on a user-defined cell threshold (`thres`).

Running `Count_matrix.sh`
Use the following command to run the script with the required input parameters:
```bash
bash Count_matrix.sh genomeDir=$genomeDir genomeFastaFile=$genomeFastaFile sjdbGTFfile=$sjdbGTFfile your_prefix=$your_prefix thres=$thres
```
#### Input Parameters:
genomeDir: Directory for STAR genome index files.
genomeFastaFile: Path to the genome FASTA file `your_prefix_2.fq.gz`.
sjdbGTFfile: Path to the GTF file containing gene annotations.
your_prefix: The prefix for the output files.
thres: Threshold for cell count, usually between 50,00 and 50,000.
#### Output:
`$your_prefix.counts.tsv`: Gene expression matrix.
`$your_prefix.selected`: Filtered cell barcodes based on the threshold.

### 3.2 Complex Microbial Samples
Use the `MICtools` pipeline for count matrix generation of human gut microbiome single-microbe RNA sequencing. MICtools provides modules for taxonomic annotation (`MIC-Anno`) and transcriptional matrix construction (`MIC-Bac`).

#### Barcode Selection
`Barcodes_selection.py` filters barcodes from the processed FASTQ files based on a threshold set in the script. This threshold determines the minimum UMI count per barcode (e.g., 500 in this example). Users can adjust this threshold or include all barcodes if desired.
Running the Script:
```bash
python Barcodes_selection.py
```
#### Intput:
`your_prefix_1.fq.gz`
`your_prefix_2.fq.gz`
#### Output:
`your_prefix_1_extracted.fq.gz`: Filtered barcodes and UMIs.
`your_prefix_2_extracted.fq.gz`: Filtered cDNA sequences.
Notes: 
1. The threshold can be adjusted within Barcodes_selection.py by modifying the threshold variable in the script. A default value of 500 UMIs is used in the provided example.
2. This step filters out low-quality or low-complexity barcodes (e.g., with very few UMIs), thereby enhancing the speed of downstream analysis.

#### Taxonomic Annotation
#### Using MICtools for Human Gut Microbiome and Other Samples
`MIC-Anno` uses a K-mer-based taxonomic classification strategy with Kraken2 to annotate microbes at the species and genus levels.
Running the Command:
```bash
MICtools anno --module pipeline -s your_prefix_2_extracted.fq.gz -r [kraken_ref_filepath] -p [your_prefix]
```
#### Intput:
`your_prefix_2_extracted.fq.gz`: FASTQ file containing extracted cDNA sequences.
`kraken_ref_filepath`: Path to the Kraken2 reference database.
For human gut microbiome samples, use the UHGG gut microbiome genome database. For other complex microbial samples, replace with an appropriate Kraken2 database for your community.
#### Output:
`your_prefix.barcode_count.txt`: Main single microbe taxonomic annotation result and the meanings of each column.
`your_prefix.genus_info.txt`: Genus info of each microbe.
`your_prefix.species_info.txt`: Species info of each microbe.
`your_prefix_genus.pdf`: Pieplot of sample genus composition.
`your_prefix_species.pdf`: Pieplot of sample specie composition.

#### For Bovine Rumen Fluid Samples
Taxonomic annotation for bovine rumen fluid samples can utilize a Kraken2-based gOTUs (genomic operational taxonomic units) database. The Bovine Gastrointestinal Microbial Genome Map (BGMGM) is available for download from Figshare at https://figshare.com/articles/dataset/Microbiome_single-cell_transcriptomics_reveal_functional_heterogeneity_of_metabolic_niches_covering_more_than_2_500_species_in_the_rumen/24844344. 
#### Steps to Create a Kraken2-Based gOTUs Database:
1. Mask Ribosomal RNA Genes: Use Barrnap and bedtools to mask ribosomal RNA regions to avoid misclassification.
```bash
barrnap <reference_genome.fasta> > rRNA_annotation.gff
bedtools maskfasta -fi <reference_genome.fasta> -bed rRNA_annotation.gff -fo masked_genomes.fasta
```
2. Build the Kraken2 database:
```bash
kraken2-build --no-masking --add-to-library masked_genomes.fasta --db kraken2_gOTUs_db
```
3. Classify Reads by Kraken2:
```bash
kraken2 --db kraken2_gOTUs_db --report your_prefix_kraken2_report.txt \
        --output your_prefix_kraken2_output.txt your_prefix_2_trimmed.fastq.gz
```
4. Calculate Taxonomic Abundance with Bracken:
```bash
bracken -d kraken2_gOTUs_db -i your_prefix_kraken2_report.txt -o your_prefix_bracken_output.txt
```
5. Classify cells based on informative reads:
Cells with >50% informative reads are considered accurately annotated.
Cells below this threshold are labeled with the __like flag after the species name.

#### Construction of Count Matrix:
To generate a gene expression matrix for single microbes, use the `MIC-Bac` tool. By default, the tool filters abundant bacterial species (>3% of total barcodes), though users can opt to include all genera.
Running the Command:
```bash
MICtools bac -s your_prefix_2_extracted.fq.gz -r [microbiome_ref_folder_path] -f [gene/exon/transcript...] -p your_prefix
```
#### Input:
`your_prefix_2_extracted.fq.gz`: Filtered cDNA sequences..
`microbiome_ref_folder_path`: Path to reference genome files and annotations for the microbial community. For human gut microbiome samples, use UHGG genome database references. For Bovine Rumen Fluid Samples, the bovine gastrointestinal microbial genome database (BGMGM) is available on Figshare at https://figshare.com/articles/dataset/Microbiome_single-cell_transcriptomics_reveal_functional_heterogeneity_of_metabolic_niches_covering_more_than_2_500_species_in_the_rumen/24844344. For other complex microbial samples, replace with appropriate reference genome folders.
`-f`: Feature type in the GTF file (e.g., transcript or CDS), specified in column 3 of the GTF file.
#### Output:
`your_prefix_counts.tsv`: Single-microbe transcriptional matrix.

#### Recommendations for Other Complex Microbial Samples
Replace UHGG human gut microbiome genome database with a reference database suited for the target microbial community. For bovine rumen fluid, use databases such as BGMGM (Bovine Gastrointestinal Microbial Genome Map). For other complex microbial samples, databases like RefSeq, GenBank, or custom-curated genomes may be required. Update reference paths in all steps to point to the new Kraken2 database, genome folders, and associated annotations.

### 4. Objects Generation and Process
This section describes how to process the single-bacteria gene expression matrix and perform downstream analysis for both laboratory-cultured microbial samples and complex microbial samples.

---
### For Laboratory-Cultured Microbial Samples:
#### Reading and Filtering:
Follow the general scRNA-seq code and standards to process the single-bacteria gene expression matrix. Load the matrix into R to generate Seurat objects or use Python to generate Scanpy objects. Perform quality control (QC) to filter low-quality cells and ensure data integrity before downstream analysis.
#### Analysis of Antibiotic-Treated Samples: 
For samples treated with antibiotics (as described in the smRandom-seq paper), the following scripts can be used:
`CIP_antibiotic_treatment.py`: analyze samples treated with ciprofloxacin.
`AMP_antibiotic_treatment.py`: analyze samples treated with ampicillin.
`GO_enrichment.R`: GO enrichment analysis.
#### Use the following scripts for validation performance Results in smRandom-seq paper:
`Correlation_of_gene_expression.R`: Gene expression correlation. 
`Saturation_analysis.sh`: Saturation analysis.
`Technical_reproducibility.py`: Technical reproducibility.

### For Complex Microbial Samples:
For complex microbial samples, microbiome transcriptional analysis focuses on major genera and involves dimensional reduction and clustering. Seurat (v4) can be used for these analyses, as demonstrated in the following scripts:
`microbiome_transcriptional_analysis.R`: Provides an example workflow for human feces samples, including clustering and visualization.
`Reading_and_Filtering.R`: An example for microbial further analysis of the Bovine Rumen Fluid Samples, including reading the single-cell gene expression matrix into R, creating Seurat objects, and performing quality control to filters out low-quality cells and retains high-quality cells. 
`MscT_initial_matrix.R`: An example for benchmarking and clustering on the Seurat objects of the Bovine Rumen Fluid Samples. 
`MscT_HMACs.R`: An example for identifing and isolates cells with high metabolic activity (HMACs) within the dataset. Re-clusters these cells to analyze specific metabolic functions and activities. 
`MscT_basfia.R`: An example to re-clusters cells of Basfia succiniciproducens to examine functional heterogeneity within this specific microbial species. 
`MscT_basfia_monocle.R`: An example to applies Monocle to analyze the developmental trajectory of cells in Basfia succiniciproducens, exploring cellular dynamics over pseudo-time.

