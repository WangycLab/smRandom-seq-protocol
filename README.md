# smRandom-seq Protocol Analysis Workflow

This guide provides a comprehensive, step-by-step workflow for analyzing single-bacteria total-RNA-seq data using smRandom-seq, applicable to various sample types, including laboratory-cultured microbial samples, human feces, and bovine rumen fluid.
## Contents
1. [Setup and Dependencies](#1-setup-and-dependencies)
2. [Sequencing Data Pre-processing](#2-sequencing-data-pre-processing)
3. [Count Matrix Generation](#3-count-matrix-generation)
4. [Objects Generation and Process](#4-seurat-objects-generation-and-process)
---

### 1. Setup and Dependencies
Ensure the following dependencies are installed:
#### Required for all samples：
- **Python 3.8+**: For running custom scripts.
- **R 4.0+**: For statistical analysis and visualization.
- **Samtools 1.15**: (https://github.com/samtools/samtools) For sorting and indexing BAM files.
- **bedtools v2.29.2**: For genomic data manipulation.
- **STAR 2.7.10a**:(https://github.com/alexdobin/STAR) For genome indexing and read alignment.
- **FeatureCounts 2.0.3**: (https://sourceforge.net/projects/subread) For exon-based gene expression counting.
- **Umi_tools 1.1.2**: (https://github.com/CGATOxford/UMI-tools) For UMI-based read counting.
Dependencies Based on Sample Type:
#### For Complex Microbial Samples:
- **Kraken2 2.1.2**: (https://github.com/DerrickWood/kraken2) For K-mer-based taxonomic classification.
- **Bracken 2.6.0**: (https://github.com/jenniferlu717/Bracken) For calculating taxonomic abundance after Kraken2 classification.
  
### 2. Sequencing Data Pre-processing
Process raw FastQ data with `Raw_data_preprocessing.sh` using the custom tool `anchoradp.o` to separate reads, extract barcodes and UMIs, and prepare data for downstream analyses.

About `anchoradp.o`:
Processes raw FastQ data without sequencing adapters.
Separates Read 1 and Read 2, identifies PCR adapters, and extracts key sequences.

About `Raw_data_preprocessing.sh`：
Raw_data_preprocessing.sh uses the custom tool anchoradp.o to preprocess raw FastQ files.

#### Running the Script:
```bash
bash `Raw_data_preprocessing.sh` R1.fastq.gz R2.fastq.gz your_prefix
```
#### Compress FASTQ Files:
After processing, compress the output FASTQ files to prepare for alignment.
```bash
gzip your_prefix_1.fq
gzip your_prefix_2.fq
```
#### Output:
your_prefix_1.fq.gz: Contains compressed cell barcodes (20bp) and UMIs (8bp).
your_prefix_2.fq.gz: Contains compressed cDNA sequences.
The sequence IDs in your_prefix_1.fq.gz and your_prefix_2.fq.gz are aligned。

### 3. Count Matrix Generation
The workflow for genome index creation and read alignment and count matrix generation depends on the sample type.
### For Laboratory-Cultured Microbial Samples:
Use the `Count_matrix.sh` script to perform all steps of genome index creation, read alignment, and gene expression counting. This script integrates all necessary commands for STAR indexing, alignment, feature counting, and barcode selection into a single workflow.

About `Count_matrix.sh`：
Genome Index Creation: Uses the STAR aligner to generate genome index files from the input genome FASTA and GTF annotation files. `sjdbOverhang`: Set to 122 for 150-bp paired-end reads. `genomeDir`: Directory where index files are stored.
Read Alignment: Aligns cDNA sequencing reads `your_prefix_2.fq.gz` to the indexed genome. Outputs sorted BAM files for downstream analysis.
Gene Expression Counting: Uses `featureCounts` to count gene-level expression by mapping aligned reads to exon regions in the genome annotation. Produces a BAM file containing unique gene alignments.
Sorting and Indexing BAM Files: `samtools` is used to sort and index the BAM files for fast access.
UMI-Based Read Counting: Extracts unique molecular identifiers (UMIs) and calculates gene expression per cell using `umi_tools`. Generates a wide-format gene expression matrix.
Cell Barcode Selection: Uses the Python script `selectResult.py` to filter the expression matrix based on a user-defined cell threshold (thres).

Running `Count_matrix.sh`
Use the following command to run the script with the required input parameters:
```bash
bash Count_matrix.sh genomeDir=$genomeDir genomeFastaFile=$genomeFastaFile sjdbGTFfile=$sjdbGTFfile your_prefix=$your_prefix thres=$thres
```
#### Input Parameters:
genomeDir: Directory for STAR genome index files.
genomeFastaFile: Path to the genome FASTA file `your_prefix_2.fq.gz`.
sjdbGTFfile: Path to the GTF file containing gene annotations.
your_prefix: Prefix for the sample name (e.g., your_prefix.fq.extracted).
thres: Threshold for cell count, usually between 50,00 and 50,000.
#### Output:
Expression Matrix: `$your_prefix.counts.tsv` Contains gene expression values for individual cells.
Filtered Cell Barcodes: `$your_prefix.selected` Contains selected cell barcodes based on the specified threshold.

### For Complex Microbial Samples:
Use the `MICtools` pipeline for count matrix generation of human gut microbiome single-microbe RNA sequencing. MICtools provides modules for taxonomic annotation (`MIC-Anno`), transcriptional matrix construction (`MIC-Bac`), and host-phage transcriptional relationship analysis(`MIC-Phage`).

#### Barcode Selection
`Barcodes_selection.py` filters barcodes from FASTQ files. 
```bash
python Barcodes_selection.py
```
#### Intput:
`your_prefix_1.fq.gz`
`your_prefix_2.fq.gz`
#### Output:
`your_prefix_1_extracted.fq.gz`: Filtered barcodes and UMIs.
`your_prefix_2_extracted.fq.gz`: Filtered cDNA sequences.
Adjust Threshold: Modify the threshold variable in the script to set the minimum UMI count per barcode (e.g., 500 in this example).

#### Taxonomic Annotation
#### Using MICtools for Human Gut Microbiome and Other Samples
`MIC-Anno` uses a K-mer-based taxonomic classification strategy with Kraken2 to annotate microbes at the species and genus levels.
```bash
MICtools anno --module pipeline -s your_prefix_2.fq.gz -r [kraken_ref_filepath] -p [your_prefix]
```
#### Intput:
`your_prefix_2.fq.gz`: FASTQ file containing extracted cDNA sequences.
`kraken_ref_filepath`: Path to the Kraken2 reference database.
For human gut microbiome samples, use the UHGG gut microbiome genome database.
For other complex microbial samples, replace with an appropriate Kraken2 database for your community.
#### Output:
`your_prefix.barcode_count.txt`: Main single microbe taxonomic annotation result and the meanings of each column.
`your_prefix.genus_info.txt`: Genus info of each microbe.
`your_prefix.species_info.txt`: Species info of each microbe.
`your_prefix_genus.pdf`: Pieplot of sample genus composition.
`your_prefix_species.pdf`: Pieplot of sample specie composition.

#### For Bovine Rumen Fluid Samples, taxonomic annotation can use a Kraken2-Based gOTUs Database. 
The bovine gastrointestinal microbial genome database (Bovine Gastro Microbial Genome Map, BGMGM) is available on Figshare at https://figshare.com/articles/dataset/Microbiome_single-cell_transcriptomics_reveal_functional_heterogeneity_of_metabolic_niches_covering_more_than_2_500_species_in_the_rumen/24844344. 
Create a Kraken2-Based gOTUs Database
Mask ribosomal RNA genes using Barrnap and bedtools to avoid misclassification:
```bash
barrnap <reference_genome.fasta> > rRNA_annotation.gff
bedtools maskfasta -fi <reference_genome.fasta> -bed rRNA_annotation.gff -fo masked_genomes.fasta
```
Build the Kraken2 database:
```bash
kraken2-build --no-masking --add-to-library masked_genomes.fasta --db kraken2_gOTUs_db
```
Classify Reads by Kraken2:
```bash
kraken2 --db kraken2_gOTUs_db --report your_prefix_kraken2_report.txt \
        --output your_prefix_kraken2_output.txt your_prefix_2_trimmed.fastq.gz
```
Calculate Taxonomic Abundance with Bracken:
```bash
bracken -d kraken2_gOTUs_db -i your_prefix_kraken2_report.txt -o your_prefix_bracken_output.txt
```
Classify cells based on informative reads:
Cells with >50% informative reads are considered accurately annotated.
Cells below this threshold are labeled with the __like flag after the species name.

#### Use `MIC-Bac` to construct single microbe transcriptional matrix:
Abundant Species Filtering: Retains abundant bacterial species (>3% of total barcodes), though users can include all genera if desired.
Gene Expression Matrix Construction:STAR (v2.7.10a), featureCounts (v2.0.3), and umi_tools (v1.1.2) are used to generate bacterial gene expression matrices with UHGG as a reference.
```bash
$ MICtools bac -s your_prefix_2.fq.gz -r [microbiome_ref_folder_path] -f [gene/exon/transcript...]  -p [output_prefix]
```
Note: The feature type of gtf file (the input of premeter '-f') differs in different gtf files column#3. Here in the following example gtf file, we can identify 'transcript' and 'CDS' are main feature types, so we can set '-f transcript' or '-f CDS' to process the analysis.
#### Use `MIC-Phage` to construct host-phage transcriptional relationship matrix:
```bash
$ MICtools phage -s your_prefix_2.fq.gz -rr [sortmerna_rRNA_reference_datasets_path] -pr [phage_ref_folder_path] -f [gene/exon/transcript...]  -p [output_prefix]
```

### Bovine Rumen Fluid Samples:
Follow the steps below for quality control, read alignment, gene abundance calculation, and taxonomic classification.
#### Quality Control: 
Trim reverse reads: Forward reads are kept untrimmed. Raw reverse reads are trimmed using Trimmomatic with the following settings:
```bash
trimmomatic PE -phred33 your_prefix_1.fq.gz your_prefix_2.fq.gz \
your_prefix_1_trimmed.fastq.gz your_prefix_1_unpaired.fastq.gz your_prefix_2_trimmed.fastq.gz your_prefix_2_unpaired.fastq.gz \
SLIDINGWINDOW:4:15 MINLEN:50
```
Discards reads shorter than 50 bp after trimming.
Filter barcodes: Assign clean reads to cells based on barcodes and establish valid cell counts. Set read count thresholds based on sequencing depth, sample quality, and cell requirements (e.g., 2,000, 2,500, 5,000, or 6,000 reads per cell).
#### Gene Abundance Calculation:
Map Clean Reads to Gene Catalogs:
```bash
bwa mem -t N BGMGM.fasta your_prefix_2_trimmed.fastq.gz > your_prefix_aligned.sam
```
Extract Mapped Reads:
```bash
bedtools bamtobed -i your_prefix_aligned.sam > your_prefix_aligned.bed
```
Sorting and Indexing:
```bash
samtools sort -@ N your_prefix_aligned.bam -o your_prefix_aligned.sorted.bam
samtools index your_prefix_aligned.sorted.bam
```
#### Taxonomy Determination:
Create a Kraken2-Based gOTUs Database:
Mask Ribosomal RNA Genes: Use Barrnap and bedtools to mask rRNA genes in the genomes to avoid misclassification.
Kraken2 Database Creation:
```bash
kraken2-build --no-masking --add-to-library gOTU_genomes.fasta --db kraken2_gOTUs_db
```
Classify Reads by Kraken2:
```bash
kraken2 --db kraken2_gOTUs_db --report your_prefix_kraken2_report.txt \
        --output your_prefix_kraken2_output.txt your_prefix_2_trimmed.fastq.gz
```
Calculate Taxonomic Abundance with Bracken: Normalize read counts across taxonomic levels to estimate relative abundance.
```bash
bracken -d kraken2_gOTUs_db -i your_prefix_kraken2_report.txt -o your_prefix_bracken_output.txt
```
Cell Classification: Cells with >50% informative reads are considered accurately annotated; cells below this threshold are labeled with the __like flag after the species name.

### 4. Objects Generation and Process
### For Laboratory-Cultured Microbial Samples:
#### Reading and Filtering:
Follow the general scRNA-seq code and standards to read the single-bacteria gene expression matrix into R to generate Seurat objects or python to generate Scanpy objects and perform quality control of the count matrix.
#### For samples treated with antibiotics in smRandom-seq paper, use the following analysis scripts: 
CIP Treatment: `CIP_antibiotic_treatment.py`
AMP Treatment: `AMP_antibiotic_treatment.py`
GO Enrichment Analysis: `GO_enrichment.R`
#### Use the following scripts for validation performance Results in smRandom-seq paper:
Gene Expression Correlation: `Correlation_of_gene_expression.R`
Saturation Analysis: `Saturation_analysis.sh`
Technical Reproducibility: `Technical_reproducibility.py`

### For Human Feces Samples:
#### Microbiome transcriptional analysis: 
Based on major genera, Seurat (v4) performs dimensional reduction and graph-based clustering. The R script MIC-Bac/`microbiome_transcriptional_analysis.R` is an example for microbial further analysis.
#### Host-phage  relationship analysis: 
Top 20 phages per genus are compared to the UHGG genus reference genome using blastn. Phage-host correlation is validated by ≥3 kb alignment length and ≥80% genome identity. The R script MIC-Phage/`analysis_pipeline.R` is an example for host-phage transcriptonal relationship further analysis.

### For Bovine Rumen Fluid Samples:
#### Reading and Filtering:
Use `Reading_and_Filtering.R` to read the single-cell gene expression matrix into R, creates Seurat objects, and performs quality control to filters out low-quality cells and retains high-quality cells.
#### Benchmarking and Clustering:
Use `MscT_initial_matrix.R` to performs benchmarking and clustering on the Seurat objects.
#### High Metabolic Activity Cells Analysis:
Use `MscT_HMACs.R` to identifies and isolates cells with high metabolic activity (HMACs) within the dataset. Re-clusters these cells to analyze specific metabolic functions and activities.
#### Basfia succiniciproducens Re-clustering:
Use `MscT_basfia.R` to re-clusters cells of Basfia succiniciproducens to examine functional heterogeneity within this specific microbial species.
#### Pseudo-time Analysis for Basfia succiniciproducens:
Use `MscT_basfia_monocle.R` to applies Monocle to analyze the developmental trajectory of cells in Basfia succiniciproducens, exploring cellular dynamics over pseudo-time.

