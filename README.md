# smRandom-seq Protocol Analysis Workflow

This guide provides a comprehensive, step-by-step workflow for analyzing single-bacteria total-RNA-seq data using smRandom-seq, applicable to various sample types, including laboratory-cultured bacterial samples, human feces, and bovine rumen fluid.
## Contents
1. [Setup and Dependencies](#1-setup-and-dependencies)
2. [Sequencing Data Pre-processing](#2-sequencing-data-pre-processing)
3. [Count Matrix Generation](#3-count-matrix-generation)
4. [Seurat Objects Generation and Process](#4-seurat-objects-generation-and-process)
---

### 1. Setup and Dependencies

Ensure the following dependencies are installed:
Required for all samples：
- **Python 3.8+**: For running custom scripts.
- **R 4.0+**: For statistical analysis and visualization.
- **SAMtools**: For sorting and indexing BAM files.
- **bedtools**: For genomic data manipulation.
- **STAR**: For genome indexing and read alignment.
- **featureCounts**: For exon-based gene expression counting.
- **UMI-tools**: For UMI-based read counting.
Dependencies Based on Sample Type:
For Human Feces Samples:
- **MICtools**: Includes modules for taxonomic annotation (MIC-Anno), transcriptional matrix construction (MIC-Bac), and host-phage transcriptional analysis (MIC-Phage).
- **Kraken2**: For K-mer-based taxonomic classification using the UHGG microbiome genome database.
- **Bracken**: For calculating taxonomic abundance after Kraken2 classification.
- **SortMeRNA**: For filtering rRNA and tRNA reads.
- **BLAST (blastn)**: For phage-host correlation analysis.
For Bovine Rumen Fluid Samples:
- **Kraken2**: For K-mer-based taxonomic classification using the UHGG microbiome genome database.
- **Bracken**: For calculating taxonomic abundance after Kraken2 classification.

### 2. Sequencing Data Pre-processing
Process raw FastQ data with `Raw_data_preprocessing.sh` using the custom tool `anchoradp.o` to separate reads, extract barcodes and UMIs, and prepare data for downstream analyses.

About `anchoradp.o`:
Processes raw FastQ data without sequencing adapters.
Separates Read 1 and Read 2, identifies PCR adapters, and extracts key sequences.
Read 1: Extracts the 20-bp cell barcode and 8-bp UMI.
Read 2: Extracts the cDNA sequence.
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
your_prefix_1.fq.gz: Contains compressed cell barcodes and UMIs.
your_prefix_2.fq.gz: Contains compressed cDNA sequences.

### 3. Count Matrix Generation
The workflow for genome index creation and read alignment and count matrix generation depends on the sample type.
### For Laboratory-Cultured Bacterial Samples:
Use the `Count_matrix.sh` script to perform all steps of genome index creation, read alignment, and gene expression counting. This script integrates all necessary commands for STAR indexing, alignment, feature counting, and barcode selection into a single workflow.
Running `Count_matrix.sh`
Use the following command to run the script with the required input parameters:
```bash
bash Count_matrix.sh genomeDir=$genomeDir genomeFastaFile=$genomeFastaFile sjdbGTFfile=$sjdbGTFfile sample=$sample thres=$thres
```
#### Input Parameters:
genomeDir: Directory for STAR genome index files.
genomeFastaFile: Path to the genome FASTA file.
sjdbGTFfile: Path to the GTF file containing gene annotations.
sample: Prefix for the sample name (e.g., sample_2.fq.extracted and sample_1.fq.extracted).
thres: Threshold for cell count, usually between 50,00 and 50,000.
#### Output:
$sample.counts.tsv
#### Filter barcodes based on UMI and gene count levels using `selectResult.py` to create a refined count matrix:
```bash
python selectResult.py $sample.counts.tsv $thres $sample.selected
```
#### Output:

### For Human Feces Samples:
Use the `MICtools` pipeline for count matrix generation of human gut microbiome single-microbe RNA sequencing. MICtools provides modules for taxonomic annotation (`MIC-Anno`), transcriptional matrix construction (`MIC-Bac`), and host-phage transcriptional relationship analysis(`MIC-Phage`).
#### Barcodes selection: Barcode and gene count scatter plots are analyzed to set a threshold for barcode selection. Barcodes beyond the inflection point are filtered out.
#### Use `MIC-Anno` to process taxonomic annotation for each microbe: MIC-Anno uses a K-mer-based taxonomic classification strategy with Kraken2 and the UHGG gut microbiome genome database. 
```bash
MICtools anno --module pipeline -s your_prefix_2.fq.gz -r [kraken_ref_filepath] -p [your_prefix]
```
#### Output:
your_prefix.barcode_count.txt: Main single microbe taxonomic annotation result and the meanings of each column.
your_prefix.genus_info.txt: Genus info of each microbe.
your_prefix.species_info.txt: Species info of each microbe.
your_prefix_genus.pdf: Pieplot of sample genus composition.
your_prefix_species.pdf: Pieplot of sample specie composition.
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

### 4. Seurat Objects Generation and Process
### For Laboratory-Cultured Bacterial Samples:
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

