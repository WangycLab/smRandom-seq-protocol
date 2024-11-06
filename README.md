# smRandom-seq Protocol Analysis Workflow

This guide provides a comprehensive, step-by-step workflow for analyzing single-bacteria total-RNA-seq data using smRandom-seq, applicable to various sample types, including laboratory-cultured bacterial samples, human feces, and bovine rumen fluid. Each command is detailed, with explanations of dependencies and optimization suggestions.

## Contents
1. [Setup and Dependencies](#1-setup-and-dependencies)
2. [Sequencing Data Pre-processing](#2-sequencing-data-pre-processing)
3. [Genome Index Creation and Read Alignment for Laboratory-Cultured Bacterial Samples](#3-genome-index-creation-and-read-alignment-for-laboratory-cultured-bacterial-samples)
4. [Selecting Cell Barcodes (For Non-Laboratory Samples)](#4-selecting-cell-barcodes-for-non-laboratory-samples)
5. [Validation of smRandom-seq Results](#5-validation-of-smrandom-seq-results)
6. [Downstream Analysis for Antibiotic-Treated Samples](#6-downstream-analysis-for-antibiotic-treated-samples)

---

### 1. Setup and Dependencies

Ensure the following dependencies are installed:
- **Python 3.8+**: For running custom scripts.
- **R 4.0+**: For statistical analysis and visualization.
- **SAMtools**: For sorting and indexing BAM files.
- **bedtools**: For genomic data manipulation.
- **STAR**: For genome indexing and read alignment.
- **featureCounts**: For exon-based gene expression counting.
- **UMI-tools**: For UMI-based read counting.

#### Installation Example:
```bash
conda install python=3.8 r-base=4.0 samtools bedtools star subread umi-tools
