# Data Analysis for Life Sciences - BENG628

This repository contains R code and workflows from the BENG628 lectures on data analysis techniques applied to life science and medicine.

## Course Description

BENG628 focuses on equipping students with the skills to analyze life sciences data using R. The lectures cover a range of topics, including data wrangling, visualization, statistical analysis, and specific applications in genomics, transcriptomics, and other areas.

## Repository Contents

* **Lecture 1: Introduction to R and Data Exploration**
    * `BENG628_Lecture1.Rproj`: RStudio project file.
    * `BENG628_Lecture1.R`: R script with code for importing data, basic data types, control structures, and initial data exploration.
    * `data/`: Folder containing the `Clinical_Data_Discovery_Cohort.csv` dataset used in the lecture.
* **Lecture 2: Exploratory Data Analysis for Life Sciences**
    * `BENG628_Lecture2.R`: R script demonstrating data wrangling with dplyr, ggplot2 visualizations, and an introduction to survival analysis.
    * `data/`: Folder containing the `Global_Health_Statistics.csv` dataset.
* **Lecture 3: Advanced Data Wrangling and Visualization**
    * `BENG628_Lecture3.R`: R script covering advanced data manipulation techniques, including reshaping, string manipulation, and handling missing data. Also includes more advanced ggplot2 visualizations.
* **Lecture 4: Data Visualization and Basic Statistical Analysis**
    * `BENG628_Lecture4.R`: R script focusing on data visualization principles and techniques using ggplot2, as well as basic statistical summaries and functions.
* **Lecture 5: RNA-seq Analysis for Rare Diseases**
    * `BENG628_Lecture5.R`: R script demonstrating RNA-seq analysis using DESeq2 for rare genetic disorders, with a case study on Facioscapulohumeral Muscular Dystrophy (FSHD). Includes differential expression analysis, volcano plot visualization, and gene search.
* **Lecture 6: Working with Genomic Data**
    * `BENG628_Lecture6.R`: R script demonstrating how to access and analyze genomic data from resources like GEO (Gene Expression Omnibus) and Ensembl.
* **Lecture 7: Advanced Genomic Data Analysis**
    * `BENG628_Lecture7.R`: R script covering more advanced genomic data analysis techniques, including differential expression analysis, methylation analysis, and integration of multiple data types.
* **Lecture 8: Genome-Wide Association Studies (GWAS)**
    * `BENG628_Lecture8.R`: R script focusing on the analysis and visualization of GWAS data, including using tools like locuszoomr.
* **output/`: Folder to store any output files generated by the scripts.

## R Packages Used

* **Core Tidyverse:**
    * dplyr: Data manipulation
    * tidyr: Data reshaping
    * ggplot2: Data visualization
    * readr: Data import
    * purrr: Functional programming
    * tibble: Modern data frames
    * stringr: String manipulation
* **Other Key Packages:**
    * survival: Survival analysis
    * survminer: Visualization for survival analysis
    * Seurat: Single-cell RNA-seq analysis
    * DESeq2: RNA-seq differential expression analysis
    * ChIPseeker: ChIP-seq data annotation
    * edgeR: RNA-seq differential expression analysis
    * TCGAbiolinks: Accessing and analyzing TCGA data
    * SummarizedExperiment: Bioconductor class for genomic data
    * here: Project-relative paths
    * ggpubr: Publication-ready plots
    * rstatix: Statistical analysis
    * pheatmap: Heatmaps
    * ComplexHeatmap: Advanced heatmaps
    * biomaRt: Accessing biological databases
    * viridis: Colorblind-friendly palettes
    * factoextra: PCA visualization
    * naniar: Missing data visualization
    * mice: Imputation of missing data
    * randomForest: Random forest modeling
    * ranger: Fast random forest implementation
    * clusterProfiler: Enrichment analysis
    * org.Hs.eg.db: Human gene annotation
    * STRINGdb: Protein-protein interaction analysis
    * FactoMineR: Multivariate analysis
    * cluster: Clustering analysis
    * MVN: Multivariate normality tests
    * lubridate: Working with dates and times
    * airway: Example RNA-seq dataset
    * GEOquery: Downloading GEO data
    * ggseqlogo: Sequence logos
    * rgl: 3D visualizations
    * rentrez: Accessing NCBI databases
    * Biostrings: Biological sequence analysis
    * msa: Multiple sequence alignment
    * ape: Phylogenetic analysis
    * ggtree: Phylogenetic tree visualization
    * phangorn: Phylogenetics
    * locuszoomr: GWAS visualization
    * AnnotationHub: Accessing genomic annotations
    * ensembldb: Working with Ensembl databases
    * plotrix: 3D pie charts (use with caution!)
    * ggrepel: Text labels for plots
    * UCSCXenaTools: Accessing TCGA data from UCSC Xena
    * missMethyl: Methylation pathway analysis
    * IlluminaHumanMethylation450kanno.ilmn12.hg19: Methylation annotation
    * IlluminaHumanMethylation450kmanifest: Methylation manifest
    * ChAMP: Additional methylation analysis tools
    * data.table: Fast data manipulation

## How to Run the Code

1. **Clone the repository:** `git clone https://github.com/oktaykaplan/Data-Analysis-for-Life-Sciences.git`
2. **Open the RStudio project files:** For each lecture, open the corresponding `.Rproj` file in RStudio. This will set the correct working directory.
3. **Install required packages:** Make sure you have all the necessary R packages installed. You can install them using `install.packages()` or `BiocManager::install()` for Bioconductor packages.
4. **Run the R scripts:** Execute the code in the `.R` files sequentially, following the instructions and comments provided.

## Data Sources

* **Lecture 1:**
    * Clinical Data Discovery Cohort: [https://www.kaggle.com/datasets/imtkaggleteam/clinical-dataset](https://www.kaggle.com/datasets/imtkaggleteam/clinical-dataset)
* **Lecture 2:**
    * Global Health Statistics Dataset: [https://www.kaggle.com/datasets/malaiarasugraj/global-health-statistics](https://www.kaggle.com/datasets/malaiarasugraj/global-health-statistics)
* **Lecture 5:**
    * data/extracted_files/: Folder for GSE174301 RNA-seq data (to be downloaded manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174301).
    * output/Lecture5_VolcanoPlot.png: Generated volcano plot highlighting key FSHD-related genes.

      ![PCA_data](https://github.com/user-attachments/assets/9e4e5003-a47c-4660-ad78-908084d785aa)

      ![Rplot](https://github.com/user-attachments/assets/cd8bb85a-87fa-411a-935e-bb996d5e827b)

* **Lecture 6:**
    * GSE37745 (Kidney diseases): [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37745](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37745)  
* **Lecture 7:**
    * GSE11121 (Breast cancer cell lines): [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121)
    * GSE183879 (Obesity and diabetes): [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183879](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183879)
    * GSE42865 (Cancer vs. normal methylation): [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42865](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42865)
    * TCGA Prostate Adenocarcinoma (PRAD): Accessed via UCSC Xena
* **Lecture 8:**
    * GWAS Catalog data (GPR75): [https://www.ebi.ac.uk/gwas/publications/34210852](https://www.ebi.ac.uk/gwas/publications/34210852)

## Notes

* The code in this repository is provided as-is, for educational purposes.
* It is recommended to run the code in the order of the lectures.
* Some scripts may require additional setup or data download.
* Please refer to the comments within the scripts for specific instructions.

## Contributing

If you find any errors or have suggestions for improvement, feel free to submit an issue or a pull request.
