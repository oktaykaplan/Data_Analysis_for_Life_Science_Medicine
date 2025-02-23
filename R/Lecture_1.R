#Download R from https://cran.r-project.org/
#Download RStudio from https://posit.co/downloads/

#Installing Key Packages
#tidyverse - Collection of data science packages, including dplyr, tidyr, readr, ggplot2, purrr, tible and stringr
#ggplot2 - Data visualization
#dplyr - Data manipulation
#readr - Fast CSV file reading
#readxl - Excel file reading
#jsonlite - Handling JSON data
#survival - Survival analysis tools
#survminer - Visualization of survival analysis
#Seurat - Single-cell RNA-seq analysis
#DESeq2 - Differential expression analysis for RNA-Seq
#ChIPseeker - ChIP-seq data annotation
#edgeR - RNA-seq differential expression analysis
#TCGAbiolinks: A comprehensive package for accessing and analyzing data from The Cancer Genome Atlas (TCGA).
#SummarizedExperiment: A Bioconductor class for storing and working with high-throughput genomic data, often used with RNA-seq and other assays.


# Install Packages (Do this ONCE - preferably at the start of your project)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Seurat", "DESeq2", "ChIPseeker", "edgeR", "TxDb.Hsapiens.UCSC.hg38.knownGene", "TCGAbiolinks", "SummarizedExperiment"))

install.packages(c("tidyverse", "survival", "survminer", "here"), dependencies = TRUE) # Install here package



#After Installation, 
#Once all installations are complete, you should be able to load them
# Load Libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(jsonlite)
library(survival)
library(survminer)
library(Seurat)
library(DESeq2)
library(ChIPseeker)
library(edgeR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(here)
library(TCGAbiolinks)
library(SummarizedExperiment)


#Importing and Exporting Biomedical Data
#Reading Data from Different Formats
#Setting Working Directory (or using here())
#setwd("/path/to/your/project")  # OR use here() for project-relative paths

# Setting Up Working Directory
here::i_am("BENG628_Lecture1.Rproj") # Create a project file first!


# Example using here()
biomedical_data_path <- here("data", "Clinical_Data_Discovery_Cohort.csv") # Example with subfolder

# Importing Data
biomedical_data <- read_csv(biomedical_data_path)

# Exporting Data
write.csv(biomedical_data, here("output", "output_data.csv"), row.names = FALSE)
write.xlsx(data, "output_data.xlsx")


#Variables and Data Types
x <- 10  # Numeric
y <- "GeneX"  # Character
z <- TRUE  # Logical
vec <- c(1, 2, 3, 4)  # Vector

#Control Structures
for (i in 1:5) {
  print(paste("Patient ID:", i))
}
if (gene_expression > 2.0) {
  print("Upregulated gene")
} else {
  print("Normal expression")
}


# --- Variables and Data Types (Expanded) ---

age <- 30L # Integer (add L suffix)
gene_name <- "TP53"
is_cancer <- TRUE
height <- 175.5 # Double (floating point)
cancer_stage <- factor(c("I", "II", "III"), levels = c("I", "II", "III"), ordered = TRUE) # Ordered factor
patient_list <- list(age = age, gene = gene_name, is_cancer = is_cancer) # List
patient_df <- data.frame(age = c(30, 40), gene = c("TP53", "BRCA1"), stage = cancer_stage[1:2]) # Data frame

# --- Control Structures (More Examples) ---

# For loop with index and value
genes <- c("TP53", "BRCA1", "EGFR")
for (i in seq_along(genes)) { # Iterate by index
  gene <- genes[i]
  print(paste("Gene", i, ":", gene))
}

# If/else with multiple conditions
gene_expression <- 1.8
if (gene_expression > 2.0) {
  print("High expression")
} else if (gene_expression > 1.0 && gene_expression <= 2.0) { # Added condition
  print("Moderate expression")
} else {
  print("Low expression")
}

#Hands-on: Loading a Biomedical Dataset
#Download Sample Dataset
#https://www.kaggle.com/datasets/imtkaggleteam/clinical-dataset?resource=download
#Loading Data in R
library(readr)
biomedical_data <- read_csv("Clinical Data_Discovery_Cohort.csv")
head(biomedical_data)
summary(biomedical_data)
str(biomedical_data)
table(biomedical_data$`Dead or Alive`)


#Data Analysis & Visualization with ggplot2

#Exploratory Data Analysis (EDA)


# Check for missing values
sum(is.na(biomedical_data))

# Convert date columns to Date format
biomedical_data$`Specimen date` <- as.Date(biomedical_data$`Specimen date`, format="%m/%d/%Y")
biomedical_data$`Date of Death` <- as.Date(biomedical_data$`Date of Death`, format="%m/%d/%Y")
biomedical_data$`Date of Last Follow Up` <- as.Date(biomedical_data$`Date of Last Follow Up`, format="%m/%d/%Y")

# Summary statistics for survival time
biomedical_data %>% 
  group_by(`Dead or Alive`) %>% 
  summarize(
    Mean_Time = mean(Time, na.rm=TRUE),
    SD_Time = sd(Time, na.rm=TRUE)
  )

#Data Visualization: Survival Time by Status

library(ggplot2)
ggplot(biomedical_data, aes(x=`Dead or Alive`, y=Time, fill=`Dead or Alive`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title="Survival Time Across Patient Status",
       x="Patient Status",
       y="Survival Time (Days)")

#Distribution of Stages
ggplot(biomedical_data, aes(x=Stage, fill=Stage)) +
  geom_bar() +
  theme_minimal() +
  labs(title="Distribution of Cancer Stages",
       x="Stage",
       y="Count") +
  coord_flip()

#Introduction to Survival Analysis
library(survival)
library(survminer)
surv_object <- Surv(biomedical_data$Time, biomedical_data$`Dead or Alive` == "Dead")
fit <- survfit(surv_object ~ 1)
ggsurvplot(fit, data=biomedical_data, risk.table=TRUE, pval=TRUE, conf.int=TRUE,
           title="Kaplan-Meier Survival Curve")

# --- Gene Expression Analysis, Single-cell RNA-seq, ChIP-seq ---

# Gene Expression Analysis with DESeq2
count_data <- matrix(rnbinom(1000, mu = 10, size = 1), ncol = 10)
col_data <- data.frame(condition = factor(rep(c("Control", "Treated"), each = 5)))
rownames(col_data) <- colnames(count_data)
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
summary(res)


# Example: TCGAbiolinks (Brief demonstration)
# What TCGA is and  why it is important)
# Example: Download a small TCGA dataset (demonstration only)
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts" # Wrkflow type
)

GDCdownload(query) # Download the data
data <- GDCprepare(query) # Prepare the data as a SummarizedExperiment



#Single-cell RNA-seq Analysis
library(Seurat)
sce <- CreateSeuratObject(counts = matrix(rpois(200, lambda = 5), ncol = 10))
VlnPlot(sce, features = "nFeature_RNA")


#ChIP-seq Data Analysis & Correlation with Gene Expression
# Load (In each R session)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Now you can use the packages
# ChIP-seq Data Analysis
peaks <- readPeakFile("chipseq_peaks.bed")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(peaks, TxDb=txdb)
plotAnnoBar(peakAnno)

# Save Analysis Outputs
saveRDS(dds, file = here("output", "deseq2_results.rds"))
saveRDS(peakAnno, file = here("output", "chipseq_annotations.rds"))

#PCA (Principal Component Analysis) for Dimension Reduction
#Useful for exploring high-dimensional gene expression data.
# Perform PCA on gene expression matrix
gene_expr_matrix <- matrix(rnorm(1000), ncol=10)
pca_result <- prcomp(gene_expr_matrix, center = TRUE, scale. = TRUE)
biplot(pca_result)

#Differential Methylation Analysis
# Bioconductor package for methylation analysis
BiocManager::install("ChAMP")
library(ChAMP)

# --- Software Versions ---
sessionInfo() # Show R and package versions


# BENG 628 - Lecture 1 - Take-Home Task - Sample Solution
# Load required libraries
library(survival)
library(ggplot2)
library(survminer)
library(dplyr)

# Load dataset
ls("package:survival")
data("lung", package = "survival")

#Confirm the dataset exists in the package
datasets <- data(package = "survival")$results
datasets

data("veteran", package = "survival")
str(veteran)

data("aml", package = "survival")

# Inspect the structure of the 'aml' dataset
str(aml)

# Check the first few rows of the dataset
head(aml)

# Check for missing values
colSums(is.na(aml))

# Perform survival analysis using this dataset
surv_obj <- Surv(aml$time, aml$status)
fit <- survfit(surv_obj ~ x, data = aml)

# Summary of the survival fit
summary(fit)

# Check for missing values in 'x' (treatment) and impute with the most frequent category (since it's a factor)
aml$x[is.na(aml$x)] <- names(sort(table(aml$x), decreasing = TRUE))[1]

# Kaplan-Meier Survival Curve based on treatment (x)
surv_obj <- Surv(aml$time, aml$status == 1)
fit <- survfit(surv_obj ~ x, data = aml)
ggsurvplot(fit, data = aml, risk.table = TRUE, pval = TRUE, title = "Survival by Treatment")



# 1- Assess the Differences in Survival Between Groups

# Perform log-rank test to assess significance
#This will provide the chi-square statistic and p-value, 
#allowing you to evaluate whether the survival distributions between the two treatment groups are significantly different.
surv_diff <- survdiff(surv_obj ~ x, data = aml)
print(surv_diff)


#2. Cox Proportional Hazards Model
#If you want to investigate the relationship between survival and other variables (such as treatment, age, etc.), you can use the Cox proportional hazards model to estimate hazard ratios and the effect of different covariates on survival.
cox_model <- coxph(surv_obj ~ x, data = aml)
summary(cox_model)


#3. Check for Proportional Hazards Assumption
#The Cox model assumes that the hazards are proportional over time. You can check if this assumption holds by using Schoenfeld residuals or testing for interactions with time.
cox_zph <- cox.zph(cox_model)
print(cox_zph)

ggsurvplot(fit, data = aml, risk.table = TRUE, pval = TRUE, 
           conf.int = TRUE, title = "Survival by Treatment", 
           palette = c("#E69F00", "#56B4E9"))


