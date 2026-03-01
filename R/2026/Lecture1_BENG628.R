# ============================================================================
# BENG 628 - Introduction to R for Biomedical Data Analysis
# Comprehensive Lecture with Extensive Examples
# ============================================================================
# Author: Dr. Oktay I. Kaplan
# Date: 2026
# Description: Complete guide to R programming for biomedical research
# ============================================================================

# ============================================================================
# SECTION 1: GETTING STARTED WITH R AND RSTUDIO
# ============================================================================

# Download R from https://cran.r-project.org/
# Download RStudio from https://posit.co/downloads/

# Why R for Biomedical Research?
# - Specialized packages for genomics, proteomics, clinical data
# - Excellent visualization capabilities
# - Strong statistical analysis tools
# - Large bioinformatics community

# ============================================================================
# SECTION 2: INSTALLING AND LOADING PACKAGES
# ============================================================================

# ----------------------------------------------------------------------------
# 2.1 Installing Base Packages
# ----------------------------------------------------------------------------

# Install tidyverse (collection of data science packages)
# Only run once!
install.packages("tidyverse")  # Includes dplyr, ggplot2, readr, tidyr, etc.

# Install individual packages
install.packages(c(
  "ggplot2",      # Data visualization
  "dplyr",        # Data manipulation
  "readr",        # Fast CSV reading
  "readxl",       # Excel file reading
  "jsonlite",     # JSON data handling
  "survival",     # Survival analysis
  "survminer",    # Survival visualization
  "here"          # Project-relative paths
), dependencies = TRUE)

# ----------------------------------------------------------------------------
# 2.2 Installing Bioconductor Packages
# ----------------------------------------------------------------------------

# Bioconductor is a repository for bioinformatics packages
# Install BiocManager first (one-time installation)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# Install bioinformatics packages
BiocManager::install(c(
  "Seurat",                              # Single-cell RNA-seq
  "DESeq2",                              # Differential expression
  "ChIPseeker",                          # ChIP-seq annotation
  "edgeR",                               # RNA-seq DE analysis
  "TxDb.Hsapiens.UCSC.hg38.knownGene",  # Human genome annotation
  "TCGAbiolinks",                        # TCGA data access
  "SummarizedExperiment"                 # Genomic data storage
))

# ----------------------------------------------------------------------------
# 2.3 Loading Libraries
# ----------------------------------------------------------------------------

# Load libraries at the start of each R session
library(tidyverse)
library(ggplot2)
library(dplyr)
library(readr)
library(readxl)
library(jsonlite)
library(survival)
library(survminer)
library(here)

# Check loaded packages and versions
sessionInfo()

# ============================================================================
# SECTION 3: BASIC R PROGRAMMING CONCEPTS
# ============================================================================

# ----------------------------------------------------------------------------
# 3.1 Variables and Data Types
# ----------------------------------------------------------------------------

# --- NUMERIC TYPES ---

# Integer (use L suffix)
age <- 30L
patient_count <- 150L
print(paste("Age:", age))
print(paste("Class:", class(age)))

# Double (floating point)
gene_expression <- 2.45
pvalue <- 0.0001
temperature <- 37.5

print(paste("Expression level:", gene_expression))
print(paste("Class:", class(gene_expression)))

# --- CHARACTER (STRING) TYPES ---

gene_name <- "TP53"
patient_id <- "PT001"
diagnosis <- "Breast Cancer"

# String operations
full_label <- paste("Gene:", gene_name)
print(full_label)

upper_gene <- toupper(gene_name)
print(paste("Uppercase:", upper_gene))

lower_gene <- tolower(gene_name)
print(paste("Lowercase:", lower_gene))

# String concatenation
patient_info <- paste(patient_id, "-", diagnosis)
print(patient_info)

# --- LOGICAL (BOOLEAN) TYPES ---

is_cancer <- TRUE
is_healthy <- FALSE
has_mutation <- TRUE

# Logical operations
print(paste("Cancer AND Mutation:", is_cancer & has_mutation))  # AND
print(paste("Cancer OR Healthy:", is_cancer | is_healthy))       # OR
print(paste("NOT Cancer:", !is_cancer))                          # NOT

# Comparison operators
print(paste("5 > 3:", 5 > 3))
print(paste("10 == 10:", 10 == 10))
print(paste("7 != 8:", 7 != 8))

# ----------------------------------------------------------------------------
# 3.2 Vectors
# ----------------------------------------------------------------------------

# Vectors are one-dimensional arrays of the same data type

# Numeric vector
expression_levels <- c(1.2, 2.5, 3.1, 1.8, 2.9)
print("Expression levels:")
print(expression_levels)

# Character vector
gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")
print("Gene list:")
print(gene_list)

# Logical vector
mutation_status <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
print("Mutation status:")
print(mutation_status)

# Vector operations
print(paste("Mean expression:", mean(expression_levels)))
print(paste("Sum of mutations:", sum(mutation_status)))  # TRUE = 1, FALSE = 0
print(paste("Max expression:", max(expression_levels)))
print(paste("Min expression:", min(expression_levels)))
print(paste("Standard deviation:", sd(expression_levels)))

# Vector indexing (R starts at 1, not 0!)
print(paste("First gene:", gene_list[1]))
print(paste("Third gene:", gene_list[3]))
print(paste("Last gene:", gene_list[length(gene_list)]))

# Multiple elements
print("First and third genes:")
print(gene_list[c(1, 3)])

print("First three genes:")
print(gene_list[1:3])
print(gene_list[1:10])

print(gene_list[c(751, 1001)])

# Conditional indexing
high_expression <- expression_levels[expression_levels > 2.0]
print("High expression values:")
print(high_expression)

mutated_genes <- gene_list[mutation_status]
print("Genes with mutations:")
print(mutated_genes)

# Vector arithmetic
control <- c(1.0, 1.2, 0.9, 1.1, 1.0)
treated <- c(2.5, 3.1, 0.8, 4.2, 1.9)
fold_change <- treated / control
print("Fold changes:")
print(fold_change)

# ----------------------------------------------------------------------------
# 3.3 Factors
# ----------------------------------------------------------------------------

# Factors represent categorical data

# Unordered factor
cancer_type <- factor(c("Breast", "Lung", "Colon", "Breast", "Lung"))
print("Cancer types:")
print(cancer_type)
print("Levels:")
print(levels(cancer_type))

# Frequency table
print("Cancer type frequencies:")
print(table(cancer_type))

# Ordered factor (for stages, grades)
cancer_stage <- factor(
  c("I", "II", "III", "II", "IV", "I", "III"),
  levels = c("I", "II", "III", "IV"),
  ordered = TRUE
)

print("Cancer stages:")
print(cancer_stage)
print(paste("Is ordered:", is.ordered(cancer_stage)))

# Comparing ordered factors
print(paste("Stage III > Stage I:", cancer_stage[3] > cancer_stage[1]))

# ----------------------------------------------------------------------------
# 3.4 Lists
# ----------------------------------------------------------------------------

# Lists can contain different data types

patient <- list(
  id = "PT001",
  age = 45,
  genes = c("TP53", "BRCA1", "EGFR"),
  has_cancer = TRUE,
  expression = c(2.1, 3.4, 1.8),
  measurements = list(
    height = 170,
    weight = 68
  )
)

print("Patient list structure:")
str(patient)

# Access list elements
print(paste("Patient ID:", patient$id))
print(paste("Patient age:", patient$age))
print("Patient genes:")
print(patient$genes)
print(paste("Has cancer:", patient$has_cancer))

# Alternative access methods
print(patient[["age"]])
print(patient[[2]])  # Access by position

# Access nested list
print(paste("Height:", patient$measurements$height))

# ----------------------------------------------------------------------------
# 3.5 Data Frames
# ----------------------------------------------------------------------------

# Data frames are tables where each column can be a different type

clinical_data <- data.frame(
  patient_id = c("PT001", "PT002", "PT003", "PT004", "PT005"),
  age = c(45, 52, 38, 61, 44),
  sex = c("F", "M", "F", "M", "F"),
  gene = c("TP53", "BRCA1", "EGFR", "TP53", "MYC"),
  stage = factor(c("II", "III", "I", "IV", "II")),
  survival_days = c(365, 180, 730, 90, 450),
  status = c("Alive", "Dead", "Alive", "Dead", "Alive")
)

class_2026_BENG628 <- data.frame(
  student_name = c("Rabia", "Yaren", "Elif", "Aleyna", "Elif_Nur", "Seda", "Nida", "Lale", "Julide" ," Kodcu_Burak", "Mahmut"),
  age = c(18, 19, 18, 20, 18, 19, 18, 20, 18, 19, 18),
  sex = c("F", "F", "F", "F", "F", "F", "F", "F", "F", "M", "M"),
  gene = c("TP53", "BRCA1", "EGFR", "MDM2", "MYC", "WDR31", "BBS1", "ARL13B", "NRAS", "BURAK", "BBIP10"),
  status = c("Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive", "Alive")
)


class_2026_BENG628

print(class_2026_BENG628)
clinical_data1 <- data.frame(
  Rabia = c("PT001", "PT002", "PT003", "PT004", "PT005"),
  age = c(45, 52, 38, 61, 44),
  sex = c("F", "M", "F", "M", "F"),
  Turkiye = c("TP53", "BRCA1", "EGFR", "TP53", "MYC"),
  stage = factor(c("II", "III", "I", "IV", "II")),
  survival_days = c(365, 180, 730, 90, 450),
  status = c("Alive", "Dead", "Alive", "Dead", "Alive")
)

print(clinical_data1)
print("Clinical data:")
print(clinical_data)

# View structure
str(clinical_data)


colnames(clinical_data)

# Dimensions
print(paste("Number of rows:", nrow(clinical_data)))
print(paste("Number of columns:", ncol(clinical_data)))
print("Dimensions (rows, columns):")
print(dim(clinical_data))

# Column names
print("Column names:")
print(names(clinical_data))

# Access columns
print("Ages:")
print(clinical_data$age)
print(clinical_data[, "age"])
print(clinical_data[, 2])

# Access rows
print("First patient:")
print(clinical_data[1, ])

print("First three patients:")
print(clinical_data[1:3, ])

# Access specific cell
print(paste("Second patient, third column:", clinical_data[2, 3]))

# Subset data frame
female_patients <- clinical_data[clinical_data$sex == "F", ]
print("Female patients:")
print(female_patients)

high_risk <- clinical_data[clinical_data$stage %in% c("III", "IV"), ]
print("High risk patients:")
print(high_risk)

# ----------------------------------------------------------------------------
# 3.6 Matrices
# ----------------------------------------------------------------------------

# Matrices are 2D arrays of the same data type

# Create a gene expression matrix
expr_matrix <- matrix(
  c(1.2, 2.3, 1.8, 3.1, 2.5, 1.9, 2.8, 2.1, 3.5, 1.7, 2.9, 2.4),
  nrow = 3,
  ncol = 4,
  byrow = TRUE
)

# Add row and column names
rownames(expr_matrix) <- c("Gene1", "Gene2", "Gene3")
colnames(expr_matrix) <- c("Sample1", "Sample2", "Sample3", "Sample4")

print("Expression matrix:")
print(expr_matrix)

# Access elements
print(paste("Gene1, Sample1:", expr_matrix[1, 1]))
print("Gene2, all samples:")
print(expr_matrix[2, ])
print("All genes, Sample3:")
print(expr_matrix[, 3])

# Matrix operations
print("Column means:")
print(colMeans(expr_matrix))
print("Row means:")
print(rowMeans(expr_matrix))

# ============================================================================
# SECTION 4: CONTROL STRUCTURES
# ============================================================================

# ----------------------------------------------------------------------------
# 4.1 For Loops
# ----------------------------------------------------------------------------

print("=== FOR LOOP EXAMPLES ===")

# Basic for loop
genes <- c("TP53", "BRCA1", "EGFR", "MYC", "KRAS")

print("Analyzing genes:")
for (gene in genes) {
  print(paste("Analyzing gene:", gene))
}

# Loop with index
print("\nGenes with index:")
for (i in 1:length(genes)) {
  print(paste("Gene", i, ":", genes[i]))
}

# Better way using seq_along
print("\nUsing seq_along:")
for (i in seq_along(genes)) {
  print(paste("Position", i, "- Gene:", genes[i]))
}

# Practical example: Calculate fold changes
control_expression <- c(1.0, 1.2, 0.9, 1.1, 1.0)
treated_expression <- c(2.5, 3.1, 0.8, 4.2, 1.9)

print("\nCalculating fold changes:")
fold_changes <- numeric(length(genes))

for (i in 1:length(genes)) {
  fold_changes[i] <- treated_expression[i] / control_expression[i]
  print(paste(genes[i], "fold change:", round(fold_changes[i], 2)))
}

# Nested loops
print("\nNested loop example - Gene expression across conditions:")
conditions <- c("Control", "Treated", "Recovery")
for (gene in c("TP53", "BRCA1")) {
  for (condition in conditions) {
    value <- round(runif(1, 0, 5), 2)
    print(paste(gene, "-", condition, ":", value))
  }
}

# ----------------------------------------------------------------------------
# 4.2 If-Else Statements
# ----------------------------------------------------------------------------

print("\n=== IF-ELSE EXAMPLES ===")

# Simple if-else
gene_expression <- 2.8

if (gene_expression > 2.0) {
  print("Result: Upregulated gene")
} else {
  print("Result: Normal or downregulated")
}

# Multiple conditions with else if
pvalue <- 0.03
fold_change <- 2.5

print("\nDifferential expression analysis:")
if (pvalue < 0.05 & fold_change > 2.0) {
  print("Result: Significantly upregulated")
} else if (pvalue < 0.05 & fold_change < 0.5) {
  print("Result: Significantly downregulated")
} else if (pvalue < 0.05) {
  print("Result: Significant but small change")
} else {
  print("Result: Not significant")
}

# Classify gene expression levels
expression <- 1.8

print("\nExpression classification:")
if (expression > 3.0) {
  category <- "Very High"
} else if (expression > 2.0) {
  category <- "High"
} else if (expression > 1.0) {
  category <- "Moderate"
} else {
  category <- "Low"
}
print(paste("Expression category:", category))

# Cancer staging example
tumor_size <- 3.5  # cm
lymph_nodes <- 2
metastasis <- FALSE

print("\nCancer TNM staging:")

# T stage
if (tumor_size <= 2) {
  t_stage <- "T1"
} else if (tumor_size <= 5) {
  t_stage <- "T2"
} else {
  t_stage <- "T3"
}

# N stage
if (lymph_nodes == 0) {
  n_stage <- "N0"
} else if (lymph_nodes <= 3) {
  n_stage <- "N1"
} else {
  n_stage <- "N2"
}

# M stage
if (metastasis) {
  m_stage <- "M1"
} else {
  m_stage <- "M0"
}

print(paste("TNM Classification:", t_stage, n_stage, m_stage))

# ifelse() function (vectorized)
print("\nVectorized ifelse:")
ages <- c(25, 45, 65, 72, 38)
age_groups <- ifelse(ages >= 65, "Elderly", ifelse(ages >= 45, "Middle", "Young"))
print(data.frame(age = ages, group = age_groups))

# ----------------------------------------------------------------------------
# 4.3 While Loops
# ----------------------------------------------------------------------------

print("\n=== WHILE LOOP EXAMPLES ===")

# Cell division simulation
cell_count <- 1
generation <- 0
target_cells <- 1000

print("Cell division simulation:")
while (cell_count < target_cells) {
  cell_count <- cell_count * 2
  generation <- generation + 1
  if (generation <= 5 || cell_count >= target_cells) {
    print(paste("Generation", generation, ":", cell_count, "cells"))
  }
}
print(paste("Final generation:", generation))

# Drug concentration decay
print("\nDrug concentration decay:")
concentration <- 100  # mg/L
time <- 0
threshold <- 10

while (concentration > threshold) {
  concentration <- concentration * 0.8  # 20% decay per hour
  time <- time + 1
  if (time <= 5 || concentration <= threshold * 1.5) {
    print(paste("Hour", time, ":", round(concentration, 2), "mg/L"))
  }
}

# ----------------------------------------------------------------------------
# 4.4 Apply Family Functions
# ----------------------------------------------------------------------------

print("\n=== APPLY FAMILY EXAMPLES ===")

# Create sample data
expr_data <- matrix(rnorm(20, mean = 5, sd = 2), nrow = 4, ncol = 5)
rownames(expr_data) <- c("Gene1", "Gene2", "Gene3", "Gene4")
colnames(expr_data) <- c("S1", "S2", "S3", "S4", "S5")

print("Expression data:")
print(round(expr_data, 2))

# apply() - for matrices/arrays
print("\nRow means (gene expression averages):")
row_means <- apply(expr_data, 1, mean)
print(round(row_means, 2))

print("\nColumn means (sample averages):")
col_means <- apply(expr_data, 2, mean)
print(round(col_means, 2))

# lapply() - returns a list
gene_names <- c("TP53", "BRCA1", "EGFR")
print("\nlapply example - convert to lowercase:")
lower_genes <- lapply(gene_names, tolower)
print(lower_genes)

# sapply() - returns a vector/matrix
print("\nsapply example - string lengths:")
name_lengths <- sapply(gene_names, nchar)
print(name_lengths)

# ============================================================================
# SECTION 5: WORKING DIRECTORY AND FILE MANAGEMENT
# ============================================================================

print("\n=== WORKING DIRECTORY ===")

# Check current working directory
current_dir <- getwd()
print(paste("Current directory:", current_dir))

# Set working directory (use your actual path)
# setwd("/path/to/your/project")

# Using the 'here' package (RECOMMENDED)
# First, create an .Rproj file in your project folder
# Then use here::i_am() to initialize

# library(here)
# here::i_am("BENG628_Project.Rproj")

# Now you can use relative paths
# data_path <- here("data", "clinical_data.csv")
# output_path <- here("output", "results.csv")

# List files in directory
# print("Files in current directory:")
# print(list.files())

# ============================================================================
# SECTION 6: DATA IMPORT AND EXPORT
# ============================================================================

# ----------------------------------------------------------------------------
# 6.1 Creating Sample Data for Examples
# ----------------------------------------------------------------------------

# Create sample clinical data
sample_clinical <- data.frame(
  PatientID = paste0("PT", sprintf("%03d", 1:20)),
  Age = sample(30:80, 20, replace = TRUE),
  Sex = sample(c("M", "F"), 20, replace = TRUE),
  Stage = sample(c("I", "II", "III", "IV"), 20, replace = TRUE),
  Time = sample(30:1000, 20, replace = TRUE),
  Status = sample(c("Alive", "Dead"), 20, replace = TRUE),
  Gene = sample(c("TP53", "BRCA1", "EGFR", "KRAS"), 20, replace = TRUE)
)

print("Sample clinical data:")
print(head(sample_clinical))

# ----------------------------------------------------------------------------
# 6.2 Exporting Data
# ----------------------------------------------------------------------------

# Create output directory if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
}

# Export as CSV
write.csv(sample_clinical, "output/sample_clinical_data.csv", row.names = FALSE)
print("Exported to: output/sample_clinical_data.csv")

# Export using readr (faster, better)
write_csv(sample_clinical, "output/sample_clinical_readr.csv")
print("Exported using readr")

# ----------------------------------------------------------------------------
# 6.3 Importing Data
# ----------------------------------------------------------------------------

# Import CSV using base R
# clinical_data <- read.csv("data/Clinical_Data_Discovery_Cohort.csv")

# Import CSV using readr (faster, more consistent)
# clinical_data <- read_csv("data/Clinical_Data_Discovery_Cohort.csv")

# Read our sample data back
imported_data <- read_csv("output/sample_clinical_readr.csv")
print("Imported data:")
print(head(imported_data))

# ----------------------------------------------------------------------------
# 6.4 Working with Excel Files
# ----------------------------------------------------------------------------

# Install and load readxl if needed
# install.packages("readxl")
# library(readxl)

# Read Excel file
# gene_data <- read_excel("data/gene_expression.xlsx")

# Read specific sheet
# gene_data <- read_excel("data/gene_expression.xlsx", sheet = "Results")

# Read specific range
# gene_data <- read_excel("data/gene_expression.xlsx", range = "A1:D100")

# ----------------------------------------------------------------------------
# 6.5 Working with RDS Files
# ----------------------------------------------------------------------------

# Save single R object (preserves structure)
saveRDS(sample_clinical, "output/clinical_data.rds")
print("Saved as RDS file")

# Load RDS file
loaded_clinical <- readRDS("output/clinical_data.rds")
print("Loaded from RDS:")
print(head(loaded_clinical))

# Save multiple objects
analysis_results <- list(mean_age = mean(sample_clinical$Age))
save(sample_clinical, analysis_results, file = "output/all_data.RData")
print("Saved multiple objects")

# Load multiple objects
# load("output/all_data.RData")

# ============================================================================
# SECTION 7: EXPLORATORY DATA ANALYSIS (EDA)
# ============================================================================

print("\n=== EXPLORATORY DATA ANALYSIS ===")

# Use our sample data
biomedical_data <- sample_clinical

# ----------------------------------------------------------------------------
# 7.1 Initial Data Inspection
# ----------------------------------------------------------------------------

print("\n--- Initial Inspection ---")

# Dimensions
print(paste("Rows:", nrow(biomedical_data)))
print(paste("Columns:", ncol(biomedical_data)))
print("Dimensions:")
print(dim(biomedical_data))

# Column names
print("Column names:")
print(names(biomedical_data))

# Structure
print("\nData structure:")
str(biomedical_data)

# First few rows
print("\nFirst 6 rows:")
print(head(biomedical_data))

# Last few rows
print("\nLast 6 rows:")
print(tail(biomedical_data))

# Summary statistics
print("\nSummary statistics:")
print(summary(biomedical_data))

# ----------------------------------------------------------------------------
# 7.2 Checking Data Quality
# ----------------------------------------------------------------------------

print("\n--- Data Quality Checks ---")

# Check for missing values
total_na <- sum(is.na(biomedical_data))
print(paste("Total missing values:", total_na))

# Missing values per column
print("\nMissing values per column:")
print(colSums(is.na(biomedical_data)))

# Percentage of missing values
print("\nPercentage missing per column:")
missing_pct <- colSums(is.na(biomedical_data)) / nrow(biomedical_data) * 100
print(round(missing_pct, 2))

# Check for duplicates
duplicates <- sum(duplicated(biomedical_data))
print(paste("\nNumber of duplicate rows:", duplicates))

# Check for unique values in categorical columns
print("\nUnique values in Status:")
print(unique(biomedical_data$Status))

print("\nUnique values in Stage:")
print(unique(biomedical_data$Stage))

# ----------------------------------------------------------------------------
# 7.3 Data Type Conversions
# ----------------------------------------------------------------------------

print("\n--- Data Type Conversions ---")

# Convert to factor
biomedical_data$Stage <- factor(
  biomedical_data$Stage,
  levels = c("I", "II", "III", "IV"),
  ordered = TRUE
)

biomedical_data$Status <- factor(biomedical_data$Status)
biomedical_data$Sex <- factor(biomedical_data$Sex)
biomedical_data$Gene <- factor(biomedical_data$Gene)

print("Converted to factors:")
print(str(biomedical_data))

# ----------------------------------------------------------------------------
# 7.4 Summary Statistics by Groups
# ----------------------------------------------------------------------------

print("\n--- Summary Statistics by Groups ---")

# Using dplyr
library(dplyr)

# Summary by Status
print("\nStatistics by Status:")
status_summary <- biomedical_data %>%
  group_by(Status) %>%
  summarize(
    Count = n(),
    Mean_Age = mean(Age, na.rm = TRUE),
    SD_Age = sd(Age, na.rm = TRUE),
    Mean_Time = mean(Time, na.rm = TRUE),
    Median_Time = median(Time, na.rm = TRUE)
  )
print(status_summary)

# Summary by Stage
print("\nStatistics by Stage:")
stage_summary <- biomedical_data %>%
  group_by(Stage) %>%
  summarize(
    Count = n(),
    Mean_Age = mean(Age),
    Mean_Survival = mean(Time)
  )
print(stage_summary)

# Multiple groupings
print("\nStatistics by Status and Sex:")
multi_summary <- biomedical_data %>%
  group_by(Status, Sex) %>%
  summarize(
    Count = n(),
    Avg_Survival = mean(Time)
  )
print(multi_summary)

# Frequency tables
print("\nFrequency table - Status:")
print(table(biomedical_data$Status))

print("\nFrequency table - Stage:")
print(table(biomedical_data$Stage))

print("\nCross-tabulation - Status by Sex:")
print(table(biomedical_data$Status, biomedical_data$Sex))

# Proportions
print("\nProportions - Stage:")
print(prop.table(table(biomedical_data$Stage)))

# ----------------------------------------------------------------------------
# 7.5 Data Manipulation with dplyr
# ----------------------------------------------------------------------------

print("\n--- Data Manipulation Examples ---")

# Filter rows
print("\nPatients with Stage III or IV:")
high_risk <- biomedical_data %>%
  filter(Stage %in% c("III", "IV"))
print(high_risk)

print("\nElderly patients (Age >= 65):")
elderly <- biomedical_data %>%
  filter(Age >= 65)
print(head(elderly))

# Select columns
print("\nSelect specific columns:")
minimal_data <- biomedical_data %>%
  select(PatientID, Age, Stage, Time, Status)
print(head(minimal_data))

# Create new columns
print("\nAdd calculated columns:")
biomedical_data <- biomedical_data %>%
  mutate(
    Age_Group = case_when(
      Age < 50 ~ "Young",
      Age < 70 ~ "Middle",
      Age >= 70 ~ "Elderly"
    ),
    High_Risk = ifelse(Stage %in% c("III", "IV"), "Yes", "No"),
    Survival_Years = round(Time / 365, 2)
  )
print(head(biomedical_data))

# Arrange (sort) data
print("\nSort by survival time (descending):")
sorted_data <- biomedical_data %>%
  arrange(desc(Time))
print(head(sorted_data))

# Chain multiple operations
print("\nChaining operations:")
analysis_result <- biomedical_data %>%
  filter(!is.na(Time)) %>%
  group_by(Stage, Sex) %>%
  summarize(
    n = n(),
    mean_survival = mean(Time),
    median_survival = median(Time),
    sd_survival = sd(Time)
  ) %>%
  arrange(Stage, Sex)
print(analysis_result)

# ============================================================================
# SECTION 8: DATA VISUALIZATION WITH GGPLOT2
# ============================================================================

print("\n=== DATA VISUALIZATION ===")

library(ggplot2)

# ----------------------------------------------------------------------------
# 8.1 The Grammar of Graphics
# ----------------------------------------------------------------------------

# ggplot2 uses layers:
# 1. Data
# 2. Aesthetics (aes) - what variables map to what visual properties
# 3. Geometries (geom) - type of plot
# 4. Themes - overall appearance

# ----------------------------------------------------------------------------
# 8.2 Boxplots
# ----------------------------------------------------------------------------
biomedical_data

Aleyna_nida <- ggplot(biomedical_data, aes(x = Status, y = Time)) +
  geom_boxplot() +
  labs(title = "Aleyna and Nida")
print(Aleyna_nida)

ggplot(biomedical_data, aes(x = Status, y = Time)) 

# Basic boxplot
p1 <- ggplot(biomedical_data, aes(x = Status, y = Time)) +
  geom_boxplot() +
  labs(title = "Basic Boxplot")
print(p1)

# Enhanced boxplot with color
p2 <- ggplot(biomedical_data, aes(x = Status, y = Time, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Survival Time by Patient Status",
    x = "Patient Status",
    y = "Survival Time (Days)",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("Alive" = "#00BA38", "Dead" = "#F8766D"))
print(p2)


p2 <- ggplot(biomedical_data, aes(x = Status, y = Time, fill = Status)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Survival Time by Patient Status",
    x = "Patient Status",
    y = "Mahmut",
    fill = "Status"
  ) +
  scale_fill_manual(values = c("Alive" = "#FFCCCB", "Dead" = "#ADD8E6"))
print(p2)


# Boxplot by multiple groups
p3 <- ggplot(biomedical_data, aes(x = Stage, y = Time, fill = Sex)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = "Survival Time by Stage and Sex",
    x = "Cancer Stage",
    y = "Survival Time (Days)"
  )
print(p3)

# ----------------------------------------------------------------------------
# 8.3 Bar Charts
# ----------------------------------------------------------------------------

# Simple bar chart
p4 <- ggplot(biomedical_data, aes(x = Stage)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Distribution of Cancer Stages",
    x = "Stage",
    y = "Count"
  )
print(p4)

# Horizontal bar chart
p5 <- ggplot(biomedical_data, aes(x = Stage, fill = Stage)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Distribution of Cancer Stages (Horizontal)",
    x = "Stage",
    y = "Count"
  )
print(p5)

# Stacked bar chart
p6 <- ggplot(biomedical_data, aes(x = Stage, fill = Status)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  labs(
    title = "Patient Status by Stage (Stacked)",
    x = "Stage",
    y = "Count"
  )
print(p6)

# Side-by-side bar chart
p7 <- ggplot(biomedical_data, aes(x = Stage, fill = Status)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(
    title = "Patient Status by Stage (Side-by-side)",
    x = "Stage",
    y = "Count"
  )
print(p7)

# ----------------------------------------------------------------------------
# 8.4 Histograms
# ----------------------------------------------------------------------------

# Basic histogram
p8 <- ggplot(biomedical_data, aes(x = Age)) +
  geom_histogram(bins = 10, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Age Distribution",
    x = "Age (years)",
    y = "Frequency"
  )
print(p8)

# Histogram with density overlay
p9 <- ggplot(biomedical_data, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "skyblue", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  theme_minimal() +
  labs(
    title = "Age Distribution with Density Curve",
    x = "Age (years)",
    y = "Density"
  )
print(p9)

# ----------------------------------------------------------------------------
# 8.5 Scatter Plots
# ----------------------------------------------------------------------------

# Basic scatter plot
p10 <- ggplot(biomedical_data, aes(x = Age, y = Time)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Age vs Survival Time",
    x = "Age (years)",
    y = "Survival Time (Days)"
  )
print(p10)

# Scatter plot with color by group
p11 <- ggplot(biomedical_data, aes(x = Age, y = Time, color = Status)) +
  geom_point(size = 3, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Age vs Survival Time by Status",
    x = "Age (years)",
    y = "Survival Time (Days)",
    color = "Patient Status"
  )
print(p11)

# Scatter plot with regression line
p12 <- ggplot(biomedical_data, aes(x = Age, y = Time)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue") +
  theme_minimal() +
  labs(
    title = "Age vs Survival Time with Regression Line",
    x = "Age (years)",
    y = "Survival Time (Days)"
  )
print(p12)

# ----------------------------------------------------------------------------
# 8.6 Violin Plots
# ----------------------------------------------------------------------------

p13 <- ggplot(biomedical_data, aes(x = Stage, y = Time, fill = Stage)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Survival Time Distribution by Stage",
    x = "Cancer Stage",
    y = "Survival Time (Days)"
  )
print(p13)

# ----------------------------------------------------------------------------
# 8.7 Faceting (Multiple Panels)
# ----------------------------------------------------------------------------

# Facet by one variable
p14 <- ggplot(biomedical_data, aes(x = Age, y = Time)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ Stage) +
  theme_minimal() +
  labs(
    title = "Age vs Survival Time by Stage",
    x = "Age (years)",
    y = "Survival Time (Days)"
  )
print(p14)

# Facet by two variables
p15 <- ggplot(biomedical_data, aes(x = Age, y = Time, color = Status)) +
  geom_point() +
  facet_grid(Sex ~ Stage) +
  theme_minimal() +
  labs(
    title = "Age vs Survival Time by Sex and Stage",
    x = "Age (years)",
    y = "Survival Time (Days)"
  )
print(p15)

# ============================================================================
# SECTION 9: SURVIVAL ANALYSIS
# ============================================================================

print("\n=== SURVIVAL ANALYSIS ===")

library(survival)
library(survminer)

# ----------------------------------------------------------------------------
# 9.1 Creating Survival Objects
# ----------------------------------------------------------------------------

# Create survival object
# Time: time to event or censoring
# Status: 1 = event occurred, 0 = censored
surv_object <- Surv(
  time = biomedical_data$Time,
  event = biomedical_data$Status == "Dead"
)

print("Survival object:")
print(head(surv_object))

# ----------------------------------------------------------------------------
# 9.2 Kaplan-Meier Survival Curves
# ----------------------------------------------------------------------------

# Overall survival (no grouping)
fit_overall <- survfit(surv_object ~ 1, data = biomedical_data)

print("\nOverall survival summary:")
print(summary(fit_overall))

# Plot overall survival
ggsurvplot(
  fit_overall,
  data = biomedical_data,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Kaplan-Meier Survival Curve",
  xlab = "Time (Days)",
  ylab = "Survival Probability"
)

# Survival by Status (this is just for demonstration)
fit_status <- survfit(surv_object ~ Status, data = biomedical_data)

ggsurvplot(
  fit_status,
  data = biomedical_data,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival by Initial Status",
  xlab = "Time (Days)",
  ylab = "Survival Probability"
)

# Survival by Sex
fit_sex <- survfit(surv_object ~ Sex, data = biomedical_data)

print("\nSurvival by Sex:")
print(summary(fit_sex))

ggsurvplot(
  fit_sex,
  data = biomedical_data,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  title = "Kaplan-Meier Curve by Sex",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  palette = c("#E69F00", "#56B4E9"),
  legend.labs = c("Female", "Male")
)

# Survival by Stage
fit_stage <- survfit(surv_object ~ Stage, data = biomedical_data)

ggsurvplot(
  fit_stage,
  data = biomedical_data,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  title = "Kaplan-Meier Curve by Stage",
  xlab = "Time (Days)",
  ylab = "Survival Probability"
)

# ----------------------------------------------------------------------------
# 9.3 Log-Rank Test
# ----------------------------------------------------------------------------

# Test difference between groups
surv_diff_sex <- survdiff(surv_object ~ Sex, data = biomedical_data)
print("\nLog-rank test (Sex):")
print(surv_diff_sex)

surv_diff_stage <- survdiff(surv_object ~ Stage, data = biomedical_data)
print("\nLog-rank test (Stage):")
print(surv_diff_stage)

# ----------------------------------------------------------------------------
# 9.4 Cox Proportional Hazards Model
# ----------------------------------------------------------------------------

# Univariate Cox model
cox_age <- coxph(surv_object ~ Age, data = biomedical_data)
print("\nCox model (Age):")
print(summary(cox_age))

cox_sex <- coxph(surv_object ~ Sex, data = biomedical_data)
print("\nCox model (Sex):")
print(summary(cox_sex))

# Multivariable Cox model
cox_multi <- coxph(surv_object ~ Age + Sex + Stage, data = biomedical_data)
print("\nMultivariable Cox model:")
print(summary(cox_multi))

# Forest plot of Cox model
ggforest(cox_multi, data = biomedical_data)

# ----------------------------------------------------------------------------
# 9.5 Testing Proportional Hazards Assumption
# ----------------------------------------------------------------------------

# Test proportional hazards assumption
cox_zph_test <- cox.zph(cox_multi)
print("\nProportional hazards test:")
print(cox_zph_test)

# Plot Schoenfeld residuals
plot(cox_zph_test)

# ============================================================================
# SECTION 10: USING BUILT-IN DATASETS FOR PRACTICE
# ============================================================================

print("\n=== WORKING WITH BUILT-IN DATASETS ===")

# ----------------------------------------------------------------------------
# 10.1 AML Dataset
# ----------------------------------------------------------------------------

# Load AML dataset from survival package
data("aml", package = "survival")

print("\nAML dataset structure:")
str(aml)

print("\nFirst few rows:")
print(head(aml))

# Check for missing values
print("\nMissing values:")
print(colSums(is.na(aml)))

# Create survival object
aml_surv <- Surv(aml$time, aml$status == 1)

# Fit survival model by treatment
aml_fit <- survfit(aml_surv ~ x, data = aml)

print("\nAML survival summary:")
print(summary(aml_fit))

# Plot Kaplan-Meier curve
ggsurvplot(
  aml_fit,
  data = aml,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  title = "AML Survival by Treatment",
  xlab = "Time (Weeks)",
  ylab = "Survival Probability",
  palette = c("#E69F00", "#56B4E9"),
  legend.labs = c("Maintained", "Not Maintained")
)

# Log-rank test
aml_diff <- survdiff(aml_surv ~ x, data = aml)
print("\nLog-rank test:")
print(aml_diff)

# Cox model
aml_cox <- coxph(aml_surv ~ x, data = aml)
print("\nCox model:")
print(summary(aml_cox))

# Test proportional hazards
aml_zph <- cox.zph(aml_cox)
print("\nProportional hazards test:")
print(aml_zph)

# ----------------------------------------------------------------------------
# 10.2 Lung Cancer Dataset
# ----------------------------------------------------------------------------

# Load lung dataset
data("lung", package = "survival")

print("\nLung cancer dataset structure:")
str(lung)

print("\nFirst few rows:")
print(head(lung))

# Summary
print("\nSummary:")
print(summary(lung))

# Create survival object
lung_surv <- Surv(lung$time, lung$status == 2)  # status 2 = dead

# Survival by sex
lung_fit_sex <- survfit(lung_surv ~ sex, data = lung)

ggsurvplot(
  lung_fit_sex,
  data = lung,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  title = "Lung Cancer Survival by Sex",
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  legend.labs = c("Male", "Female")
)

# Multivariable Cox model
lung_cox <- coxph(
  lung_surv ~ age + sex + ph.ecog + ph.karno + pat.karno,
  data = lung
)
print("\nLung cancer Cox model:")
print(summary(lung_cox))

# ============================================================================
# SECTION 11: GENE EXPRESSION ANALYSIS WITH DESeq2
# ============================================================================

print("\n=== GENE EXPRESSION ANALYSIS ===")

# Note: DESeq2 requires installation via Bioconductor
# BiocManager::install("DESeq2")

library(DESeq2)

# ----------------------------------------------------------------------------
# 11.1 Simulated RNA-seq Data
# ----------------------------------------------------------------------------

set.seed(123)

# Create simulated count matrix (100 genes x 10 samples)
count_data <- matrix(
  rnbinom(1000, mu = 10, size = 1),
  nrow = 100,
  ncol = 10
)

# Add row and column names
rownames(count_data) <- paste0("Gene", 1:100)
colnames(count_data) <- paste0("Sample", 1:10)

print("Count matrix (first 5 genes, first 5 samples):")
print(count_data[1:5, 1:5])

# Create sample information
col_data <- data.frame(
  condition = factor(rep(c("Control", "Treated"), each = 5)),
  row.names = colnames(count_data)
)

print("\nSample information:")
print(col_data)

# ----------------------------------------------------------------------------
# 11.2 DESeq2 Analysis
# ----------------------------------------------------------------------------

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ condition
)

print("\nDESeq2 dataset:")
print(dds)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

print("\nDESeq2 results summary:")
print(summary(res))

# View results
print("\nTop differentially expressed genes:")
print(head(res[order(res$padj), ]))

# Filter significant genes
sig_genes <- res[which(res$padj < 0.05), ]
print(paste("\nNumber of significant genes (padj < 0.05):", nrow(sig_genes)))

# ----------------------------------------------------------------------------
# 11.3 Visualization of Results
# ----------------------------------------------------------------------------

# MA plot
plotMA(res, main = "MA Plot", ylim = c(-2, 2))

# Volcano plot (manual)
plot_data <- as.data.frame(res)
plot_data$significant <- plot_data$padj < 0.05

ggplot(plot_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

# PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

# ----------------------------------------------------------------------------
# 11.4 Heatmap of Top Genes
# ----------------------------------------------------------------------------

# Get top 20 genes
top_genes <- head(order(res$padj), 20)
top_gene_names <- rownames(res)[top_genes]

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
top_counts <- normalized_counts[top_gene_names, ]

# Simple heatmap
heatmap(
  as.matrix(top_counts),
  scale = "row",
  main = "Top 20 Differentially Expressed Genes"
)

# ============================================================================
# SECTION 12: PRINCIPAL COMPONENT ANALYSIS (PCA)
# ============================================================================

print("\n=== PRINCIPAL COMPONENT ANALYSIS ===")

# ----------------------------------------------------------------------------
# 12.1 PCA on Gene Expression Data
# ----------------------------------------------------------------------------

# Create gene expression matrix
set.seed(456)
gene_expr_matrix <- matrix(
  rnorm(1000, mean = 5, sd = 2),
  nrow = 100,
  ncol = 10
)
rownames(gene_expr_matrix) <- paste0("Gene", 1:100)
colnames(gene_expr_matrix) <- paste0("Sample", 1:10)

# Perform PCA
pca_result <- prcomp(t(gene_expr_matrix), center = TRUE, scale. = TRUE)

print("\nPCA summary:")
print(summary(pca_result))

# Variance explained
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
print("\nVariance explained by each PC:")
print(head(variance_explained))

# Scree plot
plot(
  variance_explained,
  type = "b",
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Scree Plot"
)

# Biplot
biplot(
  pca_result,
  main = "PCA Biplot",
  cex = 0.7
)

# Plot PC1 vs PC2
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
pca_data$Group <- rep(c("Control", "Treated"), each = 5)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA Plot",
    x = paste0("PC1 (", round(variance_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2] * 100, 1), "%)")
  )

# ============================================================================
# SECTION 13: ADVANCED TOPICS OVERVIEW
# ============================================================================

print("\n=== ADVANCED TOPICS OVERVIEW ===")

# ----------------------------------------------------------------------------
# 13.1 Single-Cell RNA-seq with Seurat (Brief Example)
# ----------------------------------------------------------------------------

print("\n--- Single-Cell RNA-seq ---")

# Note: Seurat requires Bioconductor installation
# BiocManager::install("Seurat")

library(Seurat)

# Create small example Seurat object
set.seed(789)
counts <- matrix(rpois(200, lambda = 5), ncol = 10)
rownames(counts) <- paste0("Gene", 1:20)
colnames(counts) <- paste0("Cell", 1:10)

seurat_obj <- CreateSeuratObject(counts = counts, project = "Example")

print("\nSeurat object:")
print(seurat_obj)

# Basic QC violin plot
VlnPlot(seurat_obj, features = "nFeature_RNA")

# ----------------------------------------------------------------------------
# 13.2 TCGA Data Access (Conceptual)
# ----------------------------------------------------------------------------

print("\n--- TCGA Data Access ---")

# Note: This is a conceptual example
# Actual TCGA data download requires significant time and storage

# library(TCGAbiolinks)
# library(SummarizedExperiment)
#
# # Query TCGA data
# query <- GDCquery(
#   project = "TCGA-BRCA",
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )
#
# # Download data
# GDCdownload(query)
#
# # Prepare data
# data <- GDCprepare(query)
#
# # Access expression data
# expr_data <- assay(data)

print("TCGA data access requires TCGAbiolinks package")
print("See commented code above for example usage")

# ----------------------------------------------------------------------------
# 13.3 ChIP-seq Analysis (Conceptual)
# ----------------------------------------------------------------------------

print("\n--- ChIP-seq Analysis ---")

# Note: Requires peak files and genome annotations
# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#
# # Read peak file
# peaks <- readPeakFile("chipseq_peaks.bed")
#
# # Annotate peaks
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# peakAnno <- annotatePeak(peaks, TxDb = txdb)
#
# # Plot annotation
# plotAnnoBar(peakAnno)

print("ChIP-seq analysis requires peak files and genome annotations")
print("See commented code above for example usage")

# ============================================================================
# SECTION 14: SAVING AND LOADING ANALYSIS RESULTS
# ============================================================================

print("\n=== SAVING ANALYSIS RESULTS ===")

# Create output directory
if (!dir.exists("output")) {
  dir.create("output")
}

# Save individual objects as RDS
saveRDS(biomedical_data, file = "output/biomedical_data.rds")
saveRDS(surv_object, file = "output/survival_object.rds")
saveRDS(cox_multi, file = "output/cox_model.rds")

print("Saved analysis objects as RDS files")

# Save DESeq2 results
if (exists("res")) {
  saveRDS(res, file = "output/deseq2_results.rds")
  write.csv(as.data.frame(res), "output/deseq2_results.csv", row.names = TRUE)
  print("Saved DESeq2 results")
}

# Save workspace
save.image(file = "output/analysis_workspace.RData")
print("Saved entire workspace")

# ============================================================================
# SECTION 15: BEST PRACTICES AND TIPS
# ============================================================================

print("\n=== BEST PRACTICES ===")

cat("
BEST PRACTICES FOR R PROGRAMMING:

1. PROJECT ORGANIZATION:
   - Use RStudio Projects (.Rproj files)
   - Create folders: data/, output/, scripts/, figures/
   - Use the 'here' package for relative paths

2. CODE ORGANIZATION:
   - Comment your code extensively
   - Use meaningful variable names
   - Break long scripts into sections
   - Load all libraries at the beginning

3. DATA MANAGEMENT:
   - Always check for missing values
   - Verify data types after import
   - Create backups of original data
   - Document data transformations

4. VERSION CONTROL:
   - Use Git for version control
   - Commit frequently with clear messages
   - Keep sensitive data out of repositories

5. REPRODUCIBILITY:
   - Set random seeds for reproducible results
   - Document package versions (sessionInfo())
   - Use RMarkdown for reports
   - Share code and data when possible

6. PERFORMANCE:
   - Vectorize operations when possible
   - Use apply() family instead of loops for large datasets
   - Profile code to identify bottlenecks
   - Consider data.table for very large datasets

7. ERROR HANDLING:
   - Use try() and tryCatch() for robust code
   - Validate inputs and outputs
   - Write informative error messages

8. DOCUMENTATION:
   - Comment complex operations
   - Create README files for projects
   - Document custom functions
   - Keep analysis notebooks
")

# ============================================================================
# SECTION 16: USEFUL FUNCTIONS REFERENCE
# ============================================================================

print("\n=== USEFUL FUNCTIONS REFERENCE ===")

cat("
COMMONLY USED R FUNCTIONS:

DATA INSPECTION:
- head(), tail()         : View first/last rows
- str()                  : Structure of object
- summary()              : Summary statistics
- dim(), nrow(), ncol()  : Dimensions
- names(), colnames()    : Column names
- class(), typeof()      : Object type

DATA MANIPULATION:
- subset(), filter()     : Select rows
- select()               : Select columns
- mutate()               : Create/modify columns
- arrange()              : Sort data
- group_by()             : Group data
- summarize()            : Summary statistics

MISSING DATA:
- is.na()                : Check for NA
- complete.cases()       : Complete rows
- na.omit()              : Remove NA rows
- na.rm = TRUE           : Ignore NA in calculations

STATISTICAL FUNCTIONS:
- mean(), median()       : Central tendency
- sd(), var()            : Dispersion
- min(), max(), range()  : Range
- quantile()             : Quantiles
- cor()                  : Correlation
- t.test(), wilcox.test(): Statistical tests

SURVIVAL ANALYSIS:
- Surv()                 : Create survival object
- survfit()              : Kaplan-Meier
- coxph()                : Cox model
- survdiff()             : Log-rank test

PLOTTING:
- ggplot()               : Initialize plot
- geom_point()           : Scatter plot
- geom_line()            : Line plot
- geom_boxplot()         : Box plot
- geom_histogram()       : Histogram
- facet_wrap()           : Multiple panels
")

# ============================================================================
# SECTION 17: SESSION INFORMATION
# ============================================================================

print("\n=== SESSION INFORMATION ===")

# Print session information for reproducibility
sessionInfo()

# ============================================================================
# END OF LECTURE
# ============================================================================

cat("\n
================================================================================
CONGRATULATIONS! You've completed the BENG 628 R Programming Lecture!
================================================================================

NEXT STEPS:
1. Practice with your own biomedical datasets
2. Explore additional packages relevant to your research
3. Join R communities (Stack Overflow, RStudio Community, Bioconductor)
4. Read package documentation and vignettes
5. Take online courses (Coursera, DataCamp, etc.)

USEFUL RESOURCES:
- R for Data Science: https://r4ds.had.co.nz/
- Bioconductor: https://www.bioconductor.org/
- RStudio Cheatsheets: https://www.rstudio.com/resources/cheatsheets/
- CRAN Task Views: https://cran.r-project.org/web/views/

HAPPY CODING!
================================================================================
")