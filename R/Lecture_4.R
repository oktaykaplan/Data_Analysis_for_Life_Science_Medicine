# BENG 628 - Lecture 4: Exploratory Data Analysis for Life Sciences
#We will expand on foundational R skills by covering essential data visualization and basic statistical analysis, which are crucial for data science. Here are some topics that would complement what I've taught so far:

#1. Data Visualization with ggplot2
#Introduce ggplot2, a powerful visualization package in R.
#Basics of ggplot(), aes(), and geom_ functions.
#Common plots:
#Scatter plots: geom_point()
#Bar plots: geom_bar()
#Line charts: geom_line()
#Histograms: geom_histogram()
#Boxplots: geom_boxplot()
#Customizing plots with labels, themes, and colors.
#2. Handling Missing Data
#Identify missing values: is.na()
#Remove missing data: na.omit(), drop_na()
#Replace missing values: mutate() with ifelse(), replace_na()
#Summarize missing values: summary(), colSums(is.na(data))
#3. Basic Data Transformation with tidyr
#Reshaping data: pivot_longer(), pivot_wider() (We covered in the lecture 2 and remember it again)
#Separating and uniting columns: separate(), unite()
#Handling categorical variables: factor(), recode()
#4. Basic Statistical Summaries and Functions
#Measures of central tendency: mean(), median(), mode()
#Measures of spread: sd(), var(), quantile(), IQR()
#Simple hypothesis testing: t.test(), cor()
#Creating frequency tables: table(), prop.table()
#5. Writing Functions in R
#Creating user-defined functions using function()
#Using apply(), sapply(), lapply() for iteration
#Vectorized operations vs. loops

# Advanced Data Handling & Visualization in R for Life Sciences

#Objective: Visualize gene expression data from different tissue types
#Let's generate gene expression 

# Load necessary library
library(ggplot2)


# Simulated gene expression data
gene_expression <- data.frame(
  Gene = rep(c("BRCA1", "TP53", "MYC"), each = 5),
  Tissue = rep(c("Liver", "Lung", "Brain", "Kidney", "Heart"), times = 3),
  Expression = c(12, 8, 15, 7, 10, 5, 3, 6, 8, 4, 14, 10, 12, 11, 9)
)
view(gene_expression)
# Simulated gene expression data 2
gene_kinesin_expression <- data.frame(
  Gene = rep(c("KIF1A", "KIF2A", "KIF3A", "KIF5B", "KIF11", "KIF14", "KIF21A", "KIF22", "BRCA1", "TP53", "MYC"), each = 5),
  Tissue = rep(c("Liver", "Lung", "Brain", "Kidney", "Heart"), times = 11), 
  Expression = sample(5:20, 55, replace = TRUE) 
)

view(gene_kinesin_expression)

# Create a boxplot for gene expression
ggplot(gene_expression, aes(x = Tissue, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Across Tissues",
       x = "Tissue Type",
       y = "Expression Level")


# Create a boxplot for gene expression
ggplot(gene_kinesin_expression, aes(x = Tissue, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Across Tissues",
       x = "Tissue Type",
       y = "Expression Level")


#2. Handling Missing Data (DNA Methylation Example)
#Objective: Identify and handle missing CpG methylation values
#Dataset: Simulated CpG site methylation levels (some missing values).
# Simulated DNA methylation dataset
methylation <- data.frame(
  CpG_Site = paste0("CpG_", 1:10),
  Methylation_Level = c(85, 78, NA, 67, 90, 72, 88, NA, 80, 75)
)

# Check missing values
print(methylation)
print(sum(is.na(methylation$Methylation_Level)))  # Count missing values


# Replace NA values with mean methylation level
methylation$Methylation_Level <- ifelse(is.na(methylation$Methylation_Level),
                                        mean(methylation$Methylation_Level, na.rm = TRUE),
                                        methylation$Methylation_Level)
# Print cleaned dataset
print(methylation)
#Exercise: Instead of using the mean, replace missing values with the median and observe the difference.

#3. Data Reshaping with tidyr (GWAS Example)
# Objective: Transform GWAS (genome-wide association study) results from wide to long format.
#Dataset: Simulated SNP association p-values across two conditions.

# Load library
library(tidyr)

# Simulated GWAS data
gwas_results <- data.frame(
  SNP = c("rs123", "rs456", "rs789"),
  Condition_A = c(0.02, 0.05, 0.001),
  Condition_B = c(0.03, 0.08, 0.0005)
)

# Convert from wide to long format
gwas_long <- pivot_longer(gwas_results, cols = starts_with("Condition"), 
                          names_to = "Condition", values_to = "p_value")

print(gwas_long)

#Exercise: Convert the long format back into wide format using pivot_wider().


#4. Basic Statistical Summaries (Genetic Variant Example)
#Objective: Analyze the allele frequency distribution of a genetic variant.
#Dataset: Simulated allele frequency data for a genetic variant in a population.
# Simulated allele frequencies
allele_freq <- c(0.12, 0.34, 0.45, 0.29, 0.21, 0.38, 0.25, 0.32, 0.41, 0.37)

# Compute statistics
mean_allele <- mean(allele_freq)
sd_allele <- sd(allele_freq)
median_allele <- median(allele_freq)
quantiles_allele <- quantile(allele_freq)

# Print results
print(paste("Mean Allele Frequency:", mean_allele))
print(paste("Standard Deviation:", sd_allele))
print(paste("Median:", median_allele))
print("Quantiles:")
print(quantiles_allele)

#Exercise: Create a histogram of the allele frequency distribution.

#5. Writing Functions in R (Hardy-Weinberg Equilibrium Calculator)
#Objective: Create a function to calculate Hardy-Weinberg Equilibrium (HWE) proportions.
# Function to calculate Hardy-Weinberg proportions
hwe_calc <- function(p) {
  q <- 1 - p  # Recessive allele frequency
  AA <- p^2    # Homozygous dominant
  Aa <- 2 * p * q  # Heterozygous
  aa <- q^2    # Homozygous recessive
  
  return(data.frame(Genotype = c("AA", "Aa", "aa"),
                    Frequency = c(AA, Aa, aa)))
}

# Test function
hwe_calc(0.6)

#Let's use real-world datasets from Ensembl, GEO, or other genomic databases
#Gene Expression Analysis (GEO Data)
# Objective: Load and visualize gene expression data from GEO (Gene Expression Omnibus).
#Dataset: Microarray or RNA-seq data from GEO.
# Tools: GEOquery (for downloading data), ggplot2 (for visualization).
install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

#Download and Load Gene Expression Data
#Example using GSE37745 (a dataset related to gene expression in kidney diseases).
options(timeout = 300)  # Set timeout to 300 seconds globally
gse <- getGEO("GSE37745", GSEMatrix = TRUE) 
data <- exprs(gse[[1]])  # Extract expression values

# Convert to a data frame
gene_data <- data.frame(Gene = rownames(data), Expression = data[, 1])  # Take the first sample
head(gene_data)

#Visualize Gene Expression
library(ggplot2)
library(dplyr)
library(tibble)

set.seed(123)  # For reproducibility
data <- as.data.frame(matrix(runif(200, min = 5, max = 20), nrow = 20))

# Convert to a proper data frame
data <- as.data.frame(as.matrix(data))

# Add gene names
rownames(data) <- paste0("Gene", 1:nrow(data))

# Convert row names to column and compute mean expression
Mean_Expression_gene <- data %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Mean_Expression = rowMeans(across(where(is.numeric)))) 

# Ensure Mean_Expression is numeric
Mean_Expression_gene$Mean_Expression <- as.numeric(Mean_Expression_gene$Mean_Expression)

# Select top 20 highly expressed genes
top_genes <- Mean_Expression_gene %>%
  arrange(desc(Mean_Expression)) %>%
  dplyr::slice(1:20)

#Check if Mean_Expression exists in top_genes
colnames(top_genes)
str(top_genes)

# Plot the mean expression of top 20 genes
ggplot(top_genes, aes(x = reorder(Gene, Expression), y = Expression)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Highly Expressed Genes (Mean Across Samples)",
       x = "Gene", y = "Expression Level")

#Genetic Variant Analysis (Ensembl Data)
#Objective: Retrieve genetic variant annotations from Ensembl.
#Dataset: Human SNP annotation from Ensembl BioMart.
#Tools: biomaRt (for querying Ensembl), dplyr (for filtering variants).
BiocManager::install("biomaRt", force = TRUE)

library(biomaRt)
#Connect to Ensembl and Query Variants
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# Retrieve variant information for a specific gene (e.g., TP53)
variants <- getBM(attributes = c("refsnp_id", "chr_name", "start_position", 
                                 "allele", "consequence_type"),
                  filters = "ensembl_gene", 
                  values = "ENSG00000141510", 
                  mart = ensembl)

head(variants)
#Filter and Analyze Variants
library(dplyr)

# Filter for only missense mutations
missense_variants <- variants %>% 
  filter(consequence_type == "missense_variant")

print(missense_variants)
#Modify the script to retrieve variants for a different gene (e.g., BRCA1 or CFTR)

#GWAS Summary Statistics (UK Biobank)
# Objective: Analyze genome-wide association study (GWAS) data from UK Biobank.
# Dataset: UK Biobank GWAS summary statistics (Public Data).
# Tools: data.table (for handling large datasets), ggplot2 (for visualization).
#Load a GWAS Summary File
#Example: UK Biobank GWAS results for height (hypothetical).
library(data.table)

# Load the GWAS summary statistics file (example)
gwas <- fread("https://www.ebi.ac.uk/gwas/api/search/downloads/alternative")

# Check structure
head(gwas)
colnames(gwas)


#Manhattan Plot (GWAS Visualization)
library(ggplot2)

# Rename columns for easier reference
gwas <- gwas %>%
  dplyr::rename(chr = CHR_ID, p_value = `P-VALUE`)

# Convert chromosome column to numeric (handle "X" and "Y" if present)
gwas$chr <- as.character(gwas$chr)
gwas$chr[gwas$chr == "X"] <- "23"
gwas$chr[gwas$chr == "Y"] <- "24"
gwas$chr[gwas$chr == "MT"] <- "25"
gwas$chr <- as.numeric(gwas$chr)

# Convert p_value to numeric
gwas$p_value <- as.numeric(gwas$p_value)

# Remove NAs
gwas <- na.omit(gwas)

# Filter significant SNPs (p-value < 5e-8)
significant_snps <- gwas[gwas$p_value < 5e-8, ]

# Plot Manhattan plot
ggplot(significant_snps, aes(x = chr, y = -log10(p_value))) +
  geom_point(aes(color = as.factor(chr))) +
  theme_minimal() +
  labs(title = "Manhattan Plot of GWAS Results",
       x = "Chromosome",
       y = "-log10(P-value)")


#Exercise: Find the top 10 most significant SNPs and plot their effect sizes.

#Single-Cell RNA-Seq Data (Human Brain Atlas)
#Objective: Explore single-cell RNA-seq expression of brain cell types.
# Dataset: Human Cell Atlas (single-cell RNA-seq).
#Tools: Seurat (for single-cell analysis).

#Install and Load Seurat
BiocManager::install("Seurat")
library(Seurat)

#Load Single-Cell RNA-Seq Data
#Example dataset: PBMC (Peripheral Blood Mononuclear Cells)
pbmc <- Read10X_h5(filename = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5")

pbmc <- CreateSeuratObject(counts = pbmc, project = "PBMC3K")
pbmc

#Exercise: Cluster the data using FindClusters() and visualize with DimPlot().


#Create a Methylation Dataset

# Generate random methylation levels (higher methylation may indicate lower expression)
set.seed(123)  # For reproducibility
gene_kinesin_expression$Methylation <- sample(20:80, 55, replace = TRUE)

# View the updated dataset
View(gene_kinesin_expression)

#Correlate Expression and Methylation
#Let's calculate the correlation between gene expression and methylation levels.
correlation_result <- cor.test(gene_kinesin_expression$Expression, gene_kinesin_expression$Methylation)

# Print correlation coefficient and p-value
print(correlation_result)

#Visualize Expression vs. Methylation
library(ggplot2)

ggplot(gene_kinesin_expression, aes(x = Methylation, y = Expression, color = Gene)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  theme_minimal() +
  labs(title = "Expression vs. Methylation Pattern",
       x = "Methylation (%)",
       y = "Gene Expression (arbitrary units)")

#Cluster Genes Based on Methylation-Expression Patterns
#Cluster genes using hierarchical clustering (heatmap) and k-means clustering.

library(pheatmap)
#Hierarchical Clustering with Heatmap
# Prepare data matrix
#Groups genes based on similarity in methylation and expression levels.
clustering_data <- gene_kinesin_expression %>%
  dplyr::select(Gene, Expression, Methylation) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(Expression = mean(Expression), Methylation = mean(Methylation)) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# Create heatmap
pheatmap(clustering_data, scale = "row", clustering_distance_rows = "euclidean",
         clustering_method = "complete", main = "Gene Clustering by Methylation & Expression")

#Compare Methylation Between Tissues (Brain vs. Liver)
#Compare methylation levels between Brain and Liver using a boxplot.
#Helps compare how methylation patterns differ in these two tissues.
ggplot(gene_kinesin_expression %>% filter(Tissue %in% c("Brain", "Liver")), 
       aes(x = Tissue, y = Methylation, fill = Tissue)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Methylation Levels in Brain vs. Liver",
       y = "Methylation (%)",
       x = "Tissue")
#Retrieve real-world methylation data from GEO (Gene Expression Omnibus) using the GEOquery package.
#Download Methylation Data from GEO 
#Loads real methylation levels from GEO (e.g., for cancer samples).
#Perform differential methylation analysis (e.g., between cancer and normal samples)?
# Integrate real expression data (e.g., RNA-seq from TCGA)?
library(GEOquery)
# Example dataset (adjust for your gene of interest)
geo_data <- getGEO("GSE42865", GSEMatrix = TRUE)
methylation_data <- exprs(geo_data[[1]])

# View methylation data
head(methylation_data)


