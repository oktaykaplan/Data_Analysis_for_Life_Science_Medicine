# BENG 628 - Lecture 3: Exploratory Data Analysis for Life Sciences

#Lecture Goals
#II. Advanced Data Wrangling
#Multivariate Data Exploration
# Advanced Outlier Detection
# Feature Engineering & Transformations
#Time-Series & Longitudinal Data Handling
#Reproducible Research in R
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

#Data Wrangling and Preparation:
#pivot_longer() and pivot_wider():  These are crucial for reshaping data between wide and long formats.  
#Many statistical analyses and visualizations require data in a specific format
wide_data <- data.frame(patient = 1:3, week1 = c(10, 12, 11), week2 = c(15, 14, 16))
view(wide_data)

long_data <- wide_data %>% pivot_longer(cols = starts_with("week"), names_to = "week", values_to = "measurement")
view(long_data)

# Long to wide
long_data %>% pivot_wider(names_from = week, values_from = measurement)


#String Manipulation (stringr package)
#Working with text data is very common in biology and medicine.  
#Scientists often need to clean, extract, or manipulate strings

#Extracting parts of strings
gene_names <- c("TP53_mutant", "KRAS_wildtype", "EGFR_mutant")
view(gene_names)
#Let's extract using mutant or wildtype
gene_types <- str_extract(gene_names, "(mutant|wildtype)")
view(gene_types)

# Replacing patterns
dna_sequence <- "ATGCCTAG"
view(dna_sequence)
#Replace a specific someting with another
mutated_sequence <- str_replace(dna_sequence, "CCT", "GGT")
view(mutated_sequence)

#Splitting strings
patient_info <- "John Doe, 35, Male"
view(patient_info)
info_parts <- str_split(patient_info, ", ")
view(info_parts)

library(mice)
#Beyond just filtering, scientists need strategies for dealing with missing data
# Example: Simple imputation (replace NAs with the mean) - not always best!
data_with_na <- data.frame(values = c(1, 2, NA, 4, 5))
mean_values <- mean(data_with_na$values, na.rm = TRUE)
data_with_na$values[is.na(data_with_na$values)] <- mean_values

# More advanced imputation (using mice)
imputed_data <- mice(data_with_na, m = 5, method = "pmm", printFlag = FALSE) # Create 5 imputed datasets
final_data <- complete(imputed_data, 1) # Choose the first imputed dataset

# Create sample gene expression data
set.seed(123) # for reproducibility

# Generate sample data with 5 genes, 3 conditions, 2 replicates
gene_expr <- data.frame(
  ENSG00000001 = rnorm(6, mean = 100, sd = 10),
  ENSG00000002 = rnorm(6, mean = 500, sd = 50),
  ENSG00000003 = rnorm(6, mean = 250, sd = 25),
  ENSG00000004 = rnorm(6, mean = 1000, sd = 100),
  ENSG00000005 = rnorm(6, mean = 50, sd = 5),
  condition = rep(c("control", "treatment_A", "treatment_B"), each = 2),
  replicate = rep(1:2, 3)
)

# Add time_point as numeric variable
gene_expr <- gene_expr %>%
  mutate(
    time_point = case_when(
      condition == "control" ~ 0,
      condition == "treatment_A" ~ 24,
      condition == "treatment_B" ~ 48,
      TRUE ~ NA_real_
    ),
    condition = factor(condition, levels = c("control", "treatment_A", "treatment_B"))
  )

# Transform to long format (keeping all metadata)
long_data <- gene_expr %>%
  pivot_longer(
    cols = starts_with("ENSG"),
    names_to = "gene_id",
    values_to = "expression"
  )

# Grouped summary by condition only (time_point is redundant)
grouped_summary <- gene_expr %>%
  group_by(condition) %>%
  summarise(across(starts_with("ENSG"),
                   list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}"))

print(grouped_summary)

# Technical replicates summary
summarized_expr <- long_data %>%
  group_by(gene_id, condition) %>%
  summarise(
    mean_expr = mean(expression),
    cv = sd(expression)/mean_expr * 100,
    n_detected = sum(expression > 0),
    .groups = "drop"
  )

# Visualization 1: Boxplot of expression by condition
ggplot(long_data, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~gene_id, scales = "free_y") +
  labs(title = "Gene Expression Distribution by Condition",
       x = "Condition",
       y = "Expression Level") +
  theme_minimal()

# Visualization 2: Line plot of expression over time
time_summary <- long_data %>%
  group_by(gene_id, time_point) %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr = sd(expression),
    .groups = "drop"
  )

ggplot(time_summary, aes(x = time_point, y = mean_expr, color = gene_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), width = 3) +
  labs(title = "Gene Expression Trajectories Over Time",
       x = "Time Point (hours)",
       y = "Mean Expression Level",
       color = "Gene ID") +
  theme_minimal()

# Visualization 3: Heatmap of mean expression
ggplot(summarized_expr, aes(x = condition, y = gene_id, fill = mean_expr)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Mean Gene Expression Heatmap",
       x = "Condition",
       y = "Gene ID",
       fill = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(summarized_expr) # Print the summary results


#Hands-on Exercise 1 (30 min)
#Working with real RNA-seq count data:
# Load example RNA-seq dataset
data("airway", package = "airway")

# Convert to DESeq2 object
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# Normalize counts
normalized_counts <- counts(dds, normalized=TRUE)

# Exercise tasks:
# 1. Convert normalized counts to tidy format
# 2. Calculate mean expression per condition
# 3. Identify highly variable genes
#Part 2: Visualization Techniques and Best Practices
#I. Advanced ggplot2 (45 min)
#A. Complex Visualizations
# Enhanced volcano plot
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulated), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("gray", "red", "blue")) +
  theme_minimal() +
  labs(title = "Differential Expression Analysis",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

# Advanced heatmap
expr_matrix <- gene_expr %>%
  select(-c(condition, time_point)) %>%
  as.matrix()

pheatmap(expr_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = FALSE,
         annotation = annotation_df)
#B. Custom Themes and Color Palettes
# Create publication-ready theme
theme_publication <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      axis.text = element_text(color = "black"),
      legend.position = "bottom"
    )
}

# Color-blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Hands-on Exercise 2 (45 min)
#Create publication-quality figures:

#  Generate a multi-panel figure showing:

#  Expression distributions
#PCA plot
#Heatmap of top genes


#Customize themes and colors
#Add appropriate labels and legends

#Part 3: Statistical Analysis and Interpretation
#I. Advanced Statistical Methods (25 min)
# Multiple testing approaches

results_df %>%
  mutate(
    bonferroni = p.adjust(pvalue, method = "bonferroni"),
    bh = p.adjust(pvalue, method = "BH"),
    qvalue = qvalue(pvalue)$qvalues
  )

# Effect size calculations
gene_expr %>%
  group_by(gene_id) %>%
  summarise(
    cohens_d = cohens_d(expression ~ condition)$effsize,
    cliff_delta = cliff.delta(expression ~ condition)$estimate
  )
#II. Biological Interpretation (20 min)
# Option 1: If you have a data frame with gene IDs
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

# Example ENSEMBL IDs
ensembl_ids <- c("ENSG0000013579", "ENSG0000024680", "ENSG0000011223") 

# Load Ensembl database with human dataset
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
# Try a Different Mirror Manually
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "asia")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
# Use useMart() Instead of useEnsembl()
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")


# Check if dataset is properly set
print(listDatasets(mart))  # Ensure "hsapiens_gene_ensembl" is listed

# Convert ENSEMBL to Entrez IDs
entrez_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                    filters = "ensembl_gene_id", 
                    values = ensembl_ids, 
                    mart = mart)

# Extract the Entrez IDs
significant_genes <- na.omit(entrez_ids$entrezgene_id)  # Remove NAs

# Perform GO enrichment analysis 
ego <- enrichGO(gene = significant_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP") 
print(ego)
nalysis
library(STRING)
string_db <- STRINGdb$new()
string_interactions <- string_db$get_interactions(significant_genes)

# Load dataset
# Assume we have a dataset "life_science_data.csv" with biological measurements
data <- read_csv("life_science_data.csv")

# Preview dataset
head(data)
summary(data)
str(data)

glimpse(data)  # Alternative to str() for a better overview

# Handling missing values
missing_values <- colSums(is.na(data))  # Count missing values per column
missing_values

data_clean <- na.omit(data)  # Remove rows with missing values

# Basic statistics
summary(data_clean)

# Distribution of numerical variables
num_cols <- data_clean %>% select(where(is.numeric))

ggplot(stack(num_cols), aes(values)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.6) +
  facet_wrap(~ind, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Numerical Variables")

# Boxplots to check for outliers
ggplot(stack(num_cols), aes(x = ind, y = values)) +
  geom_boxplot(fill = "orange", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Boxplots of Numerical Variables", x = "Variable", y = "Value") +
  coord_flip()

# Correlation matrix
cor_matrix <- cor(num_cols, use = "complete.obs")
cor_matrix

# Visualizing correlation matrix
library(ggcorrplot)
ggcorrplot(cor_matrix, lab = TRUE, outline.color = "white")

# Relationship between two selected variables
# Assume "GeneExpression" and "ProteinLevels" are relevant variables
ggplot(data_clean, aes(x = GeneExpression, y = ProteinLevels)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", col = "red") +
  theme_minimal() +
  labs(title = "Scatter Plot of Gene Expression vs. Protein Levels",
       x = "Gene Expression", y = "Protein Levels")

# Principal Component Analysis (PCA)
scaled_data <- scale(num_cols)
pca_res <- prcomp(scaled_data, center = TRUE, scale. = TRUE)

# Scree plot
ggplot(data.frame(PC = seq_along(pca_res$sdev), Variance = pca_res$sdev^2),
       aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Scree Plot for PCA")

# PCA Biplot
autoplot(pca_res, data = data_clean, colour = 'CategoryVariable') +
  theme_minimal() +
  labs(title = "PCA Biplot")

# Clustering analysis
library(cluster)
kmeans_res <- kmeans(num_cols, centers = 3)  # Assume 3 clusters

data_clean$cluster <- as.factor(kmeans_res$cluster)

ggplot(data_clean, aes(x = GeneExpression, y = ProteinLevels, color = cluster)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "K-Means Clustering of Gene Expression and Protein Levels")

# Summary of clusters
aggregate(num_cols, by = list(cluster = kmeans_res$cluster), mean)


#Multivariate Data Exploration (PCA, Clustering)
#When datasets have multiple features, we need ways to visualize and simplify the data.
# PCA - Principal Component Analysis

library(FactoMineR)
library(factoextra)

# Simulated gene expression data (7 patients, 5 genes)
#The PCA (Principal Component Analysis) biplot is a powerful way to visualize the relationships between variables (in this case, genes) and observations (in this case, patients). 
#Let’s break down what this plot shows and how to interpret it:
#Axes (PC1 and PC2):
#The axes in the plot represent the principal components (PCs). These are the directions in the data that capture the most variance.
#PC1 (horizontal axis): This axis captures the most variance in the data. It represents the largest spread of the data.
#PC2 (vertical axis): This axis captures the second-largest variance, orthogonal (at a 90-degree angle) to PC1.

gene_data <- data.frame(
  Patient = paste0("P", 1:7),
  ARL13B = c(10, 15, 20, 18, 25, 30, 28),
  IFT88 = c(5, 8, 6, 12, 10, 15, 18),
  IFT81 = c(50, 55, 60, 58, 65, 70, 75),
  BBS1 = c(100, 120, 110, 105, 130, 140, 135)
)

# Perform PCA (excluding Patient column)
pca_result <- PCA(gene_data[, -1], scale.unit = TRUE, graph = FALSE)

#Grouping Variable as Factor: The habillage argument is set to factor(gene_data$Patient) 
#to ensure that the Patient column is treated as a categorical variable for color grouping.
# Create PCA biplot
# If patients are clustered together in the plot, it indicates that they have similar gene expression profiles across the selected genes.
#If ARL13B has a long arrow pointing rightward along PC1, and IFT88 has a short arrow, it means GeneA is more important in distinguishing patients along PC1.
#Dim1 (89.5%): The horizontal axis represents the first principal component (PC1).
#Dim2 (6.8%): The vertical axis represents the second principal component (PC2). It explains 6.8% of the variance. Together, PC1 and PC2 explain 96.3% of the variance, suggesting they provide a good representation of the data.
#The labels (IFT88, IFT81, ARL13B, BBS1) represent the original variables.
# P1 (red circles) tend to cluster on the left side, while P4 (blue X's) and P6 (light blue diamonds) are more spread out on the right. This suggests that these groups have different characteristics based on the original variables.
fviz_pca_biplot(pca_result, 
                label = "var", 
                habillage = factor(gene_data$Patient),  # Ensure 'habillage' is a factor
                addEllipses = TRUE, 
                title = "PCA Biplot")



#Clustering - Finding Subgroups in Data
library(cluster)

# Hierarchical clustering
dist_matrix <- dist(gene_data[, -1])  # Compute distance matrix
hc <- hclust(dist_matrix, method = "ward.D2")  # Hierarchical clustering

# Plot dendrogram
plot(hc, labels = gene_data$Patient, main = "Hierarchical Clustering of Patients")


#Advanced Outlier Detection
#Scientists often deal with hidden outliers that affect results.
# Mahalanobis Distance for Multivariate Outliers

library(MVN)

# Compute Mahalanobis distance
m_dist <- mahalanobis(gene_data[, -1], colMeans(gene_data[, -1]), cov(gene_data[, -1]))

# Set a threshold (95% confidence)
threshold <- qchisq(0.95, df = ncol(gene_data) - 1)
# Identify outliers
outliers <- which(m_dist > threshold)
print(outliers)

#Feature Engineering & Transformations
#Raw data often needs transformation before statistical modeling.
# Normalization & Log Transformation
# Min-Max Scaling
gene_data_scaled <- as.data.frame(scale(gene_data[, -1], center = FALSE, scale = apply(gene_data[, -1], 2, max)))

# Log transformation (when data is skewed)
gene_data$GeneA_log <- log1p(gene_data$GeneA)


#Handling Time-Series & Longitudinal Data
#Scientists often collect repeated measures over time (e.g., patient monitoring).
# Working with Time-Series Data

library(lubridate)

# Simulated patient data
time_data <- data.frame(
  Patient = rep(paste0("P", 1:3), each = 4),
  Timepoint = rep(seq(as.Date("2024-01-01"), by = "month", length.out = 4), times = 3),
  BloodPressure = c(120, 122, 125, 123, 130, 128, 132, 135, 110, 115, 118, 120)
)

# Convert to time-series format
time_data$Timepoint <- as.Date(time_data$Timepoint)

# Line plot
ggplot(time_data, aes(x = Timepoint, y = BloodPressure, group = Patient, color = Patient)) +
  geom_line(linewidth = 1) +  # Replace 'size' with 'linewidth'
  geom_point(size = 2) +
  labs(title = "Blood Pressure Over Time", x = "Time", y = "Blood Pressure (mmHg)") +
  theme_minimal()

