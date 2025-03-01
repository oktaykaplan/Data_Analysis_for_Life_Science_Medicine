# BENG 628 - Lecture 4: Exploratory Data Analysis for Life Sciences
# Enhanced Code with Explanations


# 1. Advanced Data Visualization with ggplot2 ----------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)


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


# Simulated gene expression data 2
gene_kinesin_expression <- data.frame(
  Gene = rep(c("KIF1A", "KIF2A", "KIF3A", "KIF5B", "KIF11", "KIF14", "KIF21A", "KIF22", "BRCA1", "TP53", "MYC"), each = 5),
  Tissue = rep(c("Liver", "Lung", "Brain", "Kidney", "Heart"), times = 11), 
  Expression = sample(5:20, 55, replace = TRUE) 
)

View(gene_kinesin_expression)


# Gene Expression Visualization -------------------------------------------------
# Improved boxplot with better styling and orientation
ggplot(gene_kinesin_expression, aes(x = Gene, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.color = "red", alpha = 0.8) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression Distribution Across Tissues",
       subtitle = "Kinesin Family Genes and Cancer Markers",
       x = "Gene Symbol",
       y = "Expression Level (RPKM)",
       caption = "Simulated data for educational purposes")


#2. Handling Missing Data (DNA Methylation Example) --------------------------------------------
# Improved methylation analysis with multiple imputation strategies
#Objective: Identify and handle missing CpG methylation values
#Dataset: Simulated CpG site methylation levels (some missing values).
# Simulated DNA methylation dataset
# Simulated DNA methylation dataset
methylation <- data.frame(
  CpG_Site = paste0("CpG_", 1:10),
  Methylation_Level = c(85, 78, NA, 67, 90, 72, 88, NA, 80, 75)
)

# Check missing values
print(methylation)
print(sum(is.na(methylation$Methylation_Level)))  # Count missing values

# Imputation methods: mean and median
methylation_enhanced <- methylation %>%
  mutate(
    Mean_imputed = ifelse(is.na(Methylation_Level),
                          mean(Methylation_Level, na.rm = TRUE),
                          Methylation_Level),
    Median_imputed = ifelse(is.na(Methylation_Level),
                            median(Methylation_Level, na.rm = TRUE),
                            Methylation_Level)
  )

# Print original and imputed data
print("Original data with missing values:")
print(methylation, row.names = FALSE)

print("\nImputation results:")
print(methylation_enhanced, row.names = FALSE)

#ChIP-Seq Data Simulation and Correlation with Methylation and Gene Expression
#Now let's generate simulated ChIP-Seq data for binding of a protein to specific genomic regions and calculate its correlation with both methylation and gene expression.
#Simulated ChIP-Seq Data
# Simulated ChIP-Seq data for a few genes
set.seed(123)
chipseq_data <- data.frame(
  Gene = rep(c("KIF1A", "KIF2A", "KIF3A", "KIF5B", "KIF11", "KIF14", "KIF21A", "KIF22", "BRCA1", "TP53", "MYC"), each = 5),
  Tissue = rep(c("Liver", "Lung", "Brain", "Kidney", "Heart"), times = 11),
  ChIP_Signal = runif(55, 0, 10)  # Simulated ChIP-Seq signal
)

# Merge ChIP-Seq data with gene expression data
merged_data <- merge(gene_kinesin_expression, chipseq_data, by = c("Gene", "Tissue"))

# Calculate the correlation between gene expression, methylation, and ChIP-Seq data
# Create a combined data frame with expression, methylation, and ChIP-Seq values
combined_data <- merged_data %>%
  left_join(methylation_enhanced %>% select(CpG_Site, Mean_imputed), by = c("Gene" = "CpG_Site"))

# Calculate the correlation between methylation, expression, and ChIP-Seq
cor_matrix_combined <- cor(combined_data %>% select(Expression, Mean_imputed, ChIP_Signal), method = "pearson")

# Print the correlation matrix
print(cor_matrix_combined)


#Visualizing Correlations with Heatmap
#Let's visualize the correlations between gene expression, methylation, and ChIP-Seq data in a heatmap.
# Reshape the correlation matrix for heatmap visualization
library(reshape2)

cor_melted_combined <- melt(cor_matrix_combined)

# Create the heatmap of correlations
ggplot(cor_melted_combined, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap: Gene Expression, Methylation, and ChIP-Seq",
       x = "Variables",
       y = "Variables",
       fill = "Correlation")

#Step-by-Step Process for Correlation Visualization of Specific Genes:
#Step 1: Filter Data for Specific Genes
# Choose specific genes of interest for the correlation analysis
selected_genes <- c("KIF1A", "BRCA1", "TP53", "MYC")

# Filter the merged data to include only the selected genes
filtered_data <- merged_data %>%
  filter(Gene %in% selected_genes)

# Merge methylation data (imputed values) for the selected genes
filtered_combined_data <- filtered_data %>%
  left_join(methylation_enhanced %>% select(CpG_Site, Mean_imputed), by = c("Gene" = "CpG_Site"))

# View the filtered dataset
print(filtered_combined_data)

#Step 2: Calculate Correlations Between Gene Expression, Methylation, and ChIP-Seq Data
#Calculate the correlation matrix between the expression, methylation, and ChIP-Seq signal for these selected genes.
cor_matrix_filtered <- cor(filtered_combined_data %>% select(Expression, Mean_imputed, ChIP_Signal), method = "pearson")

# Print the correlation matrix for selected genes
print(cor_matrix_filtered)

#Step 3: Visualize the Correlation Matrix in a Heatmap
#Now, visualize the correlation matrix with a heatmap for the selected genes.
# Reshape the correlation matrix for heatmap visualization
library(reshape2)

cor_melted_filtered <- melt(cor_matrix_filtered)

#You will see how KIF1A, BRCA1, TP53, and MYC correlate with methylation and ChIP-Seq data.
#The heatmap will help identify strong positive or negative correlations between these variables.
# Create the heatmap of correlations for selected genes
ggplot(cor_melted_filtered, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap: Gene Expression, Methylation, and ChIP-Seq (Selected Genes)",
       x = "Variables",
       y = "Variables",
       fill = "Correlation")



# 3. Advanced Data Reshaping ---------------------------------------------------
# Enhanced GWAS analysis with multiple conditions

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


gwas_long_expanded <- gwas_long %>%
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "NS"
  ))

ggplot(gwas_long_expanded, aes(x = SNP, y = -log10(p_value), fill = Condition)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = Significance), position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme_classic() +
  labs(title = "GWAS Results Across Experimental Conditions",
       x = "SNP ID", 
       y = "-log10(p-value)")

# 4. Comprehensive Statistical Analysis ----------------------------------------
# Enhanced allele frequency analysis with visualization
#Objective: Analyze the allele frequency distribution of a genetic variant.
#Dataset: Simulated allele frequency data for a genetic variant in a population.
# Simulated allele frequencies
allele_freq <- c(0.12, 0.34, 0.45, 0.29, 0.21, 0.38, 0.25, 0.32, 0.41, 0.37)

allele_df <- data.frame(Frequency = allele_freq)
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


ggplot(allele_df, aes(x = Frequency)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean_allele, color = "red", linetype = "dashed") +
  annotate("text", x = mean_allele + 0.02, y = 3, 
           label = paste("Mean:", round(mean_allele, 2)), color = "red") +
  scale_x_continuous(limits = c(0, 0.5)) +
  labs(title = "Allele Frequency Distribution in Population",
       subtitle = "Simulated Genetic Variant Data",
       x = "Allele Frequency",
       y = "Count") +
  theme_linedraw()


#5.Enhanced HWE Function with Chi-Square Test
hwe_calc_with_test <- function(p, observed_counts = NULL) {
  if(p < 0 | p > 1) stop("Allele frequency must be between 0 and 1")
  
  q <- 1 - p
  AA <- p^2
  Aa <- 2*p*q
  aa <- q^2
  
  expected_counts <- c(AA, Aa, aa) * 1000  # Expected counts based on a sample of 1000
  
  results <- data.frame(
    Genotype = c("AA", "Aa", "aa"),
    Expected_Frequency = c(AA, Aa, aa),
    Expected_Count = round(expected_counts, 2)
  )
  
  # Perform chi-square test if observed counts are provided
  if (!is.null(observed_counts)) {
    if (length(observed_counts) != 3) stop("Observed counts must have exactly 3 values for AA, Aa, and aa")
    
    chisq_test <- chisq.test(observed_counts, p = c(AA, Aa, aa))
    
    results$Observed_Count <- observed_counts
    results$Chi_Square <- round(chisq_test$statistic, 3)
    results$P_Value <- round(chisq_test$p.value, 5)
  }
  
  return(results)
}

# Example usage
observed_data <- c(350, 450, 200)  # Example observed counts for AA, Aa, and aa
hwe_results <- hwe_calc_with_test(0.6, observed_data)

# Display results using knitr
knitr::kable(hwe_results, caption = "Hardy-Weinberg Equilibrium with Chi-Square Test")


# 6. Real-World GEO Data Integration ------------------------------------------
library(GEOquery)
library(tibble)

# Enhanced GEO analysis with full processing pipeline
gse <- getGEO("GSE37745", GSEMatrix = TRUE)
expression_matrix <- exprs(gse[[1]])
pheno_data <- pData(phenoData(gse[[1]]))

# Create annotated data frame
geo_df <- expression_matrix %>%
  as.data.frame() %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
  left_join(pheno_data %>% select(Sample = geo_accession, Tissue = source_name_ch1),
            by = "Sample")
#Ensure Tissue has meaningful label
unique(geo_df$Tissue)

# Visualization of GEO data
geo_df %>%
  group_by(GeneID, Tissue) %>%
  summarize(Mean_Expression = mean(Expression), .groups = "drop") %>%
  slice_max(Mean_Expression, n = 20) %>%
  ggplot(aes(x = reorder(GeneID, Mean_Expression), y = Mean_Expression, fill = Tissue)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(title = "Top Expressed Genes in Kidney Tissue Samples",
       subtitle = "GSE37745 Dataset Analysis",
       x = "Gene ID",
       y = "Normalized Expression Level") +
  theme_minimal(base_size = 12)

#Distribution of Gene Expression Levels
ggplot(geo_df, aes(x = Expression, fill = Tissue)) +
  geom_density(alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(title = "Distribution of Gene Expression",
       x = "Expression Level",
       y = "Density") +
  theme_minimal()

#Principal Component Analysis (PCA)
library(FactoMineR)
library(factoextra)
library(tidyverse)

# Extract sample metadata before transformation
sample_metadata <- geo_df %>%
  select(Sample, Tissue) %>%
  distinct()

# Convert data to wide format
pca_data <- geo_df %>%
  select(-Tissue) %>%  # Remove categorical variable
  pivot_wider(names_from = GeneID, values_from = Expression) %>%
  column_to_rownames("Sample")

# Perform PCA
pca_result <- PCA(pca_data, graph = FALSE)

# Ensure sample metadata is aligned with PCA row names
sample_metadata <- sample_metadata %>%
  filter(Sample %in% rownames(pca_data)) %>%
  arrange(match(Sample, rownames(pca_data)))

# Convert Tissue column to a factor (important for visualization)
sample_metadata$Tissue <- as.factor(sample_metadata$Tissue)

# Check if metadata and PCA rownames match
identical(rownames(pca_data), sample_metadata$Sample)  # Should return TRUE

# Visualize PCA with tissue-based coloring
fviz_pca_ind(pca_result, label = "none", habillage = sample_metadata$Tissue,
             addEllipses = TRUE, ellipse.level = 0.95,
             title = "PCA of Gene Expression Data")

#Hierarchical Clustering and Heatmap
#Identify groups of co-expressed genes or tissue-specific clustering.

library(pheatmap)
library(tidyverse)

# Aggregate data: Mean expression per gene per sample
heatmap_data <- geo_df %>%
  group_by(Sample, GeneID) %>%
  summarize(Average_Expression = mean(Expression), .groups = "drop") %>%
  pivot_wider(names_from = GeneID, values_from = Average_Expression) %>%
  column_to_rownames("Sample")

# Scale the data for better visualization
heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix <- scale(heatmap_matrix)  # Normalize expression values

# Generate heatmap
pheatmap(heatmap_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Gene Expression (GSE37745)")
#Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# Select top 100 most variable genes for a smaller dataset
top_genes <- geo_df %>%
  group_by(GeneID) %>%
  summarize(Variance = var(Expression, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Variance)) %>%
  slice_head(n = 100) %>%
  pull(GeneID)

# Filter for only the top variable genes
heatmap_data_small <- geo_df %>%
  filter(GeneID %in% top_genes) %>%
  group_by(Sample, GeneID) %>%
  summarize(Average_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = GeneID, values_from = Average_Expression) %>%
  column_to_rownames("Sample")

# Convert to matrix
heatmap_matrix_small <- as.matrix(heatmap_data_small)
heatmap_matrix_small <- scale(heatmap_matrix_small)  # Normalize expression

# Generate heatmap with a smaller dataset
pheatmap(heatmap_matrix_small,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Top 100 Variable Genes (GSE37745)")


#Statistical Comparisons
#a) Differential Gene Expression Analysis
#Identify genes that are significantly upregulated or downregulated between tissues.
library(limma)

rm(expression_matrix)  # Remove large data objects
gc()  # Force garbage collection

# Convert to expression matrix
expression_matrix <- geo_df %>%
  pivot_wider(names_from = Sample, values_from = Expression) %>%
  column_to_rownames("GeneID") %>%
  as.matrix()
#Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# Create design matrix for tissue comparison
design <- model.matrix(~ 0 + geo_df$Tissue)
colnames(design) <- levels(geo_df$Tissue)

# Install bigstatsr if not already installed
#Use bigstatsr for Memory Management
#bigstatsr allows you to work with very large datasets without having to load everything into memory.
#install.packages("bigstatsr")
library(bigstatsr)
library(dplyr)
library(tidyr)

# Calculate the variance for each gene
gene_variances <- geo_df %>%
  group_by(GeneID) %>%
  summarize(variance = var(Expression), .groups = 'drop')

# Filter genes with variance above a certain threshold (e.g., 0.5)
filtered_genes <- gene_variances %>%
  filter(variance > 0.5) %>%
  pull(GeneID)

# Subset the original data to only include genes with sufficient variance
filtered_geo_df <- geo_df %>%
  filter(GeneID %in% filtered_genes)

# Aggregate expression values by averaging across samples for each gene
filtered_geo_df_unique <- filtered_geo_df %>%
  group_by(GeneID, Sample) %>%
  summarize(Average_Expression = mean(Expression), .groups = 'drop')

# Pivot to wide format
filtered_matrix <- filtered_geo_df_unique %>%
  pivot_wider(names_from = Sample, values_from = Average_Expression) %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

# Clean Sample and GeneID values
filtered_geo_df$Sample <- trimws(filtered_geo_df$Sample)
filtered_geo_df$GeneID <- trimws(filtered_geo_df$GeneID)

# View the first few rows of the matrix to check the structure
head(filtered_matrix)
# Generate heatmap of expression matrix
pheatmap(filtered_matrix,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = FALSE, 
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Filtered Gene Expression")

#b) Volcano Plot of Differentially Expressed Genes
#Visualize significantly differentially expressed genes.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
# Assuming you already have an expression matrix (filtered_matrix) and a design matrix (design)
# Fit the linear model for differential expression
# Example: assuming you have 6 samples with two conditions (Group1 and Group2)
# This design matrix will describe how the samples are grouped
group <- factor(c("Group1", "Group1", "Group1", "Group2", "Group2", "Group2"))

# Check the column names (samples) in filtered_matrix
colnames(filtered_matrix)
# Check the row names (samples) in design matrix
rownames(design)
# Check the number of columns in the filtered_matrix (samples)
dim(filtered_matrix)

# Assuming the design matrix has the correct number of samples, update the row names to match the sample identifiers
rownames(design) <- colnames(filtered_matrix)[1:6]  # Adjust this based on your specific setup

# Check if the row names of the design matrix now match the column names of filtered_matrix
all(rownames(design) == colnames(filtered_matrix)[1:6])

# Subset filtered_matrix to match the design matrix
filtered_matrix_subset <- filtered_matrix[, rownames(design)]

# Fit the linear model for differential expression
fit <- lmFit(filtered_matrix_subset, design)
# Apply empirical Bayes moderation
fit <- eBayes(fit)

# Check the number of rows in the design matrix (samples)
dim(design)

# Create the design matrix (ensure your design matches your actual sample groupings)
design <- model.matrix(~ group)


# View the design matrix to check
print(design)

# Extract results: logFC (log fold change) and p-values
results <- topTable(fit, coef = 1, n = Inf)  # Replace coef=1 with the appropriate coefficient if needed

# Now you can use EnhancedVolcano
EnhancedVolcano(results,
                lab = rownames(results),
                x = "logFC",
                y = "P.Value",
                title = "Differential Expression Analysis",
                pCutoff = 0.05,
                FCcutoff = 1,
                colAlpha = 0.6,
                legendPosition = "right")

#Functional Enrichment Analysis
#a) Gene Ontology (GO) and KEGG Pathway Analysis
#Check biological functions and pathways enriched in differentially expressed genes.
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
gene_list <- bitr(rownames(results), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Enrichment
go_results <- enrichGO(gene = gene_list$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr")
dotplot(go_results)

# KEGG Pathway Analysis
kegg_results <- enrichKEGG(gene = gene_list$ENTREZID, organism = "hsa", pAdjustMethod = "fdr")
dotplot(kegg_results)

#Machine Learning Approaches
#a) Classification of Tissues Using Gene Expression
library(caret)
library(randomForest)

# Prepare data
train_data <- geo_df %>%
  pivot_wider(names_from = GeneID, values_from = Expression) %>%
  mutate(Tissue = as.factor(Tissue))

# Train random forest classifier
rf_model <- randomForest(Tissue ~ ., data = train_data, ntree = 500)
print(rf_model)

# Feature importance
varImpPlot(rf_model)

#b) Clustering Analysis to Identify Sample Subtypes
#Use k-means clustering to explore subgroups in gene expression patterns.
set.seed(123)
kmeans_result <- kmeans(expression_matrix, centers = 3)
fviz_cluster(kmeans_result, data = expression_matrix)


# 9. Single-Cell RNA-Seq Analysis ---------------------------------------------
library(Seurat)
library(patchwork)

# Enhanced single-cell analysis workflow
pbmc <- NormalizeData(pbmc) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# Visualization panel
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
p2 <- FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E"))
p1 + p2 + plot_annotation(title = "PBMC Single-Cell Clustering",
                          subtitle = "10x Genomics Dataset Analysis")

# 10. Interactive HTML Report -------------------------------------------------
library(plotly)

# Convert ggplot to interactive plot
ggplotly(
  gene_kinesin_expression %>%
    ggplot(aes(x = Tissue, y = Expression, color = Gene)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45))
)

#Here are real-world datasets and code implementations for our tasks using publicly available multi-omics data. 
#I'll use data from TCGA (The Cancer Genome Atlas), ENCODE, and GEO

options(java.parameters = "-Xmx12g")  # Allow Java-based tools to use 12GB RAM
Sys.setenv(R_MAX_VSIZE = 50 * 1024^3)  # Increase max vector memory to 50GB
#Check Available Memory
gc()  # Run garbage collection to free memory
library(TCGAbiolinks)
library(DESeq2)

# Query for RNA-seq data (BRCA: Breast Cancer), limiting by sample type
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"  # Limit to a specific sample type
)

# Download the subset of data
GDCdownload(query)

# Prepare the data after downloading
data <- GDCprepare(query)

# Check the data (this should now be smaller)
head(data)


# Download data in smaller batches (you can specify smaller subsets based on conditions)
GDCdownload(query)

#Since our query retrieved different types of gene expression measurements (e.g., FPKM, TPM, etc.), we need to choose the appropriate assay. 
# Select the appropriate assay (e.g., FPKM values)
expression_data <- assay(data, "fpkm_unstrand")

# Select raw count data (if available)
expression_data <- assay(data, "unstranded")

# Check for missing values in the `paper_BRCA_Subtype_PAM50` column
sum(is.na(colData(data)$paper_BRCA_Subtype_PAM50))

# Remove rows with missing PAM50 subtypes from colData and the corresponding count data
colData_filtered <- colData(data)[!is.na(colData(data)$paper_BRCA_Subtype_PAM50),]
# Convert `paper_BRCA_Subtype_PAM50` to a factor
colData_filtered$paper_BRCA_Subtype_PAM50 <- factor(colData_filtered$paper_BRCA_Subtype_PAM50)

# Prepare DESeqDataSet using filtered and factorized data
dds <- DESeqDataSetFromMatrix(
  countData = expression_data_filtered,  # Filtered expression data
  colData = colData_filtered,            # Filtered sample information
  design = ~ paper_BRCA_Subtype_PAM50    # Use PAM50 subtypes for differential expression
)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Prepare DESeqDataSet using raw counts
dds <- DESeqDataSetFromMatrix(
  countData = expression_data,  # Raw count data
  colData = colData(data),      # Sample information (e.g., clinical data)
  design = ~ paper_BRCA_Subtype_PAM50  # Use PAM50 subtypes for differential expression
)

# Perform DESeq2 analysis
dds <- DESeq(dds)

#Check the results of the DESeq2 analysis
# Check results for differential expression
res <- results(dds)
head(res)

library(EnhancedVolcano)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'pvalue',
  title = 'Differential Expression - PAM50 Subtypes',
  pCutoff = 1e-10,
  FCcutoff = 2
)

