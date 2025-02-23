#BENG 628 - Lecture 5: Data Visualization in R with ggplot2
#Lecture Objectives:
#Understand when and why to use different types of data visualizations.
#Gain hands-on experience with ggplot2 for effective data presentation.
#Explore various plot types for categorical and numerical data.
#Customize plots to improve readability and aesthetics.
#Introduction to Data Visualization
#Why visualization is important in life sciences
#Principles of good visualization (clarity, accuracy, storytelling)
#Quick recap of ggplot2 syntax

#Categorical Data Visualizations
#Bar Plots: Showing frequencies and proportions
#Stacked and Grouped Bar Plots
#Pie Charts (when to avoid them)

#Bar Plots for Gene Expression Categories
#Dataset: Differentially Expressed Genes (DEG) in RNA-seq
#Task: Compare the number of upregulated vs. downregulated genes.
#To perform categorical data visualizations (specifically bar plots comparing upregulated vs. downregulated genes) from an RNA-seq dataset like GSE88509, you don't need the normalized count matrices (FPKM or TPM) directly.  Those are for looking at expression levels.  For identifying categories (upregulated/downregulated), you need the results of a differential expression analysis.  Here's the workflow and what you'll need:

#Differential Expression Analysis:

#https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE88509
#o perform categorical data visualizations (specifically bar plots comparing upregulated vs. downregulated genes) from an RNA-seq dataset like GSE88509, you don't need the normalized count matrices (FPKM or TPM) directly.
#Raw Counts Matrix (GSE88509_raw_counts_GRCh38.p13_NCBI.tsv.gz): This is the essential input for differential expression analysis. You must start with raw counts.
#Gene Annotation (Human.GRCh38.p13.annot.tsv.gz): You'll need this to connect the gene IDs in the raw counts matrix to gene names or symbols.
#Experimental Design/Metadata (GSE88509_family.soft.gz): This file (or similar metadata) is absolutely critical. It tells you which samples are controls, which are treatments, and any other relevant experimental conditions. Without this, you can't do differential expression.
#You'll use a tool like DESeq2, edgeR, or limma in R (Bioconductor) to perform the differential expression analysis. 1   These tools take the raw counts, metadata, and gene annotation as input and produce a table of results.   

#Load libraries
library(DESeq2)
library(tidyverse)
library(data.table)
library(GEOquery)
library(edgeR)
library(limma)
library(Biobase)
library(ggplot2)

#let's provide some concrete examples to illustrate those visualization principles
# Clarity
#Good Example:

# Sample data
data <- data.frame(
  Genes = c("ARL13B", "WDR31", "WDR54", "BBS1"),
  Expressiom_level = c(25, 40, 30, 55)
)

# Create a bar plot
#Why it's good: The plot has a clear title, informative axis labels, distinct bars, and a clean background, making it easy to understand the data.
ggplot(data, aes(x = Genes, y = Expressiom_level)) +
  geom_col(fill = "steelblue") +  # Use a clear color
  labs(title = "Gene_Expression",  # Concise and informative title
       x = "Genes",              # Clear axis labels
       y = "Expressiom_level") +
  theme_minimal()                # A clean theme

#Bad Example:
# Create a 3D pie chart (generally avoid these!)
#Why it's bad: 3D pie charts distort the proportions, 
#making it difficult to compare categories accurately. 
#The bright, distracting colors don't add any useful information.
library(plotrix) # Needed for 3D pie charts
pie3D(data$Expressiom_level, labels = data$Genes, explode = 0.1,
      main = "Gene_Expression (3D Pie Chart)",
      col = rainbow(length(data$Genes))) # Too many bright colors


#Accuracy
# Sample data with a clear trend
x <- 1:10
y <- x + rnorm(10, mean = 0, sd = 1)  # Add some noise
data <- data.frame(x, y)

#Good Example
# Create a scatter plot with a trend line
#Why it's good: The scatter plot accurately represents the data points,
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add trend line
  labs(title = "Relationship between X and Y",
       x = "X",
       y = "Y") +
  theme_minimal()

# Same data as above

#Bad Example
# Create a bar plot with manipulated y-axis
#Why it's bad: Manipulating the y-axis scale can create a misleading impression of the data. 
#In this case, the differences between the bars appear larger than they actually are.
ggplot(data, aes(x = factor(x), y = y)) +  # Treat x as a factor
  geom_col() +
  ylim(0, max(y) * 1.5) +  # Exaggerate the y-axis scale
  labs(title = "Relationship between X and Y (Misleading)",
       x = "X",
       y = "Y") +
  theme_minimal()

#Storytelling

#Good Example
# Sample gene expression data
genes <- c("ARL13B", "WDR31", "WDR54", "BBS1")
expression_control <- c(10, 12, 8, 15)
expression_treated <- c(5, 20, 10, 12)
data <- data.frame(Gene = genes, Control = expression_control, Treated = expression_treated)

# Reshape data for plotting

data_long <- data %>%
  tidyr::pivot_longer(cols = c("Control", "Treated"), names_to = "Condition", values_to = "Expression")

# Create a grouped bar plot
#Storytelling: "This plot shows how the expression of four genes changes in response to treatment. 
#WDR31 shows a dramatic increase in expression, while ARL13B is significantly downregulated. 
#This suggests that WDR31 may be involved in the treatment's mechanism of action, while ARL13B might be suppressed.
ggplot(data_long, aes(x = Gene, y = Expression, fill = Condition)) +
  geom_col(position = "dodge") +  # Grouped bars
  labs(title = "Gene Expression Changes After Treatment",
       x = "Gene",
       y = "Expression Level",
       fill = "Condition") +
  theme_minimal()

#Bad Example
# Same data as above
# Create a simple bar plot without grouping
#Why it lacks storytelling: This plot only shows the control group, providing no information about the treatment effect. 
#It doesn't tell a complete story about the data.
ggplot(data, aes(x = Gene, y = Control)) +
  geom_col(fill = "blue") +
  labs(title = "Gene Expression (Control Group)",
       x = "Gene",
       y = "Expression Level") +
  theme_minimal()

#Choosing the Right Plot


#Interactive Visualizations (plotly, ggvis):
library(plotly)
#plotly:A simple scatter plot with tooltips that show data values when hovering over points
plot_ly(data = mtcars, x = ~wt, y = ~mpg, type = 'scatter', mode = 'markers',
        hoverinfo = 'text',
        text = ~paste("Car: ", rownames(mtcars), "<br>Weight: ", wt, "<br>MPG: ", mpg))



# Download and parse the metadata from GEO
#GSE11121: This dataset includes RNA-seq data from breast cancer cell lines treated with different drugs.
gse <- getGEO("GSE11121", GSEMatrix = TRUE)  # Loads the dataset
# Extract the sample information (metadata)
sample_info <- pData(gse[[1]])  # Extract metadata from the first list element
# Extract the count data (expression matrix)
count_data <- exprs(gse_obese[[1]])  # Extract expression data from the first element of the list
#Inspect Column Names to Find Treatment Information
colnames(sample_info)  # Check available metadata columns

#Check characteristics_ch1 Columns for Treatment Info
#Run this command to preview the data in each column
#Look for keywords like treatment, condition, control, experiment, or disease state.
sample_info %>% dplyr::select(starts_with("characteristics_ch1")) %>% head()
#Look at description Columns:
sample_info %>% dplyr::select(starts_with("description")) %>% head()
#the "description" columns mainly contain metadata about the experiment setup, such as sample origin, sequencing details, and replicate numbers.
sample_info %>% dplyr::select(source_name_ch1) %>% head()
sample_info %>% dplyr::select(starts_with("characteristics_ch1")) %>% head()
sample_info %>% dplyr::select(starts_with("relation")) %>% head()
sample_info %>% dplyr::select(starts_with("supplementary_file")) %>% head()

# Correct way to include row names as a column
sample_info_selected <- sample_info %>%
  mutate(Sample = rownames(.)) %>%  # Add row names as a new column called 'Sample'
  dplyr::select(Sample, Condition = characteristics_ch1.6)
sample_info_selected
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121


# Download and parse the metadata from GEO
#GSE183879:RNA sequencing of human skeletal muscle biopsies from individuals with obesity with and without type 2 diabetes.
#Conditions:
#Obese with type 2 diabetes
#Obese without type 2 diabetes
#Lean controls
#Compare gene expression between the three groups to identify genes associated with obesity and/or diabetes.
# Load necessary libraries
library(GEOquery)
library(limma)
library(dplyr)
library(data.table)
library(DESeq2)

# Load the dataset
gse_obese <- getGEO("GSE183879", GSEMatrix = TRUE)

# Extract expression data (first element of the list)
expr_data <- exprs(gse_obese[[1]])

# Verify the structure of the GSE object
str(gse_obese)
head(expr_data)

# Extract phenoData
pheno_data <- pData(gse_obese[[1]])

# Check unique conditions in the column 'characteristics_ch1.2'
unique_conditions <- unique(pheno_data$characteristics_ch1.2)
print(unique_conditions)

# Subset the data for treated samples
treated_samples <- pheno_data[pheno_data$characteristics_ch1.2 == "treatment: Treated with 1 µM G1 for 24 hours", ]
treated_expression <- gse_obese[[1]][, rownames(treated_samples)]

# Check the number of features in the ExpressionSet
dim(treated_expression)

# Download the raw data from GEO FTP
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183879/suppl/GSE183879_all_compare.txt.gz"
download.file(url, destfile = "GSE183879_all_compare.txt.gz")

# Read the raw data into R
raw_data <- fread("GSE183879_all_compare.txt.gz")

# Check the structure of raw_data
head(raw_data)

# Assuming 'raw_data' contains count data in columns 2 to 7 (adjust based on actual column structure)
count_data <- raw_data[, 2:7]

# Check the structure of count data
str(count_data)

# Assuming 'col_data' contains sample metadata (condition, etc.)
col_data <- data.frame(
  row.names = colnames(count_data),
  condition = pheno_data$characteristics_ch1.2  # Adjust this based on actual condition column name
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Extract normalized counts
normalized_data <- counts(dds, normalized = TRUE)

# Check the top 6 differentially expressed genes (if you want)
res <- results(dds)
head(res)
# Remove rows with missing values (NA) for padj or log2FoldChange
res_clean <- as.data.frame(res) %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))

#let's check how many genes are actually significant
# Check number of significant genes
sum(res_clean$padj < 0.05, na.rm = TRUE)

# Look at the range of adjusted p-values
summary(res_clean$padj)
# Look at the distribution of fold changes
summary(res_clean$log2FoldChange)

#Sample size/statistical power:
#How many replicates do you have per condition?
# Check sample sizes
table(col_data$condition)
# Check design
design(dds)

# Look at the contrast being tested
resultsNames(dds)

#If you're not seeing any red points, this suggests one of several possibilities:
#Your differential expression analysis didn't find any significant genes
#There might be an issue with the experimental design or contrast
#The p-value adjustment might be too stringent

#Your differential expression analysis didn't find any significant genes
#The analysis is structured correctly but the small sample size (3 vs 3) is likely making it hard 
#to detect statistically significant differences after multiple testing correction.
# Create the volcano plot again
ggplot(res_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", title = "Volcano plot")

# Generate PCA plot to check sample clustering
vst_data <- vst(dds)
plotPCA(vst_data, intgroup="condition")

# Check sample-to-sample distances
sampleDists <- dist(t(assay(vst_data)))
plot(hclust(sampleDists))

# Look at raw p-values before adjustment
res_raw <- results(dds)
hist(res_raw$pvalue, breaks=50)

#Try a less stringent analysis:
# Use different p-value adjustment method
res_alt <- results(dds, alpha=0.1, pAdjustMethod="BH")
sum(res_alt$padj < 0.1, na.rm=TRUE)

# Look at genes with raw p < 0.05 even if not significant after adjustment
res_raw <- results(dds)
sum(res_raw$pvalue < 0.05, na.rm=TRUE)

# Get top genes by p-value regardless of significance
top_genes <- head(res_raw[order(res_raw$pvalue),], 20)
print(top_genes)

#Modify the visualization to be more informative
# Create a more detailed volcano plot
res_df <- as.data.frame(res_raw)
res_df$significant <- ifelse(res_df$pvalue < 0.05, "raw p<0.05", "not significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.5) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 p-value", 
       title = "Volcano plot (using raw p-values)")


# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

# Create sample dataset for student performance
set.seed(123)
student_data <- data.frame(
  student_id = 1:100,
  math_score = rnorm(100, mean = 75, sd = 10),
  science_score = rnorm(100, mean = 70, sd = 12),
  study_hours = runif(100, 1, 8),
  group = sample(c("A", "B", "C"), 100, replace = TRUE)
)

# Create time series data for stock prices
dates <- seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "day")
stock_data <- data.frame(
  date = dates,
  price = cumsum(rnorm(length(dates), mean = 0.1, sd = 1)) + 100,
  volume = runif(length(dates), 1000, 5000)
)

# Create correlation matrix for heatmap
subjects <- c("Math", "Science", "English", "History", "Art")
correlation_matrix <- matrix(runif(25, -1, 1), nrow = 5)
diag(correlation_matrix) <- 1
colnames(correlation_matrix) <- subjects
rownames(correlation_matrix) <- subjects
correlation_df <- as.data.frame(correlation_matrix) %>%
  rownames_to_column("Subject1") %>%
  pivot_longer(-Subject1, names_to = "Subject2", values_to = "Correlation")

# 1. Histogram
histogram_plot <- ggplot(student_data, aes(x = math_score)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "white") +
  labs(title = "Distribution of Math Scores",
       x = "Math Score",
       y = "Count") +
  theme_minimal()

# 2. Boxplot
boxplot <- ggplot(student_data, aes(x = group, y = math_score, fill = group)) +
  geom_boxplot() +
  labs(title = "Math Scores by Group",
       x = "Group",
       y = "Math Score") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# 3. Violin Plot
violin_plot <- ggplot(student_data, aes(x = group, y = math_score, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Math Score Distribution by Group",
       x = "Group",
       y = "Math Score") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# 4. Scatter Plot with Trend Line
scatter_plot <- ggplot(student_data, aes(x = study_hours, y = math_score)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Relationship between Study Hours and Math Scores",
       x = "Study Hours",
       y = "Math Score") +
  theme_minimal()

# 5. Bubble Plot
bubble_plot <- ggplot(student_data, aes(x = math_score, y = science_score, size = study_hours, color = group)) +
  geom_point(alpha = 0.6) +
  labs(title = "Math vs Science Scores",
       x = "Math Score",
       y = "Science Score",
       size = "Study Hours",
       color = "Group") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

# 6. Time Series Line Plot
line_plot <- ggplot(stock_data, aes(x = date, y = price)) +
  geom_line(color = "steelblue") +
  labs(title = "Stock Price Over Time",
       x = "Date",
       y = "Price") +
  theme_minimal()

# 7. Area Chart
area_plot <- ggplot(stock_data, aes(x = date, y = volume)) +
  geom_area(fill = "lightblue", alpha = 0.5) +
  labs(title = "Trading Volume Over Time",
       x = "Date",
       y = "Volume") +
  theme_minimal()

# 8. Heatmap
heatmap_plot <- ggplot(correlation_df, aes(x = Subject1, y = Subject2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Subject Correlation Heatmap",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Publication-ready adjustments
publication_plot <- scatter_plot +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

# Save high-resolution plot
ggsave("publication_plot.png", publication_plot, 
       width = 8, height = 6, dpi = 300)


# Load necessary libraries
library(GEOquery)  # For downloading GEO datasets
library(ggplot2)   # For visualization
library(dplyr)     # For data manipulation

# Download and load the dataset from GEO (e.g., GSE1000)
gse <- getGEO("GSE1000", GSEMatrix = TRUE)

# Extract the expression data (for example, the first expression set)
expr_data <- exprs(gse[[1]])

# Check the first few rows of the expression data
head(expr_data)

# Convert the first gene's expression values (column) into a data frame
gene_expression <- data.frame(expr_data = expr_data[1, ])


#Histograms can help you understand the distribution of gene expression values for a particular gene or set of genes.
#Histogram of expression for the first gene in the dataset
ggplot(gene_expression, aes(x = expr_data)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  labs(title = "Gene Expression Distribution", x = "Gene Expression", y = "Frequency")

#Boxplots: Comparing Distributions Across Groups
#To compare distributions between different conditions (e.g., control vs treatment), boxplots are idea
# Assuming 'pheno_data' contains sample metadata, and you are comparing two conditions
# Create a data frame where each row corresponds to a sample
# Convert the first gene's expression values (column) into a data frame
# Convert the first gene's expression values into a data frame with proper column names
# Ensure pheno_data contains the correct condition column
pheno_data <- pData(gse[[1]])  # Extract phenotype data
gene_expression_df$Condition <- factor(pheno_data$condition)  # Replace 'condition' with the actual column name
nrow(gene_expression_df)  # Should match the number of samples
nrow(pheno_data)  # Should also match the number of samples
colnames(pheno_data)  # See available column names

# Ensure we use the correct metadata column for conditions
colnames(pheno_data)  # Check available columns
# Example: If "source_name_ch1" contains conditions
pheno_data$Condition <- factor(pheno_data$source_name_ch1)  
# Create a dataframe for the first gene's expression values
gene_expression_df <- data.frame(
  Sample = colnames(expr_data), 
  Expression = expr_data[1, ],
  Condition = pheno_data$Condition  # Ensure correct mapping
)

#Boxplots: Comparing distributions across groups
ggplot(gene_expression_df, aes(x = Condition, y = Expression)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Gene Expression Comparison", x = "Condition", y = "Gene Expression") +
  theme_minimal()

#Violin Plots: A Richer Alternative to Boxplots
ggplot(gene_expression_df, aes(x = Condition, y = Expression)) +
  geom_violin(fill = "lightblue") +
  labs(title = "Gene Expression Distribution", x = "Condition", y = "Gene Expression") +
  theme_classic()

colnames(gene_expression_df)
# Extract expression data for the second gene
second_gene_expression <- data.frame(expr_data = expr_data[2, ])  # Assuming the second gene is in row 2
# Add this second gene expression data to your existing data frame
gene_expression_df$Second_Gene_Expression <- as.numeric(second_gene_expression[1, ])
#Scatter Plots: Exploring Correlations
ggplot(gene_expression_df, aes(x = Expression, y = Second_Gene_Expression)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot of Gene Expression", x = "Gene 1 Expression", y = "Gene 2 Expression") +
  theme_minimal()

# Extract expression data for the second gene
second_gene_expression <- data.frame(expr_data = expr_data[2, ])  # Replace with the correct row index for the second gene

# Add this data to your gene_expression_df
gene_expression_df$Another_Gene_Expression <- as.numeric(second_gene_expression[1, ])

#Adding Trend Lines (geom_smooth)
ggplot(gene_expression_df, aes(x = Expression, y = Another_Gene_Expression)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Scatter Plot with Trend Line", x = "Gene 1 Expression", y = "Gene 2 Expression") +
  theme_classic()

#Bubble Plots: Adding a Third Variable
#To create a bubble plot, you need to scale the size of the points based on a third variable. 
#For example, let's add another feature (e.g., a measure of gene 3 expression).
# Assuming 'Expression3' is the third gene expression
# Assuming 'expr_data' is a matrix or data frame with expression values for multiple genes
gene_expression_df$Expression3 <- expr_data[3, ]  # Replace with actual gene expression data

ggplot(gene_expression_df, aes(x = Expression, y = Another_Gene_Expression, size = Expression3)) +
  geom_point(color = "blue") +
  labs(title = "Bubble Plot of Gene Expression", x = "Gene 1 Expression", y = "Gene 2 Expression") +
  theme_minimal()


#Line Plots: Trends Over Time
# Assuming you have a 'Time' column in your metadata
# (Make sure to replace this with your actual time data)
gene_expression_df$Time <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)  # Example time points

ggplot(gene_expression_df, aes(x = Time, y = Expression)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Gene Expression Over Time", x = "Time", y = "Gene Expression") +
  theme_minimal()

#Area Charts
#If you want to visualize the cumulative sum of expression values over time or any other measure, you can use an area chart
ggplot(gene_expression_df, aes(x = Time, y = Expression)) +
  geom_area(fill = "skyblue", alpha = 0.5) +
  labs(title = "Area Chart of Gene Expression Over Time", x = "Time", y = "Gene Expression") +
  theme_minimal()

#Heatmaps for Multivariate Data
#If you want to visualize multiple genes’ expressions in a heatmap, you can use geom_tile() for each gene across different conditions:
# Assuming you have expression data for multiple genes in a matrix form
# Create a matrix for heatmap (Example with first 5 genes)
# Assuming expr_data is a matrix with gene expressions
# Assuming expr_data is a matrix with gene expression data
heatmap_data <- expr_data[1:5, ]  # Select the first 5 genes for example
colnames(heatmap_data) <- paste("Sample", 1:ncol(heatmap_data), sep = "_")  # Sample names as columns
heatmap_data <- as.data.frame(t(heatmap_data))  # Transpose so that samples are rows

# Convert data to long format
library(reshape2)
heatmap_data_long <- melt(heatmap_data, varnames = c("Sample", "Gene"), value.name = "Expression")

# Now the data frame should have the correct structure: Sample, Gene, Expression
head(heatmap_data_long)  # Check the structure of the long format data

# Plot the heatmap
ggplot(heatmap_data_long, aes(x = Sample, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Heatmap of Gene Expression", x = "Sample", y = "Gene") +
  theme_minimal()



#Heatmap Improvements (pheatmap, ComplexHeatmap):
#Clustering: Show how clustering can reveal patterns.
# Sample data (replace with your actual data)
data <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(data) <- paste0("Gene", 1:10)
colnames(data) <- paste0("Sample", 1:10)

#Explanation: Explain how the clustering algorithm groups similar rows and columns together
# pheatmap example
library(pheatmap)
pheatmap(data, scale="row", clustering_distance_row="euclidean", clustering_method="complete")

# ComplexHeatmap example
library(ComplexHeatmap)
#row_km = 2: This argument specifies that you want to cluster the rows into 2 groups
#column_km = 3: This argument specifies that you want to cluster the columns into 3 groups. You can adjust these numbers based on your data and how many clusters you expect to see
#name = "Expression": This argument sets the title of the legend, making it clear what the colors in the heatmap represent.
#Clustering: The clustering in ComplexHeatmap is hierarchical by default.
#Scaling: If you want to scale your data (e.g., to z-scores) before plotting, you can do that using the scale()
Heatmap(data, row_km = 2,  # Number of row clusters
        column_km = 3,  # Number of column clusters
        name = "Expression" # Title of the legend
) 


#Exercises and Practice:  Add a section with specific exercises for students to apply what they've learned.  This could include:
#Creating different plot types with provided datasets.
#Customizing existing plots.
#Critiquing visualizations from published papers.



#Numerical Data Visualizations

#Histograms: Understanding distributions

#Boxplots: Comparing distributions across groups

#Violin Plots: A richer alternative to boxplots

#Relationship Between Variables

#Scatter Plots: Exploring correlations 

#Adding Trend Lines (geom_smooth)

#Bubble Plots: Adding a third variable

#Time-Series Data Visualization

#Line Plots: Trends over time

#Area Charts

#Heatmaps for Multivariate Data

#Using geom_tile() for heatmaps

#Customizing color scales

#Customizing ggplot2 for Publication-Ready Figures

#Changing themes (theme_minimal(), theme_classic())

#Adjusting labels and legends (labs(), theme())

#Saving high-resolution images (ggsave())






