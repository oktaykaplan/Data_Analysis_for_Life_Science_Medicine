# Load libraries
library(rentrez)
library(Biostrings)
#Search for Genome Accession Numbers
#We first need to get the GenBank accession numbers for C. elegans and Caenorhabditis inopinata.

# Search for genome sequences of C. elegans
celegans_search <- entrez_search(db = "nucleotide", term = "Caenorhabditis elegans[Organism] AND genome", retmax = 5)

# Search for genome sequences of C. inopinata
cinopinata_search <- entrez_search(db = "nucleotide", term = "Caenorhabditis inopinata[Organism] AND genome", retmax = 5)

# Retrieve GenBank accession IDs
celegans_ids <- celegans_search$ids
cinopinata_ids <- cinopinata_search$ids

print(celegans_ids)   # Check IDs
print(cinopinata_ids) # Check IDs

#Fetch FASTA Sequences
#Now, we retrieve the genome sequences in FASTA format using these accession numbers.
# Fetch FASTA sequences for C. elegans
celegans_fasta <- entrez_fetch(db = "nucleotide", id = celegans_ids[1], rettype = "fasta")

# Fetch FASTA sequences for C. inopinata
cinopinata_fasta <- entrez_fetch(db = "nucleotide", id = cinopinata_ids[1], rettype = "fasta")

# Print a preview
cat(substr(celegans_fasta, 1, 500))  # Preview first 500 characters
cat(substr(cinopinata_fasta, 1, 500))

#Save FASTA Files Locally
#To analyze these sequences further, save them as FASTA files.
write(celegans_fasta, file = "celegans_genome.fasta")
write(cinopinata_fasta, file = "cinopinata_genome.fasta")

# Confirm files exist
file.exists("celegans_genome.fasta")
file.exists("cinopinata_genome.fasta")

#Load Sequences into R for Analysis
#Once saved, you can use Biostrings to load and manipulate the sequences.
# Load the FASTA files
celegans_seq <- readDNAStringSet("celegans_genome.fasta")
cinopinata_seq <- readDNAStringSet("cinopinata_genome.fasta")

# Check sequence lengths
length(celegans_seq)
length(cinopinata_seq)

# Print sequence summaries
summary(celegans_seq)
summary(cinopinata_seq)
#Next Steps
#Align the sequences using MUSCLE or msa.
#Construct a phylogenetic tree to compare evolutionary relationships.
#Identify conserved regions between the two species.
#We can extend this with alignment or phylogenetic analysis! 
# Load libraries
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(phangorn)
library(ggseqlogo)
library(rgl)

#Perform Multiple Sequence Alignment (MSA) using MUSCLE
#We will align the sequences using the msa package, which supports MUSCLE.

# Combine sequences into one set
sequences <- DNAStringSet(c(celegans_seq, cinopinata_seq))


# Perform multiple sequence alignment (MSA) using MUSCLE
#Use MUSCLE via Command Line (Faster than R package)
#Try running MUSCLE directly in the terminal instead of using msa(), which can be slow due to R's overhead. 
#Download and install MUSCLE from here. Then, run: muscle -in input.fasta -out output.aln -maxiters 2 -diags
alignment <- msa(sequences, method = "Muscle")
#Run MUSCLE with Lower Iterations in R
#This reduces computation time significantly
alignment <- msa(sequences, method = "Muscle", maxiters = 2)
library(parallel)
options(mc.cores = detectCores())  # Use all available cores



# Print the alignment results
print(alignment)

# Convert to a sequence alignment object compatible with phylogenetic tools
aligned_sequences <- msaConvert(alignment, type = "seqinr::alignment")

#Construct a Phylogenetic Tree
#Now, we use the phangorn and ape packages to build a tree from the alignment.

# Convert alignment to phylogenetic format
phyDat_seq <- phyDat(aligned_sequences, type = "DNA")

# Compute pairwise distance matrix
dist_matrix <- dist.ml(phyDat_seq)

# Build the tree using Neighbor-Joining (NJ) method
tree <- nj(dist_matrix)

# Root the tree at midpoint (optional)
tree <- midpoint(tree)

# Plot the phylogenetic tree
plot(tree, main = "Phylogenetic Tree of C. elegans and C. inopinata")

#For a better visualization, we can use ggtree
# Visualize the phylogenetic tree with ggtree
ggtree(tree) +
  geom_tiplab(size = 4) +  # Add labels
  theme_minimal()

#Identify Conserved Regions
# Compute consensus sequence from alignment
consensus <- consensusString(alignment)

# Print the consensus sequence
print(consensus)

# Identify highly conserved regions (e.g., 80% conservation threshold)
conserved_regions <- which(consensusMatrix(alignment)["A",] >= 0.8 | 
                             consensusMatrix(alignment)["T",] >= 0.8 |
                             consensusMatrix(alignment)["G",] >= 0.8 |
                             consensusMatrix(alignment)["C",] >= 0.8)

# Print conserved positions
print(conserved_regions)

#Save Alignment and Phylogenetic Tree
# Save alignment to FASTA
writeXStringSet(DNAStringSet(alignment), file = "aligned_sequences.fasta")

# Save tree to file
write.tree(tree, file = "phylogenetic_tree.nwk")


# Load data from GEO
library(dplyr)
library(GEOquery)
library(DESeq2)

# 1. Download raw counts from GEO supplementary files
# Example for GSE11121:
#GSE11121 is a microarray dataset from GEO, typically associated with gene expression profiling 
#(in this case, likely using the GPL96 platform, which is Affymetrix Human Genome U133A)
# 2. Get metadata
# Load GSE11121
gse <- getGEO("GSE11121", GSEMatrix = TRUE)

# Access expression data
expr_data <- exprs(gse[[1]])

# Access and inspect phenotype data
pheno <- pData(gse[[1]])
head(pheno[, c("title", "source_name_ch1", "characteristics_ch1")])
print(colnames(pheno))  # Check available columns


unique(metadata$condition)
table(metadata$condition)
#Typically, exprs(gse[[1]]) from GEOquery’s series matrix (e.g., GSE11121_series_matrix.txt.gz) is RMA-normalized and log2-transformed. 
#Your expr_data being raw suggests
head(expr_data[, 1:6])  # First 6 samples
summary(expr_data)      # Range of values

#We need to log2-transform expr_data before running lmFit() to get correct logFC values:
#Log2-Transform expr_data
#log2FC will reflect true log2 fold changes (e.g., -1 = 2-fold down), 
#aligning with microarray norms and biological interpretation.
expr_data_log2 <- log2(expr_data + 1)  # Add 1 to avoid log(0)


#If grade:ch1 Doesn’t Work: If unique(metadata$condition) shows only one value or all NA, try another column:
metadata <- select(pheno, "title", "source_name_ch1", "grade:ch1")
colnames(metadata)[colnames(metadata) == "grade:ch1"] <- "condition"
head(metadata)
unique(metadata$condition)
table(metadata$condition)


# Save the data (optional)
write.table(expr_data, file = "GSE11121_counts.txt", sep = "\t", quote = FALSE)
write.table(metadata, file = "GSE11121_metadata.txt", sep = "\t", quote = FALSE)

# Check row and column names
rownames(metadata) # Should be GSM IDs
colnames(expr_data) # Should be GSM IDs


library(limma)

# Fit linear model with log2-transformed data
fit <- lmFit(expr_data_log2, design)
contrast.matrix <- makeContrasts("grade2-grade1", "grade3-grade1", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Results for grade 2 vs 1
results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
head(results)

#This will show counts for grades "1", "2", "3" (e.g., 50, 100, 50), ensuring the contrast is meaningful.
table(metadata$condition)

library(ggplot2)
#Volcano Plot
degs <- results
colnames(degs)[colnames(degs) == "logFC"] <- "log2FC"
colnames(degs)[colnames(degs) == "adj.P.Val"] <- "padj"

ggplot(degs, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.05, "red", "gray")), alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-value)", title = "Volcano Plot: Grade 2 vs 1") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")

#Verify Transformation (optional):
head(expr_data_log2[, 1:6])  # Should now be ~0–15

#If expr_data isn’t log2-transformed and the logFC values are inflated, we might need to log-transform expr_data before lmFit()
expr_data_log2 <- log2(expr_data + 1)  # Add 1 to avoid log(0)
fit <- lmFit(expr_data_log2, design)


#Epigenetics: DNA Methylation & Expression
#Objective: Investigate relationships between methylation and gene expression.
#Dataset: GSE42865 (cancer vs. normal methylation)
#Load methylation β-values, handle missing data	mice, tidyr
#Identify differentially methylated regions (DMRs)	limma, minfi
#Integrate with TCGA expression data (e.g., UCSCXenaTools)	Data merging, dplyr
#Correlation analysis (methylation vs. expression)	corrplot, ggplot2
#Deliverable: Heatmap of DMRs, scatter plots with trend lines, and pathway analysis.
# Load Required Libraries
library(GEOquery)         # For downloading GSE42865
library(minfi)            # For methylation data processing
library(limma)            # For differential methylation analysis
library(mice)             # For handling missing data
library(tidyr)           # For data tidying
library(dplyr)           # For data manipulation
library(UCSCXenaTools)    # For TCGA expression data
library(corrplot)         # For correlation analysis
library(ggplot2)          # For visualization
library(pheatmap)         # For heatmap
library(missMethyl)       # For pathway analysis
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # Annotation
library(IlluminaHumanMethylation450kmanifest) 
library(ChAMP)            # For additional methylation analysis

# Set working directory (adjust as needed)
setwd("your/path/here")

# --- Load Methylation Data ---

# Download and load GSE42865 from GEO
gse <- getGEO("GSE42865", GSEMatrix = TRUE)
methyl_data <- gse[[1]]

# Extract β-values (methylation levels)
beta_values <- exprs(methyl_data)

# Extract phenotype data
pData <- pData(methyl_data)
ls()
str(sample_groups)
table(sample_groups)
head(pData)  # Check column names
colnames(pData)

# Handle missing data
beta_values_df <- as.data.frame(beta_values)
if (any(is.na(beta_values_df))) {
  imputed_data <- mice(beta_values_df, m = 5, method = "pmm", seed = 123)
  beta_values_complete <- complete(imputed_data)
} else {
  beta_values_complete <- beta_values_df
}


# --- Differential Methylation Analysis ---

# Create sample groups based on 'genotype/variation' column
sample_groups <- factor(ifelse(grepl("LMNA mutation", pData$`genotype/variation:ch1`), "prostate cancer", "normal"))
length(sample_groups)  # Verify the number of samples
table(sample_groups)   # Check the distribution of sample groups

# Design matrix for limma
design <- model.matrix(~ sample_groups)
rownames(design) <- colnames(beta_values_complete)  # Ensure row names match sample names

# Check dimensions to ensure the design matrix and beta_values are aligned
dim(beta_values_complete)
dim(design)

# Ensure that design matches the order of beta values
design <- design[match(colnames(beta_values_complete), rownames(design)), ]

# Fit model using limma and eBayes for differential methylation
fit <- lmFit(beta_values_complete, design)
fit <- eBayes(fit)
#Examine the Data for Variability: Check whether the methylation values in your dataset have substantial variability across samples
# Check for zero variance in your beta values
apply(beta_values_complete, 1, var)  # Apply variance check for each CpG
diff_meth_results <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

# Select significantly differentially methylated CpGs
dm_CpGs <- diff_meth_results %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.5)

# --- Annotate CpGs with Genes ---

# Load annotation data
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Check column names in dm_CpGs
colnames(dm_CpGs)

# Check column names in annotation
colnames(annotation)
# Check for duplicate values in both dataframes
any(duplicated(dm_CpGs$ID))
any(duplicated(annotation$Name))

# Assuming that the CpG identifiers are available in 'annotation$Name', you could add it to dm_CpGs
dm_CpGs$ID <- rownames(dm_CpGs)  # Or use an appropriate column that holds CpG identifiers in 'dm_CpGs'

# Merge CpG data with annotation to get gene information
dm_annotated <- merge(dm_CpGs, annotation, by.x = "ID", by.y = "Name")

# Extract unique gene symbols
dm_genes <- unique(dm_annotated$UCSC_RefGene_Name)

# --- Load Gene Expression Data (TCGA PRAD) ---

# Retrieve TCGA prostate cancer expression data (Ensure the correct dataset is loaded)
library(UCSCXenaTools)
library(dplyr)

# Find the correct dataset
prad_datasets <- XenaData %>%
  filter(XenaCohorts == "TCGA Prostate Adenocarcinoma") %>%
  filter(grepl("RNAseq|HiSeq", DataSubtype))  # Focus on RNA-seq
print(prad_datasets)  # Inspect to confirm

# Example: Use HiSeqV2 RNA-seq dataset
tcga_query <- XenaData %>%
  filter(XenaDatasets == "TCGA.PRAD.sampleMap/HiSeqV2")  # Adjust based on print output

# Download and prepare
tcga_data <- XenaGenerate(tcga_query) %>%
  XenaQuery() %>%
  XenaDownload(destdir = "./tcga_data") %>%  # Specify a directory
  XenaPrepare()

# Load the data
tcga_expr <- tcga_data  # Should be a matrix or data frame
# Check structure
str(tcga_expr)
head(tcga_expr)

# If there's a 'sample' column, remove it or adjust
if ("sample" %in% colnames(tcga_expr)) {
  tcga_expr <- tcga_expr %>% select(-sample)  # Drop 'sample' column
} else if (colnames(tcga_expr)[1] == "X") {
  tcga_expr <- tcga_expr %>% select(-X)  # Drop row index column if present
}

# Ensure row names are genes
rownames(tcga_expr) <- tcga_expr[, 1]  # Assuming first column is genes
tcga_expr <- tcga_expr[, -1]  # Drop gene column from data
head(tcga_expr)  # Should now be samples as columns, genes as rows


library(tidyr)
library(dplyr)
# From tcga_expr
tcga_expr_df <- as.data.frame(tcga_expr)
if ("sample" %in% colnames(tcga_expr_df)) {
  rownames(tcga_expr_df) <- tcga_expr_df$sample
  tcga_expr_df <- tcga_expr_df %>% select(-sample)
}

# Verify gene names are row names
head(tcga_expr_df)

# Convert to data frame and set gene names as row names
# Tidy the data 
tcga_expr_tidy <- tcga_expr_df %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

head(tcga_expr_tidy)

#Process GSE42865 Methylation and Integrate

# Methylation data
beta_values_tidy <- beta_values_complete %>%
  as_tibble() %>%
  mutate(CpG = rownames(beta_values_complete)) %>%
  pivot_longer(-CpG, names_to = "Sample", values_to = "Beta")

# Map CpGs to genes
library(minfi)
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
cpg_to_gene <- as.data.frame(annotation[, c("Name", "UCSC_RefGene_Name")])
colnames(cpg_to_gene) <- c("CpG", "Gene")
cpg_to_gene <- cpg_to_gene %>% filter(!is.na(Gene))

# Aggregate methylation by gene
methyl_gene <- beta_values_tidy %>%
  left_join(cpg_to_gene, by = "CpG") %>%
  group_by(Gene, Sample) %>%
  summarise(Avg_Beta = mean(Beta, na.rm = TRUE))

# Merge with TCGA expression
combined_data <- methyl_gene %>%
  inner_join(tcga_expr_tidy, by = "Gene")

head(combined_data)



