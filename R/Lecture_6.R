#Rare Diseases: Uncovering Mutations in Genetic Disorders
#Case Study: Identifying Novel Genes in Rare Genetic Disorders
#Background: Rare genetic disorders often have limited diagnostic options due to the unknown functions of many involved genes. RNA-seq can provide insights into gene expression patterns and help uncover new mutations.
#RNA-seq Application: Researchers performed RNA-seq on patient samples with an undiagnosed rare genetic disorder. Through the analysis, they were able to identify novel mutations and uncover how these mutations lead to disruptions in normal gene expression and cellular function.
#Outcome: This approach led to the identification of a new genetic mutation causing a rare disorder. Furthermore, the study provided insights into the molecular mechanisms at play and suggested potential treatment options, such as gene therapy.
#Key Takeaways:
#RNA-seq can be a powerful tool for uncovering the molecular basis of rare and previously undiagnosed genetic disorders.
#It can help identify novel mutations and suggest potential therapeutic approaches.

# Load libraries
library(DESeq2)
library(airway)
library(ggplot2)
library(biomaRt)

# Load and Prepare Real RNA-seq Data
data("airway")
counts <- assay(airway)  # Gene expression counts
colData <- colData(airway)  # Sample metadata

# Subset to 6 samples (3 patients, 3 controls)
set.seed(123)  # For reproducibility
sample_subset <- sample(colnames(counts), 6)
counts_subset <- counts[, sample_subset]
colData_subset <- colData[sample_subset, ]

# Simulate patient vs. control groups
colData_subset$condition <- factor(c("patient", "patient", "patient", "control", "control", "control"),
                                   levels = c("control", "patient"))

# View the metadata
print(colData_subset)

# Differential Expression Analysis with DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                              colData = colData_subset,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "patient", "control"))
res <- res[order(res$padj), ]
head(res)

# Filter Significant Overexpressed Genes
sig_genes <- subset(res, padj < 0.05 & log2FoldChange > 1)
cat("Number of overexpressed genes:", nrow(sig_genes), "\n")

# Convert ENSG to HGNC using biomaRt
# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get HGNC symbols for significant genes
ensg_ids <- rownames(sig_genes)
gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensg_ids,
                   mart = ensembl)

# Merge HGNC symbols with results
sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df$ensembl_gene_id <- rownames(sig_genes_df)
sig_genes_df <- merge(sig_genes_df, gene_info, by = "ensembl_gene_id", all.x = TRUE)

# Handle cases where HGNC is missing (keep ENSG as fallback)
sig_genes_df$hgnc_symbol <- ifelse(is.na(sig_genes_df$hgnc_symbol) | sig_genes_df$hgnc_symbol == "",
                                   sig_genes_df$ensembl_gene_id, sig_genes_df$hgnc_symbol)

# Prioritize Genes
sig_genes_df$score <- sig_genes_df$log2FoldChange * -log10(sig_genes_df$padj)
sig_genes_df <- sig_genes_df[order(-sig_genes_df$score), ]

# View top 5 prioritized genes with HGNC symbols
cat("Top 5 prioritized overexpressed genes:\n")
print(head(sig_genes_df[, c("hgnc_symbol", "ensembl_gene_id", "log2FoldChange", "padj", "score")], 5))


# Visualize Results with a Volcano Plot
res_df <- as.data.frame(res)
res_df$ensembl_gene_id <- rownames(res_df)
res_df <- merge(res_df, gene_info, by = "ensembl_gene_id", all.x = TRUE)
res_df$hgnc_symbol <- ifelse(is.na(res_df$hgnc_symbol) | res_df$hgnc_symbol == "",
                             res_df$ensembl_gene_id, res_df$hgnc_symbol)
res_df$significant <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Overexpressed", "Not Significant")
top_genes <- head(sig_genes_df$ensembl_gene_id, 15)
res_df$color <- ifelse(res_df$ensembl_gene_id %in% top_genes, "Top Priority", res_df$significant)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red", "blue"), 
                     labels = c("Not Significant", "Overexpressed", "Top Priority")) +
  geom_text(data = subset(res_df, ensembl_gene_id %in% top_genes), 
            aes(label = hgnc_symbol), 
            vjust = -0.5, hjust = 0.5, size = 3, color = "black") +
  labs(title = "Volcano Plot: Patient vs. Control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  guides(color = guide_legend(title = "Gene Status"))

# Case Study Outcome Simulation
cat("\nCase Study Outcome:\n")
cat("Students identified", nrow(sig_genes), "overexpressed genes in patient samples.\n")
cat("After prioritization, the top gene (", sig_genes_df$hgnc_symbol[1], 
    ", ", sig_genes_df$ensembl_gene_id[1], ") showed a log2 fold change of", 
    round(sig_genes_df$log2FoldChange[1], 2), "and adjusted p-value of", 
    sig_genes_df$padj[1], ".\n")
cat("This suggests a potential novel mutation driving overexpression, disrupting normal cellular function.\n")
cat("Proposed treatment: Gene therapy to silence the overexpressed gene.\n")


#Provide Patient Data: Please provide the gene expression data for the patients from 
#Cummings et al. (2017) in a similar format (e.g., a matrix or data frame with gene IDs and expression values).
library(recount3)
library(DESeq2)
library(preprocessCore)
packageVersion("recount3")

#Manually check available projects
recount3::available_projects()
#gtex is no longer listed as an available project in recount3.

library(recount3)
options(recount3_url = "https://sciserver.org/public-data/recount3/data")  # Alternative URL
human_projects <- available_projects()
gtex_rse <- create_rse(subset(human_projects, project == "MUSCLE_SKELETAL" & project_home == "data_sources/gtex"))


# Define the GTEx IDs (You should have this vector ready from the previous steps)
#All the GTEx IDs were provided in a more manageable format, split into three chunks
gtex_ids <- c("G16537.GTEX-OHPL-1626.2", "G16558.GTEX-P4PP-1626.2", "G16582.GTEX-OHPM-1626.3", "G16596.GTEX-OIZH-1626.2", "G16598.GTEX-P4PQ-1626.2", "G16687.GTEX-OOBJ-1626.3", "G18069.GTEX-OHPK-1626-SM-2YUN3.1", "G18373.GTEX-QCQG-2126-SM-2S1P8.1", "G18402.GTEX-Q2AH-1826-SM-2S1Q2.1", "G18404.GTEX-QDVN-2426-SM-2S1Q4.1", "G18447.GTEX-QV44-2026-SM-2S1RD.1", "G19326.GTEX-S32W-2326-SM-2XCAW.1", "G19355.GTEX-S7SF-2026-SM-3K2AS.1", "G19468.GTEX-POYW-0526-SM-2XCEY.2", "G19491.GTEX-SNMC-1426-SM-2XCFM.2", "G20855.GTEX-O5YT-1626-SM-32PK6.1", "G20912.GTEX-SUCS-1626-SM-32PLS.2", "G21042.GTEX-TML8-1826-SM-32QOR.1", "G21087.GTEX-TKQ2-0826-SM-33HB6.1", "G21092.GTEX-TMMY-0426-SM-33HBB.1", "G22936.GTEX-U4B1-1626-SM-3DB8N.1", "G22959.GTEX-UJHI-1726-SM-3DB9B.1", "G22964.GTEX-U3ZM-1226-SM-3DB9G.1", "G23319.GTEX-R55D-0626-SM-3GAD5.1", "G23324.GTEX-Q734-2026-SM-3GADA.1", "G23343.GTEX-UJMC-1826-SM-3GADT.1", "G23350.GTEX-T5JW-1826-SM-3GAE1.1", "G23437.GTEX-R53T-1826-SM-3GIJX.1", "G25647.GTEX-WRHK-1626-SM-3MJFH.1", "G25653.GTEX-WRHU-0826-SM-3MJFN.1", "G25657.GTEX-WOFM-1326-SM-3MJFR.1", "G25900.GTEX-POMQ-1926-SM-3NB1Y.2", "G25923.GTEX-WY7C-2526-SM-3NB2N.2", "G25950.GTEX-WXYG-2526-SM-3NB3F.2", "G26541.GTEX-WZTO-0826-SM-3NM8Q.1", "G26590.GTEX-T5JC-0626-SM-3NMA6.1", "G26643.GTEX-WHPG-2226-SM-3NMBO.1", "G29343.GTEX-WHSB-1826-SM-3TW8M.1", "G32952.GTEX-U3ZG-0326-SM-47JXN.3", "G32956.GTEX-X638-0326-SM-47JY1.3", "G32988.GTEX-X88G-0326-SM-47JZ4.3", "G33005.GTEX-XUYS-0326-SM-47JX2.3", "G34336.GTEX-XYKS-2426-SM-4AT43.2", "G34379.GTEX-XOTO-0526-SM-4B662.3", "G34401.GTEX-XPT6-2026-SM-4B64V.3", "G34459.GTEX-XQ8I-0626-SM-4BOPT.7", "G34487.GTEX-XUW1-0826-SM-4BOP6.7", "G34523.GTEX-XUJ4-2626-SM-4BOQ3.8", "G34707.GTEX-XUZC-2126-SM-4BRW8.1", "G34751.GTEX-XV7Q-2926-SM-4BRUL.1", "G35292.GTEX-U3ZH-1926-SM-4DXTR.7", "G35593.GTEX-X4XY-0626-SM-4E3IN.4", "G35606.GTEX-XBED-2626-SM-4E3J5.4", "G40775.GTEX-VUSG-2626-SM-4KKZI.2", "G42368.GTEX-XBEC-0626.2", "G42376.GTEX-XBEW-1026.2", "G46557.GTEX-Y114-2526.2", "G46837.GTEX-Y8LW-2026.1", "G46924.GTEX-Y5V5-2526.1", "G46925.GTEX-Y8E4-1026.1", "G47050.GTEX-Y5LM-2126.3", "G49116.GTEX-Y3IK-2626.1", "G49331.GTEX-YEC3-2126.2", "G49368.GTEX-ZV6S-2126.2", "G49499.GTEX-ZQUD-1726.2", "G49500.GTEX-ZQG8-1226.2", "G49505.GTEX-ZTX8-1626.2", "G52281.GTEX-Y8E5-0326.1", "G52338.GTEX-ZVZP-2526.1", "G52366.GTEX-ZPIC-2526.1", "G52382.GTEX-ZPCL-2026.1", "G52423.GTEX-ZT9X-1826.1", "G52456.GTEX-111YS-2326.3", "G52511.GTEX-11EMC-2626.3", "G52613.GTEX-1211K-2126.2", "G52634.GTEX-YBZK-0326.2", "G53040.GTEX-11I78-2426.3", "G53845.GTEX-131XF-2326.2") 
gtex_id1s <- c("G56104.GTEX-11VI4-1926.2", "G56118.GTEX-11WQC-2626.2", "G56120.GTEX-11XUK-2226.3", "G56168.GTEX-Z9EW-1726.2", "G56219.GTEX-12ZZX-0326.2", "G56346.GTEX-ZYWO-2626.3", "G56350.GTEX-ZZ64-1526.3", "G58300.GTEX-113JC-2726.2", "G58373.GTEX-ZAKK-0326.2", "G58504.GTEX-117YX-2526.2", "G58531.GTEX-12ZZY-0626.2", "G59139.GTEX-11ZTT-2626.6", "G59232.GTEX-11ZVC-2726.3", "G59253.GTEX-12BJ1-2526.3", "G59256.GTEX-12C56-1926.3", "G59396.GTEX-ZYFG-2426.4", "G59469.GTEX-1122O-2426.4", "G59538.GTEX-11WQK-0726.2", "G59578.GTEX-11P81-2526.2", "G59646.GTEX-ZY6K-2026.3", "G59671.GTEX-111CU-2026.3", "G59754.GTEX-11EM3-2126.2", "G59774.GTEX-ZVZO-0326.2", "G59785.GTEX-11NSD-2026.2", "G59863.GTEX-ZDYS-1726.2", "G59944.GTEX-131XG-2326.2", "G59960.GTEX-ZC5H-0326.2", "G59995.GTEX-1399R-2526.1", "G60093.GTEX-13JUV-2326.2", "G60097.GTEX-13FHO-0726.2", "G60107.GTEX-13NZ9-0626.2", "G60109.GTEX-139YR-2526.2", "G60125.GTEX-13FTY-0226.2", "G60145.GTEX-13PL7-0626.2", "G60147.GTEX-13U4I-1826.2", "G60275.GTEX-13CF3-1826.2", "G60309.GTEX-13N2G-2326.2", "G60340.GTEX-13O61-2326.2", "G60380.GTEX-1399Q-2426.3", "G60418.GTEX-13NZB-2626.3", "G60480.GTEX-139UW-2626.2", "G60481.GTEX-13FXS-0326.2", "G60551.GTEX-133LE-2026.2", "G60589.GTEX-13FTW-2326.2", "G60652.GTEX-1339X-2426.2", "G60659.GTEX-13D11-2526.2", "G60667.GTEX-132NY-0726.2", "G60681.GTEX-1399S-2726.2", "G60704.GTEX-13OVG-2126.2", "G60844.GTEX-1399U-2526.2", "G60873.GTEX-13FH7-2126.2", "G60893.GTEX-13N11-2726.2", "G60911.GTEX-13OVH-0626.2", "G60925.GTEX-13OVI-1726.2", "G60943.GTEX-13OW6-0626.2", "G61036.GTEX-13VXT-0326.2", "G61055.GTEX-13QBU-2426.2", "G61092.GTEX-144GM-2026.2", "G61098.GTEX-13W46-0726.2", "G61114.GTEX-144GN-2426.2", "G61138.GTEX-YB5E-2226.3", "G61143.GTEX-YB5K-2326.3", "G61154.GTEX-YEC4-2226.3", "G61166.GTEX-Y5V6-2626.3", "G61185.GTEX-YF7O-2526.3", "G61197.GTEX-Y9LG-1926.3")
gtex_id12s <- c("G61849.GTEX-145LV-2326.2", "G61884.GTEX-13111-2226.2", "G61944.GTEX-12WSJ-1726.2", "G61950.GTEX-12WSN-2526.2", "G61955.GTEX-1314G-1726.2", "G62186.GTEX-145ME-2026.2", "G62216.GTEX-1497J-2626.2", "G62228.GTEX-13W3W-2626.2", "G62275.GTEX-11DXZ-2426.2", "G62501.GTEX-ZV7C-2426.2", "G62593.GTEX-146FQ-0326.2", "G62605.GTEX-145MN-2426.2", "G62612.GTEX-147F3-0226.2", "G62660.GTEX-ZTPG-0126.2", "G62709.GTEX-145LT-1626.2", "G62774.GTEX-13YAN-0526.3", "G62899.GTEX-139D8-0726.3", "G63046.GTEX-11EQ9-2126.2", "G63072.GTEX-132AR-1")


# Get available human projects
human_projects <- available_projects()
# Check the first few rows of human_projects
head(human_projects)
# Look at unique project names
unique(human_projects$project)

# Get the correct URL for GTEx data
gtex_url <- "https://recount-opendata.s3.amazonaws.com/recount3/release/human/gtex"

# Get GTEx projects
gtex_projects <- available_projects(recount3_url = gtex_url)

# View the GTEx projects
head(gtex_projects)


# Access GTEx data specifically
project_homes <- available_homes("human")
gtex_projects <- available_projects(recount3_url = project_homes["gtex"])
# Check the structure of gtex_info
dim(gtex_info)
head(gtex_info)

# We need to specify a single tissue/study
# Let's see what options are available
unique(gtex_info$file_source)



# Create the RSE object
gtex_project <- create_rse(gtex_info)

# Define your GTEx IDs (combine all lists)
gtex_all_ids <- c(gtex_ids, gtex_id1s, gtex_id12s)

# Access GTEx project (assuming skeletal muscle data)
gtex_project <- create_rse(project = "GTEX")

# Extract metadata to match your IDs
metadata <- colData(gtex_project)
available_ids <- metadata$external_id  # Adjust column name based on recount3 metadata

# Filter for your specific IDs
matched_ids <- gtex_all_ids[gtex_all_ids %in% available_ids]
if (length(matched_ids) < length(gtex_all_ids)) {
  warning("Some IDs not found in recount3 GTEx data!")
  print(setdiff(gtex_all_ids, matched_ids))
}

# Subset the RangedSummarizedExperiment object
gtex_subset <- gtex_project[, matched_ids]

# Extract gene expression counts
gene_counts <- assays(gtex_subset)$counts

# Save to file
write.csv(as.data.frame(gene_counts), "gtex_subset_gene_counts.csv", row.names = TRUE)

# Convert to DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = assays(rse_gene_subset)$counts, 
  colData = colData(rse_gene_subset), 
  design = ~1
)

# Normalize the counts using quantile normalization
norm_counts <- normalize.quantiles(assay(dds))

# Convert the normalized data back to a matrix
assay(dds) <- norm_counts

# Run DESeq2 analysis (if needed)
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Display summary of the results
summary(res)


#################################################Facioscapulohumeral Muscular Dystrophy (FSHD):FSHD is a rare neuromuscular disorder with complex genetic mechanisms.############################
#Studies have used RNA-seq to investigate the role of DUX4 gene expression in FSHD pathogenesis.
#You can find RNA-seq data related to FSHD on NCBI GEO. Search for "FSHD RNA-seq" and look for studies with patient and control samples.
#Facioscapulohumeral muscular dystrophy (FSHD) is a genetic muscle disorder linked to the aberrant expression of the DUX4 gene. Research has identified that DUX4 resides within the D4Z4 macrosatellite repeat array on chromosome 4q35. In FSHD, a reduction in the number of these repeats leads to the inappropriate expression of the DUX4 protein, which is toxic to muscle cells. 
#To investigate the role of DUX4 in FSHD pathogenesis, several studies have utilized RNA sequencing (RNA-seq) to compare gene expression profiles between FSHD patients and healthy controls. These studies aim to elucidate the molecular mechanisms underlying FSHD and identify potential therapeutic targets.
#The National Center for Biotechnology Information's Gene Expression Omnibus (NCBI GEO) is a valuable resource for accessing such datasets. By searching for "FSHD RNA-seq" on NCBI GEO, you can find studies that have employed RNA-seq to analyze gene expression in FSHD patient samples compared to controls. These datasets can provide insights into the differential expression of genes, including DUX4, and help further our understanding of FSHD pathogenesis.
#Download Manually 
#Go to the GSE116710 page on NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108426
#GSE108426 and GSE184643

options(timeout = 600)
#Visit the GSE174301 page on the NCBI GEO website: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174301
#Download the GSE174301_RAW.tar file (36.3 Mb).
#This TAR archive likely contains the processed count data files (RSEM output) that you need.
#The GEO page indicates "Processed data provided as supplementary file", this is where the count data will be.
#This is the best place to get the data, as it is the processed count data.

# Load required libraries
library(DESeq2)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)

# Set working directory - modify this to your path
setwd("/home/kaplanlab/Downloads")

# Extract tar file if not already extracted
if (!dir.exists("extracted_files")) {
  untar("GSE174301_RAW.tar", exdir = "extracted_files")
}

# List all RSEM gene results files
file_list <- list.files("extracted_files", pattern = "*rsem.genes.results.gz", full.names = TRUE)
print(paste("Found", length(file_list), "files"))

# Extract sample names from file paths for better labeling
sample_names <- gsub(".*/(.*?)_rsem\\.genes\\.results\\.gz", "\\1", file_list)
print("All sample names:")
print(sample_names)

# Create an empty matrix to store expected counts
gene_ids <- read.table(gzfile(file_list[1]), header = TRUE, row.names = 1, sep = "\t")
merged_counts <- matrix(0, nrow = nrow(gene_ids), ncol = length(file_list))
rownames(merged_counts) <- rownames(gene_ids)
colnames(merged_counts) <- sample_names

# Loop through each file and extract the expected counts
for (i in 1:length(file_list)) {
  data <- read.table(gzfile(file_list[i]), header = TRUE, row.names = 1, sep = "\t")
  # Extract only the 'expected_count' column and round to integers for DESeq2
  merged_counts[, i] <- round(data$expected_count)
}

# Filter low-expression genes (remove genes with zero counts in more than 50% of samples)
threshold <- ncol(merged_counts) * 0.5
filtered_counts <- merged_counts[rowSums(merged_counts > 0) > threshold, ]
print(paste("Genes before filtering:", nrow(merged_counts)))
print(paste("Genes after filtering low expression:", nrow(filtered_counts)))

# Inspect sample names more carefully to understand the pattern
print("Sample name patterns analysis:")
sample_patterns <- list()
for (name in sample_names) {
  # Extract potential condition identifiers
  if (grepl("Control", name)) sample_patterns$Control <- c(sample_patterns$Control, name)
  if (grepl("FSHD1", name)) sample_patterns$FSHD1 <- c(sample_patterns$FSHD1, name)
  if (grepl("FSHD2", name)) sample_patterns$FSHD2 <- c(sample_patterns$FSHD2, name)
  
  # Extract potential time points
  for (day in c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20)) {
    day_pattern <- paste0("Day", day)
    if (grepl(day_pattern, name)) {
      if (is.null(sample_patterns[[day_pattern]])) sample_patterns[[day_pattern]] <- character(0)
      sample_patterns[[day_pattern]] <- c(sample_patterns[[day_pattern]], name)
    }
  }
}

print("Detected patterns in sample names:")
for (pattern in names(sample_patterns)) {
  print(paste(pattern, ":", length(sample_patterns[[pattern]]), "samples"))
}

# Create a proper metadata table based on the analysis
create_metadata <- function(sample_names) {
  # Initialize vectors
  condition <- character(length(sample_names))
  time_point <- character(length(sample_names))
  
  # Loop through sample names and extract information
  for (i in 1:length(sample_names)) {
    name <- sample_names[i]
    
    # Extract condition
    if (grepl("Control", name)) {
      condition[i] <- "Control"
    } else if (grepl("FSHD1", name)) {
      condition[i] <- "FSHD1"
    } else if (grepl("FSHD2", name)) {
      condition[i] <- "FSHD2"
    } else {
      condition[i] <- "Unknown"
    }
    
    # Extract time point using regex
    day_match <- regexpr("Day\\d+", name)
    if (day_match > 0) {
      time_point[i] <- substr(name, day_match, day_match + attr(day_match, "match.length") - 1)
    } else {
      time_point[i] <- "Unknown"
    }
  }
  
  # Create and return metadata
  metadata <- data.frame(
    sample = sample_names,
    condition = condition,
    time = time_point,
    stringsAsFactors = FALSE,
    row.names = sample_names
  )
  
  return(metadata)
}

# Generate metadata
metadata <- create_metadata(sample_names)
print("Generated metadata (first 10 rows):")
print(head(metadata, 10))

# Count samples per condition and time point
print("Sample counts by condition:")
print(table(metadata$condition))
print("Sample counts by time point:")
print(table(metadata$time))
print("Sample counts by condition and time:")
print(table(metadata$condition, metadata$time))

# Based on the counts, select the most appropriate comparison
# For example, if we have enough Control and FSHD2 samples at a specific time point
condition_time_counts <- table(metadata$condition, metadata$time)
print("Detailed condition-time combinations:")
print(condition_time_counts)

# Find the time point with most samples for Control vs FSHD2 comparison
control_fshd2_counts <- condition_time_counts[c("Control", "FSHD2"), ]
if ("Control" %in% rownames(control_fshd2_counts) && "FSHD2" %in% rownames(control_fshd2_counts)) {
  # Find time points with at least 3 samples in each group
  valid_time_points <- apply(control_fshd2_counts, 2, function(x) all(x >= 3))
  if (any(valid_time_points)) {
    target_time <- names(valid_time_points)[which(valid_time_points)[1]]
    target_cond1 <- "Control"
    target_cond2 <- "FSHD2"
    print(paste("Selected comparison:", target_cond2, "vs", target_cond1, "at", target_time))
  } else {
    # If Control vs FSHD2 doesn't have enough samples, try Control vs FSHD1
    control_fshd1_counts <- condition_time_counts[c("Control", "FSHD1"), ]
    if ("Control" %in% rownames(control_fshd1_counts) && "FSHD1" %in% rownames(control_fshd1_counts)) {
      valid_time_points <- apply(control_fshd1_counts, 2, function(x) all(x >= 3))
      if (any(valid_time_points)) {
        target_time <- names(valid_time_points)[which(valid_time_points)[1]]
        target_cond1 <- "Control"
        target_cond2 <- "FSHD1"
        print(paste("Selected comparison:", target_cond2, "vs", target_cond1, "at", target_time))
      } else {
        # If no obvious comparison, select the best available
        print("Warning: No ideal comparison found. Using the most samples available:")
        # Find the condition-time combination with most samples
        max_samples <- which(condition_time_counts == max(condition_time_counts), arr.ind = TRUE)
        if (nrow(max_samples) > 0) {
          row_idx <- max_samples[1, 1]
          col_idx <- max_samples[1, 2]
          target_cond1 <- rownames(condition_time_counts)[row_idx]
          target_time <- colnames(condition_time_counts)[col_idx]
          # Find another condition at the same time point
          other_conds <- rownames(condition_time_counts)[-row_idx]
          other_cond_counts <- condition_time_counts[other_conds, col_idx]
          if (length(other_cond_counts) > 0 && max(other_cond_counts) >= 3) {
            target_cond2 <- names(which.max(other_cond_counts))
          } else {
            # If necessary, use the second best condition
            target_cond2 <- names(sort(other_cond_counts, decreasing = TRUE)[1])
          }
          print(paste("Selected comparison:", target_cond2, "vs", target_cond1, "at", target_time))
        }
      }
    }
  }
} else {
  # Manual selection if needed
  print("Warning: Could not automatically select comparison. Using manual selection:")
  # Choose the most frequent condition
  cond_counts <- sort(table(metadata$condition), decreasing = TRUE)
  if (length(cond_counts) >= 2) {
    target_cond1 <- names(cond_counts)[1]
    target_cond2 <- names(cond_counts)[2]
    # Choose the most frequent time point for these conditions
    time_counts <- sort(table(metadata$time[metadata$condition %in% c(target_cond1, target_cond2)]), 
                        decreasing = TRUE)
    if (length(time_counts) >= 1) {
      target_time <- names(time_counts)[1]
    } else {
      target_time <- names(sort(table(metadata$time), decreasing = TRUE))[1]
    }
    print(paste("Manually selected:", target_cond2, "vs", target_cond1, "at", target_time))
  }
}

# Verify we have a valid comparison
if (!exists("target_cond1") || !exists("target_cond2") || !exists("target_time")) {
  print("ERROR: Could not determine a valid comparison. Please specify manually.")
  # As a fallback, just use the two most common conditions and time points
  cond_counts <- sort(table(metadata$condition), decreasing = TRUE)
  time_counts <- sort(table(metadata$time), decreasing = TRUE)
  if (length(cond_counts) >= 2 && length(time_counts) >= 1) {
    target_cond1 <- names(cond_counts)[1]
    target_cond2 <- names(cond_counts)[2]
    target_time <- names(time_counts)[1]
    print(paste("Fallback selection:", target_cond2, "vs", target_cond1, "at", target_time))
  } else {
    stop("Cannot proceed without valid comparison groups.")
  }
}

# Select samples for the chosen comparison
selected_samples <- which(
  (metadata$condition == target_cond1 | metadata$condition == target_cond2) &
    metadata$time == target_time
)

# Print selected samples for verification
print(paste("Selected", length(selected_samples), "samples for analysis:"))
print(metadata[selected_samples, ])

# Check if we have enough samples
if (length(selected_samples) < 2) {
  print("ERROR: Not enough samples for comparison. Will try any combination of conditions and time:")
  
  # Try a simpler approach - just compare any two conditions
  conditions <- unique(metadata$condition)
  if (length(conditions) >= 2) {
    target_cond1 <- conditions[1]
    target_cond2 <- conditions[2]
    selected_samples <- which(metadata$condition %in% c(target_cond1, target_cond2))
    print(paste("Comparing", target_cond2, "vs", target_cond1, "across all time points"))
    print(paste("Selected", length(selected_samples), "samples"))
    print(metadata[selected_samples, ])
  } else {
    stop("Cannot find any valid comparison. Please check your data.")
  }
}

# Subset counts and metadata
subset_counts <- filtered_counts[, selected_samples]
subset_metadata <- metadata[selected_samples, ]

# Convert factors for DESeq2
subset_metadata$condition <- factor(subset_metadata$condition)
if (length(unique(subset_metadata$time)) > 1) {
  subset_metadata$time <- factor(subset_metadata$time)
  # If we have both condition and time as factors, use a more complex design
  design_formula <- ~ condition + time
} else {
  # If only comparing conditions, use a simple design
  design_formula <- ~ condition
}

print("Design formula:")
print(design_formula)

# Create DDS object
dds <- DESeqDataSetFromMatrix(
  countData = subset_counts,
  colData = subset_metadata,
  design = design_formula
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
results_deseq <- results(dds, contrast = c("condition", target_cond2, target_cond1), alpha = 0.05)
summary(results_deseq)

# Filter significant DEGs (adjusted p-value < 0.05 and |log2FoldChange| > 1)
sigGenes <- subset(results_deseq, padj < 0.05 & abs(log2FoldChange) > 1)
print(paste("Number of significant DEGs:", nrow(sigGenes)))

# Display significant DEGs with HGNC symbols (no need for additional merge since rownames are HGNC)
sigGenes_df <- as.data.frame(sigGenes)
sigGenes_df$hgnc_symbol <- rownames(sigGenes_df)
print("Top differentially expressed genes (Control vs FSHD):")
print(head(sigGenes_df[, c("hgnc_symbol", "log2FoldChange", "padj")], 10))

# Check for DUX4
if ("DUX4" %in% sigGenes_df$hgnc_symbol) {
  print("DUX4 found in significant DEGs:")
  print(sigGenes_df[sigGenes_df$hgnc_symbol == "DUX4", c("hgnc_symbol", "log2FoldChange", "padj")])
} else {
  full_results_df <- as.data.frame(results_deseq)
  full_results_df$hgnc_symbol <- rownames(full_results_df)
  if ("DUX4" %in% full_results_df$hgnc_symbol) {
    print("DUX4 expression result (not significant but present):")
    print(full_results_df[full_results_df$hgnc_symbol == "DUX4", c("hgnc_symbol", "log2FoldChange", "padj")])
  } else {
    print("DUX4 not detected in this dataset.")
  }
}

# Save results to file
write.csv(as.data.frame(results_deseq), file = paste0("DESeq2_results_", target_cond2, "_vs_", target_cond1, ".csv"))

# Visualization 1: PCA plot to visualize sample clustering
vsd <- vst(dds, blind = FALSE)
if ("time" %in% colnames(subset_metadata) && length(unique(subset_metadata$time)) > 1) {
  pca_plot <- plotPCA(vsd, intgroup = c("condition", "time"))
} else {
  pca_plot <- plotPCA(vsd, intgroup = "condition")
}
print(pca_plot + theme_minimal() + ggtitle("PCA Plot of Samples"))

# Visualization 2: Volcano plot
if (nrow(results_deseq) > 0) {
  volcano_plot <- EnhancedVolcano(results_deseq,
                                  lab = rownames(results_deseq),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 3.0,
                                  labSize = 4.0,
                                  title = paste(target_cond2, 'vs', target_cond1),
                                  subtitle = 'Differential Expression')
  print(volcano_plot)
}

# Visualization 3: Heatmap of top DEGs
if (nrow(sigGenes) > 0) {
  # Get top 50 DEGs (or all if less than 50)
  top_genes <- rownames(sigGenes)[order(sigGenes$padj)[1:min(50, nrow(sigGenes))]]
  
  # Extract normalized counts for these genes
  normalized_counts <- counts(dds, normalized = TRUE)
  heatmap_data <- normalized_counts[top_genes, ]
  
  # Scale the data for better visualization
  heatmap_data_scaled <- t(scale(t(heatmap_data)))
  
  # Create annotation for the heatmap
  annotation_col <- data.frame(
    Condition = subset_metadata$condition,
    row.names = rownames(subset_metadata)
  )
  if ("time" %in% colnames(subset_metadata) && length(unique(subset_metadata$time)) > 1) {
    annotation_col$Time <- subset_metadata$time
  }
  
  # Create color palette
  ann_colors <- list(
    Condition = c(Control = "blue", FSHD1 = "green", FSHD2 = "red")
  )
  
  # Generate heatmap
  pheatmap(heatmap_data_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           fontsize_row = 8,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = paste("Top DEGs -", target_cond2, "vs", target_cond1))
}

# Alternative analysis using edgeR
# Create DGEList object
dge <- DGEList(counts = subset_counts, group = subset_metadata$condition)

# Calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Create design matrix for edgeR
design_edgeR <- model.matrix(design_formula, data = subset_metadata)

# Estimate dispersion
dge <- estimateDisp(dge, design_edgeR)

# Fit model
fit <- glmQLFit(dge, design_edgeR)

# Test for differential expression
if (ncol(design_edgeR) > 1) {
  # If there are multiple coefficients, test the condition effect
  # Find which column corresponds to the condition of interest
  coef_idx <- grep(paste0("condition", target_cond2), colnames(design_edgeR))
  if (length(coef_idx) == 0) coef_idx <- 2  # Default to second coefficient if not found
  qlf <- glmQLFTest(fit, coef = coef_idx)
} else {
  # Simple design, just test the intercept
  qlf <- glmQLFTest(fit, coef = 2)
}

results_edgeR <- topTags(qlf, n = Inf)

# Get significant genes
sig_edgeR <- results_edgeR$table[results_edgeR$table$FDR < 0.05 & abs(results_edgeR$table$logFC) > 1, ]
print(paste("edgeR: Number of significant DEGs:", nrow(sig_edgeR)))

# Save edgeR results
write.csv(as.data.frame(results_edgeR), file = paste0("edgeR_results_", target_cond2, "_vs_", target_cond1, ".csv"))

# Compare results between DESeq2 and edgeR
deseq2_genes <- rownames(sigGenes)
edger_genes <- rownames(sig_edgeR)
common_genes <- intersect(deseq2_genes, edger_genes)
library(VennDiagram)
# Venn diagram of overlap
venn_data <- list(DESeq2 = deseq2_genes, edgeR = edger_genes)

# Draw the Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("DESeq2", "edgeR"),
  filename = NULL,  # Do not save to file, just plot
  output = TRUE,
  main = "Overlap of DEGs between DESeq2 and edgeR",
  main.cex = 1.5,  # Adjust the title size
  cat.cex = 1.2,   # Adjust the category label size
  cat.pos = 0,     # Position of the category labels
  fill = c("blue", "red")  # Colors for the sets
)

# Display the plot
grid.draw(venn.plot)
# Session information for reproducibility
print(sessionInfo())


#Gene Annotation: If your gene IDs are Ensembl IDs, add annotation using biomaRt to map to symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = rownames(sigGenes),
                      mart = ensembl)
sigGenes_df <- merge(sigGenes_df, gene_symbols, by.x = "gene_id", by.y = "ensembl_gene_id")




#########################################################search for new genes that might be crucial for a phenotype#########################################################
###########################################################################################################################################################################
#Letting users to search for new genes that might be crucial for a phenotype, like FSHD, is a valuable feature, especially in exploratory RNA-seq analysis where novel discoveries are key. 


# Load required libraries
library(DESeq2)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(gplots)
library(RColorBrewer)
library(biomaRt)  # For HGNC conversion

# Set working directory - modify this to your path
setwd("/home/kaplanlab/Downloads")

# Extract tar file if not already extracted
if (!dir.exists("extracted_files")) {
  untar("GSE174301_RAW.tar", exdir = "extracted_files")
}

# List all RSEM gene results files
file_list <- list.files("extracted_files", pattern = "*rsem.genes.results.gz", full.names = TRUE)
print(paste("Found", length(file_list), "files"))

# Extract sample names from file paths for better labeling
sample_names <- gsub(".*/(.*?)_rsem\\.genes\\.results\\.gz", "\\1", file_list)
print("All sample names:")
print(sample_names)

# Create an empty matrix to store expected counts
gene_ids <- read.table(gzfile(file_list[1]), header = TRUE, row.names = 1, sep = "\t")
merged_counts <- matrix(0, nrow = nrow(gene_ids), ncol = length(file_list))
rownames(merged_counts) <- rownames(gene_ids)
colnames(merged_counts) <- sample_names

# Loop through each file and extract the expected counts
for (i in 1:length(file_list)) {
  data <- read.table(gzfile(file_list[i]), header = TRUE, row.names = 1, sep = "\t")
  merged_counts[, i] <- round(data$expected_count)
}

# Early HGNC Conversion
ensembl_ids <- sub("\\..*$", "", rownames(merged_counts))  # Remove version numbers
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = unique(ensembl_ids),
                      mart = ensembl)

# Create a mapping data frame
gene_mapping <- data.frame(
  ensembl_gene_id = ensembl_ids,
  hgnc_symbol = gene_symbols$hgnc_symbol[match(ensembl_ids, gene_symbols$ensembl_gene_id)],
  stringsAsFactors = FALSE
)
gene_mapping$hgnc_symbol[is.na(gene_mapping$hgnc_symbol) | gene_mapping$hgnc_symbol == ""] <- 
  gene_mapping$ensembl_gene_id[is.na(gene_mapping$hgnc_symbol) | gene_mapping$hgnc_symbol == ""]

# Update rownames of merged_counts with HGNC symbols
rownames(merged_counts) <- gene_mapping$hgnc_symbol
print("Gene IDs converted to HGNC symbols where available:")
print(head(gene_mapping, 10))

# Store the mapping for later use
gene_mapping_full <- gene_mapping  # Keep full mapping for gene search

# Filter low-expression genes
threshold <- ncol(merged_counts) * 0.5
filtered_counts <- merged_counts[rowSums(merged_counts > 0) > threshold, ]
print(paste("Genes before filtering:", nrow(merged_counts)))
print(paste("Genes after filtering low expression:", nrow(filtered_counts)))

# [Metadata creation and sample selection logic remains unchanged]
# Inspect sample names
print("Sample name patterns analysis:")
sample_patterns <- list()
for (name in sample_names) {
  if (grepl("Control", name)) sample_patterns$Control <- c(sample_patterns$Control, name)
  if (grepl("FSHD1", name)) sample_patterns$FSHD1 <- c(sample_patterns$FSHD1, name)
  if (grepl("FSHD2", name)) sample_patterns$FSHD2 <- c(sample_patterns$FSHD2, name)
  for (day in c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20)) {
    day_pattern <- paste0("Day", day)
    if (grepl(day_pattern, name)) {
      if (is.null(sample_patterns[[day_pattern]])) sample_patterns[[day_pattern]] <- character(0)
      sample_patterns[[day_pattern]] <- c(sample_patterns[[day_pattern]], name)
    }
  }
}

print("Detected patterns in sample names:")
for (pattern in names(sample_patterns)) {
  print(paste(pattern, ":", length(sample_patterns[[pattern]]), "samples"))
}

# Create metadata function
create_metadata <- function(sample_names) {
  condition <- character(length(sample_names))
  time_point <- character(length(sample_names))
  
  for (i in 1:length(sample_names)) {
    name <- sample_names[i]
    
    # Extract condition
    if (grepl("Control", name)) {
      condition[i] <- "Control"
    } else if (grepl("FSHD1", name)) {
      condition[i] <- "FSHD1"
    } else if (grepl("FSHD2", name)) {
      condition[i] <- "FSHD2"
    } else {
      condition[i] <- "Unknown"
    }
    
    # Extract time point with improved regex
    day_match <- regexpr("Day_\\d+", name)  # Match "Day_" followed by digits
    if (day_match > 0) {
      day_text <- substr(name, day_match, day_match + attr(day_match, "match.length") - 1)
      # Convert "Day_X" to "DayX" for consistency
      time_point[i] <- gsub("Day_", "Day", day_text)
    } else {
      time_point[i] <- "Unknown"
    }
  }
  
  metadata <- data.frame(
    sample = sample_names,
    condition = condition,
    time = time_point,
    stringsAsFactors = FALSE,
    row.names = sample_names
  )
  return(metadata)
}

# Generate metadata
metadata <- create_metadata(sample_names)
print("Generated metadata (first 20 rows):")
print(head(metadata, 20))

# Check unique values
print("Unique conditions in metadata:")
print(unique(metadata$condition))
print("Unique time points in metadata:")
print(unique(metadata$time))

# Dynamically select a comparison
condition_time_counts <- table(metadata$condition, metadata$time)
print("Detailed condition-time combinations:")
print(condition_time_counts)

# Find a valid comparison (at least 2 samples per condition at a time point)
valid_comparisons <- list()
for (time in colnames(condition_time_counts)) {
  for (cond1 in rownames(condition_time_counts)) {
    for (cond2 in rownames(condition_time_counts)) {
      if (cond1 != cond2 && condition_time_counts[cond1, time] >= 2 && condition_time_counts[cond2, time] >= 2) {
        valid_comparisons[[length(valid_comparisons) + 1]] <- list(cond1 = cond1, cond2 = cond2, time = time)
      }
    }
  }
}

if (length(valid_comparisons) > 0) {
  target_cond1 <- valid_comparisons[[1]]$cond1
  target_cond2 <- valid_comparisons[[1]]$cond2
  target_time <- valid_comparisons[[1]]$time
  print(paste("Selected comparison:", target_cond2, "vs", target_cond1, "at", target_time))
} else {
  top_conditions <- names(sort(table(metadata$condition), decreasing = TRUE)[1:2])
  target_cond1 <- top_conditions[1]
  target_cond2 <- top_conditions[2]
  target_time <- NULL
  print(paste("Fallback comparison:", target_cond2, "vs", target_cond1, "across all time points"))
}

# Select samples
if (is.null(target_time)) {
  selected_samples <- which(metadata$condition %in% c(target_cond1, target_cond2))
} else {
  selected_samples <- which((metadata$condition == target_cond1 | metadata$condition == target_cond2) & metadata$time == target_time)
}

print(paste("Selected", length(selected_samples), "samples for analysis:"))
if (length(selected_samples) < 2) {
  stop("ERROR: Insufficient samples selected. Check metadata and sample names.")
}
print(metadata[selected_samples, ])

# Subset counts and metadata
subset_counts <- filtered_counts[, selected_samples]
subset_metadata <- metadata[selected_samples, ]
subset_metadata$condition <- factor(subset_metadata$condition)

design_formula <- ~ condition
print("Design formula:")
print(design_formula)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = subset_metadata, design = design_formula)
dds <- DESeq(dds)
results_deseq <- results(dds, contrast = c("condition", target_cond2, target_cond1), alpha = 0.05)
summary(results_deseq)
# Filter significant DEGs
sigGenes <- subset(results_deseq, padj < 0.05 & abs(log2FoldChange) > 1)
print(paste("Number of significant DEGs:", nrow(sigGenes)))

# Display top DEGs
sigGenes_df <- as.data.frame(sigGenes)
sigGenes_df$hgnc_symbol <- rownames(sigGenes_df)
print("Top differentially expressed genes (Control vs FSHD):")
print(head(sigGenes_df[, c("hgnc_symbol", "log2FoldChange", "padj")], 10))


#WNT16 (log2FC = 2.94, padj = 2.94e-03): Upregulated in FSHD2. Part of the Wnt signaling pathway, which is involved in muscle development and regeneration.
#TMEM176A (log2FC = 6.15, padj = 1.23e-03): Strongly upregulated. A transmembrane protein, potentially novel in FSHD context—worth exploring further.
#RBM6 (log2FC = 1.63, padj = 2.34e-02): Upregulated. RNA-binding protein, possibly linked to RNA processing changes in FSHD.
#CAMKK1 (log2FC = 1.94, padj = 1.42e-02): Upregulated. Calcium/calmodulin-dependent kinase, involved in signaling pathways.
#HSPB6 (log2FC = 1.87, padj = 1.38e-04): Upregulated. Heat shock protein, potentially protective against muscle stress.
#PDK4 (log2FC = 1.98, padj = 3.54e-08): Upregulated with high significance. Pyruvate dehydrogenase kinase, linked to metabolic shifts in muscle.
#REXO5 (log2FC = -1.62, padj = 2.88e-05): Downregulated. An exonuclease, possibly affecting RNA stability.
#COPZ2 (log2FC = 1.54, padj = 2.49e-03): Upregulated. Coatomer protein, involved in intracellular transport.
#SOX8 (log2FC = -3.87, padj = 1.93e-13): Strongly downregulated with very high significance. A transcription factor, potentially a key regulator in FSHD.
#MAP3K9 (log2FC = -3.27, padj = 5.32e-06): Downregulated. A kinase in MAPK signaling, possibly linked to muscle pathology.
#DUX4 Absence: DUX4, a key gene in FSHD pathogenesis, isn’t in the top 10. 
#It might not meet the significance thresholds or could be absent/filtered out earlier (e.g., low expression in GSE174301).

# Load library for volcano plot
library(EnhancedVolcano)

# Prepare full results for plotting
full_results_df <- as.data.frame(results_deseq)
full_results_df$hgnc_symbol <- rownames(full_results_df)

# Define genes to highlight (top 10 from sigGenes_df + DUX4)
genes_to_highlight <- c("WNT16", "TMEM176A", "RBM6", "CAMKK1", "HSPB6", "PDK4", 
                        "REXO5", "COPZ2", "SOX8", "MAP3K9", "DUX4")

# Check presence of genes
present_genes <- genes_to_highlight[genes_to_highlight %in% full_results_df$hgnc_symbol]
missing_genes <- setdiff(genes_to_highlight, present_genes)
if (length(missing_genes) > 0) {
  print("Genes not found in DESeq2 results:")
  print(missing_genes)
}

# Add labeling column
full_results_df$label <- ifelse(full_results_df$hgnc_symbol %in% present_genes, 
                                full_results_df$hgnc_symbol, "")

# Create volcano plot
volcano_plot <- EnhancedVolcano(
  full_results_df,
  lab = full_results_df$label,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 4.0,
  title = paste(target_cond2, "vs", target_cond1, "Differential Expression"),
  subtitle = "Top DEGs and DUX4 (if present)",
  caption = "Red: padj < 0.05 & |log2FC| > 1, Grey: Non-significant",
  colAlpha = 0.5,
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey50",
  max.overlaps = 20
)

# Print the plot
print(volcano_plot)

# Display results for highlighted genes
print("Differential expression results for highlighted genes:")
highlighted_results <- full_results_df[full_results_df$hgnc_symbol %in% present_genes, 
                                       c("hgnc_symbol", "log2FoldChange", "padj")]
print(highlighted_results[order(highlighted_results$padj), ])





#Code with Gene Search and DUX4 Check

# Display top DEGs 
print(head(sigGenes_df[, c("hgnc_symbol", "log2FoldChange", "padj")], 10))

# Check for DUX4 explicitly
if ("DUX4" %in% sigGenes_df$hgnc_symbol) {
  print("DUX4 found in significant DEGs:")
  print(sigGenes_df[sigGenes_df$hgnc_symbol == "DUX4", c("hgnc_symbol", "log2FoldChange", "padj")])
} else {
  full_results_df <- as.data.frame(results_deseq)
  full_results_df$hgnc_symbol <- rownames(full_results_df)
  if ("DUX4" %in% full_results_df$hgnc_symbol) {
    print("DUX4 expression result (not significant but present):")
    print(full_results_df[full_results_df$hgnc_symbol == "DUX4", c("hgnc_symbol", "log2FoldChange", "padj")])
  } else {
    print("DUX4 not detected in DESeq2 results.")
    # Check if DUX4 was in the original data
    if ("DUX4" %in% rownames(filtered_counts)) {
      print("DUX4 present in filtered counts but not differentially expressed.")
    } else if ("DUX4" %in% rownames(merged_counts)) {
      print("DUX4 present in raw counts but filtered out.")
    } else {
      print("DUX4 absent from raw counts in GSE174301.")
    }
  }
}

# Interactive Gene Search
cat("\nSearch for a gene of interest (e.g., DUX4, TMEM176A, or an Ensembl ID like ENSG00000258399):\n")
gene_search <- readline(prompt = "Enter gene name or Ensembl ID: ")

full_results_df <- as.data.frame(results_deseq)
full_results_df$hgnc_symbol <- rownames(full_results_df)

if (gene_search %in% full_results_df$hgnc_symbol) {
  gene_result <- full_results_df[full_results_df$hgnc_symbol == gene_search, 
                                 c("hgnc_symbol", "log2FoldChange", "padj")]
  print("Gene found in results:")
  print(gene_result)
} else if (gene_search %in% gene_mapping_full$ensembl_gene_id) {
  hgnc_match <- gene_mapping_full$hgnc_symbol[gene_mapping_full$ensembl_gene_id == gene_search]
  if (hgnc_match %in% full_results_df$hgnc_symbol) {
    gene_result <- full_results_df[full_results_df$hgnc_symbol == hgnc_match, 
                                   c("hgnc_symbol", "log2FoldChange", "padj")]
    print(paste("Gene found via Ensembl ID (mapped to", hgnc_match, "):"))
    print(gene_result)
  } else {
    print(paste("Gene", gene_search, "mapped to", hgnc_match, "but not in DESeq2 results (possibly filtered out)."))
  }
} else {
  print("Gene not found in dataset. Check spelling or ID format.")
}






# Interactive Gene Search
cat("\nSearch for a gene of interest (e.g., DUX4, HNF1A, or an Ensembl ID like ENSG00000258399):\n")
gene_search <- readline(prompt = "Enter gene name or Ensembl ID: ")

# Process the search input
full_results_df <- as.data.frame(results_deseq)
full_results_df$hgnc_symbol <- rownames(full_results_df)

# Check if input matches HGNC or Ensembl ID
if (gene_search %in% full_results_df$hgnc_symbol) {
  gene_result <- full_results_df[full_results_df$hgnc_symbol == gene_search, 
                                 c("hgnc_symbol", "log2FoldChange", "padj")]
  print("Gene found in results:")
  print(gene_result)
} else if (gene_search %in% gene_mapping_full$ensembl_gene_id) {
  hgnc_match <- gene_mapping_full$hgnc_symbol[gene_mapping_full$ensembl_gene_id == gene_search]
  if (hgnc_match %in% full_results_df$hgnc_symbol) {
    gene_result <- full_results_df[full_results_df$hgnc_symbol == hgnc_match, 
                                   c("hgnc_symbol", "log2FoldChange", "padj")]
    print(paste("Gene found via Ensembl ID (mapped to", hgnc_match, "):"))
    print(gene_result)
  } else {
    print(paste("Gene", gene_search, "mapped to", hgnc_match, "but not in DESeq2 results (possibly filtered out)."))
  }
} else {
  print("Gene not found in dataset. Check spelling or ID format.")
  # Offer to search Ensembl online (simulated here)
  cat("Would you like to search Ensembl for this gene? (y/n): ")
  search_ensembl <- readline(prompt = "")
  if (tolower(search_ensembl) == "y") {
    print("Simulating Ensembl search... In practice, visit https://www.ensembl.org and enter your gene.")
    # Here you could integrate a web search if API access were available
  }
}








