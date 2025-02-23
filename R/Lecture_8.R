#Lecture_8
library(data.table)
library(locuszoomr)
library(AnnotationHub)
library(ensembldb)
library(biomaRt)

# Load the GWAS Catalog data
#https://www.ebi.ac.uk/gwas/publications/34210852
#Sequencing of 640,000 exomes identifies GPR75 variants associated with protection from obesity.
gwas_data <- read.table("/Users/sebihacevik/Downloads/GPR75_2025-02-23-rsId_rs114285050.tsv", header = TRUE, sep = "\t")
# Print the first few rows to check the data
head(gwas_data)
# Print the column names to understand the data structure
names(gwas_data)

#Extract rsID: You'll still need to extract the rsID from the STRONGEST.SNP.RISK.ALLELE column.
gwas_data$rsid <- sub("-.$", "", gwas_data$STRONGEST.SNP.RISK.ALLELE)
names(gwas_data)

#Rename columns (optional): You can rename the relevant columns to match locuszoomr expectations
colnames(gwas_data)[colnames(gwas_data) == "P.VALUE"] <- "p"
colnames(gwas_data)[colnames(gwas_data) == "OR.or.BETA"] <- "beta"
colnames(gwas_data)[colnames(gwas_data) == "X95..CI..TEXT."] <- "CI"
colnames(gwas_data)[colnames(gwas_data) == "MAPPED_GENE"] <- "mapped_gene"
names(gwas_data)

ah <- AnnotationHub()
ensDb_query <- query(ah, c("EnsDb", "Homo sapiens", "v109"))

## Check the available resources
ensDb_query

# Select the correct resource (usually the first one)
ensDb_resource <- ensDb_query[[1]]  

# Extract the file path of the downloaded database
#cache_path <- cache(ah, "AH109606")  # Change "AH109606" if needed
#ensDb <- ensembldb::EnsDb(cache_path)
## Load the database using the file path
#ensDb <- ensembldb::EnsDb(ensDb_file)

# Fetch only the required EnsDb database
ensDb <- AnnotationHub()[["AH109606"]]

keytypes(ensDb)
columns(ensDb)

# Extract unique gene symbols from your dataset
gene_symbols <- unique(gwas_data$mapped_gene)

# Query EnsDb using gene symbols
variant_info <- ensembldb::select(ensDb, 
                                  keys = gene_symbols, 
                                  keytype = "SYMBOL", 
                                  columns = c("GENEID", "SYMBOL", "SEQNAME", "GENESEQSTART", "GENESEQEND", "TXNAME"))
# View results
head(variant_info)
# Check the first few results
head(variant_info)
# Connect to Ensembl
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# Extract the rsID column from your data
rsids <- gwas_data$rsid

# Retrieve SNP positions
snp_info <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id", "associated_gene"),
                  filters = "snp_filter",
                  values = rsids,
                  mart = ensembl)

# View results
head(snp_info)

#Merge with Your GWAS Data
# Merge GWAS data with SNP positions
gwas_data_merged <- merge(gwas_data, snp_info, by.x = "rsid", by.y = "refsnp_id", all.x = TRUE)

# View the first few rows
head(gwas_data_merged)

#Prepare Data for Visualization (e.g., LocusZoom)
#If you plan to visualize the SNPs with locuszoomr, format the data

# Select relevant columns for plotting
gwas_plot_data <- gwas_data_merged[, c("rsid", "p", "chr_name", "chrom_start")]

# Rename columns for compatibility with locuszoomr

# Ensure data format is correct
gwas_plot_data <- gwas_data_merged[, c("rsid", "p", "chr_name", "chrom_start")]
colnames(gwas_plot_data) <- c("SNP", "P", "CHR", "BP")

# Convert chromosome number to character if needed
gwas_plot_data$CHR <- as.character(gwas_plot_data$CHR)

library(locuszoomr)
#List all available functions from locuszoomr
ls("package:locuszoomr")

# Define the locus argument explicitly
# Create a locus object
# Create a locus object using a named list
loc <- list(chr = unique(gwas_plot_data$CHR)[1], 
            bp = c(min(gwas_plot_data$BP, na.rm = TRUE) - 50000,  
                   max(gwas_plot_data$BP, na.rm = TRUE) + 50000))


library(ggplot2)
library(dplyr)

# Enhanced regional plot function with lead SNP highlighting and gene annotation
# First, let's add a data preparation function
prepare_gwas_data <- function(gwas_data) {
  # Remove duplicates while keeping the most significant p-value for each SNP
  gwas_data <- gwas_data %>%
    group_by(SNP) %>%
    slice(which.min(P)) %>%
    ungroup()
  
  # Sort by position
  gwas_data <- gwas_data[order(gwas_data$BP), ]
  
  return(gwas_data)
}

# Before plotting, deduplicate SNPs keeping the most significant entry
gwas_plot_data <- gwas_plot_data %>%
  group_by(SNP) %>%
  arrange(P) %>%  # Keep lowest p-value per SNP
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::filter(!is.na(CHR) & !is.na(BP))  # Remove NA positions

# Verify uniqueness
cat("Unique SNPs after deduplication:", length(unique(gwas_plot_data$SNP)), "\n")

# Get all genes in the locus region
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Define region boundaries
chromosome <- unique(gwas_plot_data$CHR)[1]
start_pos <- min(gwas_plot_data$BP) - 500000  # 500kb window
end_pos <- max(gwas_plot_data$BP) + 500000

# Get genes in region
gene_annotation <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = c("chromosome_name", "start", "end"),
  values = list(chromosome, start_pos, end_pos),
  mart = ensembl
) %>%
  dplyr::filter(hgnc_symbol != "") %>%  # Remove unnamed genes
  distinct()

# Enhanced regional plot function with improved SNP handling
create_enhanced_regional_plot <- function(gwas_data, gene_data, 
                                          significance_threshold = 5e-8,
                                          lead_snp = NULL,
                                          label_threshold = 1e-5) {
  
  # Prepare data
  lead_snp <- if(is.null(lead_snp)) gwas_data$SNP[which.min(gwas_data$P)] else lead_snp
  gwas_data <- gwas_data %>%
    mutate(
      highlight = (SNP == lead_snp),
      logP = -log10(P),
      significant = P < significance_threshold
    )
  
  # Create association plot
  p_assoc <- ggplot(gwas_data, aes(x = BP, y = logP)) +
    geom_point(aes(color = significant), alpha = 0.7) +
    geom_point(data = subset(gwas_data, highlight), 
               color = "red", size = 3, shape = 23, fill = "yellow") +
    geom_hline(yintercept = -log10(significance_threshold), 
               linetype = "dashed", color = "darkred") +
    geom_text_repel(
      data = subset(gwas_data, P < label_threshold | highlight),
      aes(label = SNP),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c("grey", "blue")) +
    theme_bw() +
    labs(x = NULL, y = "-log10(p-value)", title = lead_snp)
  
  # Create gene track
  p_genes <- ggplot(gene_data) +
    geom_segment(aes(x = start_position, xend = end_position, y = 1, yend = 1),
                 color = "darkblue", size = 5) +
    geom_text_repel(aes(x = (start_position + end_position)/2, y = 1, label = hgnc_symbol),
                    angle = 90, hjust = 0.5, direction = "y") +
    scale_y_continuous(limits = c(0.9, 1.1)) +
    theme_void() +
    theme(axis.text.x = element_text())
  
  # Combine plots
  p_assoc / p_genes + 
    plot_layout(heights = c(3, 1)) +
    scale_x_continuous(labels = scales::comma_format(scale = 1e-6, suffix = "Mb"))
}

# Load required package for text labels
library(ggrepel)

# Example usage:
p <- create_enhanced_regional_plot(
  gwas_data = gwas_plot_data,
  gene_data = variant_info,
  significance_threshold = 5e-8,
  lead_snp = "rs114285050",
  label_threshold = 1e-6  # Adjust this threshold to show more/fewer SNP labels
)

final_plot <- create_enhanced_regional_plot(
  gwas_data = gwas_plot_data,
  gene_data = gene_annotation,
  lead_snp = "rs114285050",
  label_threshold = 1e-3
)
final_plot

# Save the plot with appropriate dimensions for the gene track
#ggsave("final_regional_plot.png", final_plot, width = 10, height = 8, dpi = 300)


