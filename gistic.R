# Load libraries
library(DNAcopy)    # For segmentation (CBS method)
library(oligo)      # For reading and processing CEL files
library(pd.genomewidesnp.6)

# Set the directory containing CEL files
cel_files <- list.files(path = "/datasets/marta/TCGA/CEL_CNV_TCGA/COAD", pattern = "*/*.CEL", recursive=T, full.names = T)
cel_files <- subset(cel_files, !grepl("logs", cel_files))

# Read CEL files
snp_data <- read.celfiles(cel_files, pkgname = "pd.genomewidesnp.6")

# Normalize the data (RMA normalization is commonly used)
normalized_data <- oligo::normalize(snp_data)

# Extract expression values (log2 scale after normalization)
expr_values <- exprs(normalized_data)

# Convert to a matrix for further processing, if needed
expr_matrix <- as.matrix(expr_values)

# Obtain feature (probe) data including genomic positions
feature_data <- pData(featureData(normalized_data))
