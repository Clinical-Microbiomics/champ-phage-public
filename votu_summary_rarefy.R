#!/usr/bin/env Rscript

# vOTU rarefaction summary script 

# Launched from acme-map.
# Load tables and export into appropriate RDS format compatible with MGS objects.
library(optparse)   # command-line arguments
library(yaml)       # config file
source("acme-map4/acmeFunctions.R")

# Hardcoded paths
outputDir <- "04.rds"
AlgorithmName <- "acme-map4:votu_rarefaction:2023-10-09"
samplelistFilename <- "samplelist"

# command-line arguments                                                       
option_list <- list(
    make_option(
        c("-a", "--downsize_algorithm"),
        help = "The rarefaction algorithm",
        type = "character"
    ),
    make_option(
        c("-t", "--downsize_target"),
        help = "The rarefaction target",
        type = "character"
    ),
    make_option(
        c("-v", "--virusObjectFile"),
        help = "Path to the HVR catalog rds",
        type = "character"
    ),
    make_option(
        c("-o", "--output_prefix"),
        help = "Output prefix name of .RDS files",
        type = "character"
    ),
    make_option(
        c("-c", "--catalog_name"),
        help = "Gene catalog name",
        type = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)
outputPrefix <- args$output_prefix
virusCatalogName <- args$catalog_name
downsizeTarget <- args$downsize_target
downsizeAlgorithm <- args$downsize_algorithm
virusObjectFile <- args$virusObjectFile

# RDS export arguments and attributes
CM.Info <- c(CreationTime = format(Sys.time()),
             CreationDir  = getwd(),
             User         = unname(Sys.info()["user"]),
             geneCatalogName = virusCatalogName,
             Algorithm    = AlgorithmName,
             AcmeMapVersion = AcmeMapVersion) # AcmeMapVersion is defined in acmeFunctions.R

saveAs <- function(x) {
    name <- deparse(substitute(x))
    attr(x, "CM.Info") <- CM.Info
    fileName <- paste0(outputPrefix, ".", name, ".rds")
    comment(x) <- paste(fileName, format(Sys.time()), unname(Sys.info()["user"]))
    say2("Saving", fileName)
    saveRDS(x, file = file.path(outputDir, fileName))
}

#################################### Load MGV object ####################################
HVR <- readRDS(virusObjectFile)
say2("Version:", HVR$version)
say2("Number of vOTUs:", nrow(HVR$cmid2origin))
say2("Comment:", comment(HVR))


vOTUaccounting_path <- paste0("04.rds/", outputPrefix, ".vOTU.accounting.rds")
if (!file.exists(file.path(vOTUaccounting_path))){              
    stop("ERROR: Could not load the vOTU accounting object")
}
vOTU.accounting <- readRDS(vOTUaccounting_path)

# Rarefaction of vOTU counts
rarefied_votu_matrix <- "21.votu.rarefied/votu_matrices/votu.counts.mat.tsv"
if (!file.exists(rarefied_votu_matrix)){
    stop("ERROR: The 21.votu.rarefied/votu_matrices/votu.counts.mat.tsv could not be found!")
}
    
# Load count matrix
vOTU.dsCounts <- read.table(file = rarefied_votu_matrix, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
vOTU.dsCounts <- t(vOTU.dsCounts)

# Map votu-names to CMIDs and replace with CMIDs
cmid2origin <- HVR$cmid2origin
cmids <- cmid2origin[match(colnames(vOTU.dsCounts), cmid2origin$votu), ]$cmid
colnames(vOTU.dsCounts) <- cmids

# Read sample list
say1("Reading sample list:", samplelistFilename)
sampleNames <- read.table(samplelistFilename, stringsAsFactors = FALSE, fill = TRUE)[, 1]

# the kmaprofiling preceding this script re-orders columns/rows
# Re-order them to be consistent when re-profiling
# 1. order columns hvOTU.[0-9]* alpha-numerically
# 2. present rows in order of sampleNames
vOTU.dsCounts <- vOTU.dsCounts[sampleNames, sort(colnames(vOTU.dsCounts)), drop = FALSE]

# Get virus genome lengths
vOTU.effectiveLength <- HVR$genome_length[colnames(vOTU.dsCounts)]

# Calculate virus abundance matrix (reads/bp)
# Cannot use acmeFunctions.R:length_normalized_composition because it does not
# report the depth / scaling_factor
vOTU.dsdepth <- sweep(vOTU.dsCounts, MARGIN = 2, vOTU.effectiveLength, FUN = "/")
scaling_factor <- rowSums(vOTU.dsdepth) / 100
vOTU.dsComp <- 100 * proportions(vOTU.dsdepth, margin = 1)
vOTU.downsize.target <- getDownsizeTarget(
    downsizeTarget,
    downsizeAlgorithm,
    vOTU.accounting$vOTU_mapped
)
say2("vOTU downsize target:", vOTU.downsize.target)

# Set rows to NA if sample counts are below target
belowTarget <- vOTU.accounting$vOTU_mapped < vOTU.downsize.target
if( ! all(rownames(vOTU.accounting) == rownames(vOTU.dsCounts)) ){
    say2("ERROR: Something's wrong with accounting sample order and vOTU matrix order")
    stop("ERROR: Check that you are using the original samplelist!")
}
vOTU.dsCounts[belowTarget, ] <- NA_integer_
vOTU.dsComp[belowTarget, ] <- NA_real_
attr(vOTU.dsComp, "ScalingFactor") <- scaling_factor
attr(vOTU.dsComp, "downsizeTarget") <- vOTU.downsize.target
attr(vOTU.dsCounts, "downsizeTarget") <- vOTU.downsize.target
saveAs(vOTU.dsComp)
saveAs(vOTU.dsCounts)

# Update vOTU accouting table with results post rarefaction
vOTU_rarefied_mapped <- as.integer(rowSums(vOTU.dsCounts))
vOTU_rarefied_richness <- as.integer(rowSums(vOTU.dsComp > 0))
vOTU.accounting$dsvOTU_retained <- vOTU_rarefied_mapped
vOTU.accounting$dsvOTU_richness <- vOTU_rarefied_richness
saveAs(vOTU.accounting)
