#!/usr/bin/env Rscript

# RDS-file summary and export script

# Launched from acme-map.
# Load tables and export into appropriate RDS format compatible with MGS objects.
library(optparse)   # command-line arguments
library(yaml)       # config file
source("acme-map4/acmeFunctions.R")

# Hardcoded paths
outputDir <- "04.rds"
AlgorithmName <- "acme-map4:votu_summary:2023-10-09"
samplelistFilename <- "samplelist"

# command-line arguments                                                       
option_list <- list(
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
        help = "Virus catalog name",
        type = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)
outputPrefix <- args$output_prefix
virusCatalogName <- args$catalog_name
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

#################################### Load and export tables into RDS objects ####################################

# Quit programme if the votu count matrix do not exist.
if (!file.exists(file.path("20.votu/votu_matrices/votu.counts.mat.tsv"))){              
    stop("ERROR: Could not load the vOTU count matrix ... check the path ->  20.votu/votu_matrices/votu.counts.mat.tsv ")
}

# Read sample list
say1("Reading sample list:", samplelistFilename)
sampleNames <- read.table(samplelistFilename, stringsAsFactors = FALSE, fill = TRUE)[, 1]

# Summarise the total number of vOTU reads (including dropped reads) across samples
votuReadsFiles <- paste0("20.votu/kma/", "votu.reads.", sampleNames, ".tsv.gz")
votuReadsTotal <- sapply(votuReadsFiles, function(x) {
    x <- system2("zcat", args = c(x, "|", "wc", "-l"), stdout = TRUE)
    x <- as.integer(strsplit(x, split = " ")[[1]][1])
    return(x)
})
names(votuReadsTotal) <- sampleNames

# Profiled study count matrix
vOTU.counts <- read.table(file = "20.votu/votu_matrices/votu.counts.mat.tsv", 
                        sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
vOTU.counts <- t(vOTU.counts)

# Map votu-names to CMIDs and replace with CMIDs
cmid2origin <- HVR$cmid2origin
cmids <- cmid2origin[match( colnames(vOTU.counts), cmid2origin$votu), ]$cmid
colnames(vOTU.counts) <- cmids

# the kmaprofiling preceding this script re-orders columns/rows
# Re-order them to be consistent when re-profiling
# 1. order columns hvOTU.[0-9]* alpha-numerically
# 2. present rows in order of sampleNames
vOTU.counts <- vOTU.counts[sampleNames, sort(colnames(vOTU.counts)), drop = FALSE]

# Get virus genome lengths
vOTU.effectiveLength <- HVR$genome_length[colnames(vOTU.counts)]
saveAs(vOTU.effectiveLength)

# Calculate virus abundance matrix (reads/bp)
# Cannot use acmeFunctions.R:length_normalized_composition because it does not
# report the depth / scaling_factor
vOTU.depth <- sweep(vOTU.counts, MARGIN = 2, vOTU.effectiveLength, FUN = "/")
scaling_factor <- rowSums(vOTU.depth) / 100
vOTU.comp <- 100 * proportions(vOTU.depth, margin = 1)
attr(vOTU.comp, "ScalingFactor") <- scaling_factor
saveAs(vOTU.comp)
saveAs(vOTU.counts)

# Load Accounting object from MGS profiling (if existing)
MGSaccounting_path <- paste0("04.rds/", outputPrefix, ".accounting.rds")
if (file.exists(MGSaccounting_path)) {
    accounting <- readRDS(MGSaccounting_path)
    if( ! all(rownames(accounting) == rownames(vOTU.counts)) ){
        say2("ERROR: Something's wrong with MGS accounting sample order and vOTU matrix order")
        stop("ERROR: Check that you are using the same samplelist!")
    }
    vOTU.accounting <- accounting
    vOTU_mapped_reads <- as.integer(rowSums(vOTU.counts))
    vOTU.accounting$vOTU_mapped <- votuReadsTotal
    vOTU.accounting$vOTU_retained <- vOTU_mapped_reads
    vOTU.accounting$vOTU_mapped_proportion <- vOTU_mapped_reads / 
            (vOTU.accounting$TOTAL - vOTU.accounting$trimmed - vOTU.accounting$host)
} else {
    vOTU_mapped_reads <- as.integer(rowSums(vOTU.counts))
    vOTU.accounting <- data.frame(vOTU_mapped = votuReadsTotal, 
                                  vOTU_retained = as.integer(rowSums(vOTU.counts)),  
                                  row.names = rownames(vOTU.counts))
    vOTU.accounting$vOTU_mapped_proportion <- NA_real_                      
}

# Calculate viral bacterial ratio if MGS.comp is available
MGScomp_path <- paste0("04.rds/", outputPrefix, ".MGS.comp.rds")
if (file.exists(MGScomp_path)) {
    MGS.comp <- readRDS(MGScomp_path)
    if (!is.null(attr(MGS.comp,"ScalingFactor"))){
        scaling_factor <- attributes(MGS.comp)$ScalingFactor
        
        # Convert MGS Proportions to depth (avg. reads/bp)
        MGS.depth <- sweep(MGS.comp, MARGIN = 1, scaling_factor, FUN = "*")
        if( ! all(rownames(MGS.depth) == rownames(vOTU.depth)) ){
            say2("ERROR: Sample order of MGS and vOTU abundances do not match.")
            stop("ERROR: Are you analysing a different set of samples?")
        }
        # If either MGS or vOTU abundance is zero, the viral-bacterial ratio is NaN 
        vOTU.accounting$viral_bacterial_ratio <- rowSums(vOTU.depth) / rowSums(MGS.depth)

    } else {
        say2("The <MGS.comp> object did not have ScalingFactor attribute. No <MGS.votu.comp.rds> Export.")
        say2("If unsure inspect the attributes of MGS.comp object here -> ", MGScomp_path)
        vOTU.accounting$viral_bacterial_ratio <- NA_real_
    }
} else {
    say2("MGS.comp object not found. Viral-bacterial ratio not calculated.")
    vOTU.accounting$viral_bacterial_ratio <- NA_real_
}
saveAs(vOTU.accounting)

# Also write a tsv version for machine parsing in the Rarefaction step
vOTU.accounting$sampleid <- rownames(vOTU.accounting)
write.table(vOTU.accounting, file = paste0("20.votu/", outputPrefix, ".vOTU.accounting.tsv"),
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
