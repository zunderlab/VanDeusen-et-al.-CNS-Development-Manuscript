#Amy Van Deusen, University of Virginia
#20 September 2021
## Downsample by metadata and cluster according to .csv containing values

print("Start 01_Downsample_By_Metadata_And_Cluster.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ZunderPipelineFunctions)
library(data.table)
library(dplyr)

print("Finished loading libraries, getting input parameters")

## Input parameters =====================================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "/metadata.csv"
CONCAT.TRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
LAYOUT.FILENAME <- "/UMAP_layout.csv"
CLUSTERS.FILENAME <- "/clusters.csv"
FILENUMS.FILENAME <- "filenums.csv" 
DOWNSAMPLE.VALUES.FILENAME <- "/downsample_values.csv" #.csv with downsampling values in order of metadata + cluster (e.g. E11 cluster 1, E11 cluster 2, E11 cluster 3...P4 cluster 1, P4 cluster 2, etc.)
METADATA.TYPE <- "Age" # Typically "Age" but could be any metadata type
DOWNSAMPLE.OUT.FILENAME <- "Concat_Transformed_Subsampled.csv" #new file
SUBSAMPLE.ID.OUT.FILENAME <- "Subsample_IDs.csv" #new file
LAYOUT.OUT.FILENAME <- "UMAP_layout_Subsampled.csv" #new file
CLUSTERS.OUT.FILENAME <- "clusters_Subsampled.csv" #new file
FILENUMS.OUT.FILENAME <- "filenums_Subsampled.csv" #new file

print("Finished getting input and output parameters, reading subsampling values")

## Read in necessary functions ====================================================
#' Metadata file input
#'
#' Reads metadata file and produces format friendly to scripts in Zunder lab pipeline. Metadata
#' strategy inspired by Nowicka et al., 2017, F1000 Research
#' @param input.folder directory containing metadata file
#' @param md.filename metadata filename, defaults to metadata.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.metadata <- function(input.folder,md.filename = "/metadata.csv"){
  md <- read.csv(paste0(input.folder,md.filename))
  #make filenames character vectors
  md$file_name <- as.character(md$file_name)
  return(md)
}

#' Layout file input
#'
#' Reads uMAP or other 2D layout for plotting
#' @param input.folder directory containing layout file
#' @param layout.filename layout filename, defaults to UMAP_layout.csv as outputted by generation
#'  script
#' @importFrom utils read.csv
#' @export
read.layout <- function(input.folder,layout.filename = "/UMAP_layout.csv"){
  layout <- read.csv(paste0(input.folder,layout.filename))
  return(layout)
}

#' Cluster file input
#'
#' Reads clusters from pipeline
#' @param input.folder directory containing layout file
#' @param clusters.filename clustering replicate filename, defaults to clusters.csv as outputted by
#'  generation script
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @export
read.clusters <- function(input.folder,clusters.filename = "/clusters.csv"){
  clusters.in <- fread(paste0(input.folder,clusters.filename))
  if (min(clusters.in)==0){
    #Add 1 since clustering done in Python begins indexing at zero
    clusters <- as.vector(t(clusters.in) + 1)
  } else{
    clusters <- as.vector(t(clusters.in))
  }
  
  return(clusters)
}

## Read needed files ====================================================================================
metadata <- read.metadata(INPUT.FOLDER,METADATA.FILENAME)
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCAT.TRANSFORMED.FILENAME), stringsAsFactors = F)
layout.in <- read.layout(INPUT.FOLDER,LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)
downsample.in <- read.csv(paste0(INPUT.FOLDER, DOWNSAMPLE.VALUES.FILENAME), header = TRUE)
filenums <- read.table(FILENUMS.FILENAME, header=TRUE, sep=",", check.names=FALSE)[,]

# Prepare sample metadata (specified in METADATA.TYPE) and clusters for addition to dataframe
metadata_key <- as.character(metadata[,METADATA.TYPE])
names(metadata_key) <- 1:length(metadata[,METADATA.TYPE])
metadata.in <- metadata_key[filenums]
metadata.in <- as.data.frame(metadata.in)
colnames(metadata.in) <- METADATA.TYPE
downsample.values <- as.vector(as.matrix(downsample.in))
METADATA.ORDER <- colnames(downsample.in)

print("Finished reading needed files, generating full dataframe")

## Generate list of subsample.ids based on metadata and cluster==========================================
# Initialize lists for outputs
numcluster <- length(unique(clusters.in))
subsample.ids <- vector(mode = "list", length = numcluster)
subsample.ids.out <- vector(mode = "list")

print("Finished initializing lists, downsampling")

# For each metadata item (e.g. Age, SampleType), downsample by each cluster according to downsample_values.csv
for (m in METADATA.ORDER) {
  for (i in 1:numcluster) {
    mcount <- match(m, METADATA.ORDER)
    iteration.count <- i + numcluster*(mcount - 1)
    DOWNSAMPLE.VALUE <- downsample.values[iteration.count]
    if (DOWNSAMPLE.VALUE >= 1) {
      all.ids <- which((clusters.in == i) & (metadata.in == m))
      if (length(all.ids) >= DOWNSAMPLE.VALUE) {
        subsample.ids[[i]] <- sample(all.ids,size = DOWNSAMPLE.VALUE, replace = FALSE)
        cat("\r", "Cluster", i, "of", m, "downsampled to", DOWNSAMPLE.VALUE, "cells.")
      } else if ((length(all.ids) < DOWNSAMPLE.VALUE) & (length(all.ids) > 0)) {
        subsample.ids[[i]] <- sample(all.ids, size = length(all.ids), replace = FALSE)
        cat("\r", "Cluster", i, "of", m, "only has", DOWNSAMPLE.VALUE, "cells, therefore all cells were included in subsampling.")
      } else {
        subsample.ids[[i]] <- NA
        cat("Zero cells downsampled from cluster", i, "of", m)
      }
    } else {
      subsample.ids[[i]] <- NA
      cat("\r", "Cluster", i, "of", m, "was excluded from subsampling.")  
    }
  }  
  subsample.ids.out[[m]] <- subsample.ids
}  

#print(subsample.ids.out)
print("Finished subsampling by metdata and cluster, generating final list of subsample.ids")

#Convert list of all subsample.ids into dataframe and save .csv output
subsample.ids.unlist <- unlist(subsample.ids.out)
subsample.ids.unlist <- na.omit(subsample.ids.unlist)
subsample.ids.final <- as.data.frame(subsample.ids.unlist)
colnames(subsample.ids.final) <- "subsample.id"
subsample.ids.final <- sort(subsample.ids.final$subsample.id)
subsample.id.out.filename <- paste0(OUTPUT.FOLDER,"/",SUBSAMPLE.ID.OUT.FILENAME)
fwrite(list(subsample.ids.final),file = subsample.id.out.filename, col.names = FALSE,sep = ",")

print("Final list of subsample.ids output, subsampling expression matrix, UMAP layout, and clusters")

## Subsample expression matrix, UMAP layout, and clusters according to subsample.ids ====================
# Initialize lists for outputs
expression.out <- subsample.ids.final
layout.out <- subsample.ids.final
clusters.out <- subsample.ids.final
filenums.out <- subsample.ids.final

# Subsample expression matrix and output new .csv
expression.out <- concat.transformed[subsample.ids.final]
downsample.out.filename <- paste0(OUTPUT.FOLDER,"/",DOWNSAMPLE.OUT.FILENAME)
fwrite(expression.out, file = downsample.out.filename, row.names = FALSE)

# Subsample UMAP layout and output new .csv
layout.out <- layout.in[subsample.ids.final,]
layout.out.filename <- paste0(OUTPUT.FOLDER,"/",LAYOUT.OUT.FILENAME)
fwrite(layout.out,file = layout.out.filename)

# Subsample clusters and output new .csv
clusters.out <- clusters.in[subsample.ids.final]
clusters.out.filename <- paste0(OUTPUT.FOLDER,"/",CLUSTERS.OUT.FILENAME)
fwrite(list(clusters.out), file = clusters.out.filename, col.names = FALSE,sep = ",")

# Subsample filenums and output new .csv
filenums.out <- filenums[subsample.ids.final]
filenums.df <- as.data.frame(filenums.out)
colnames(filenums.df) <- "File"
filenums.out.filename <- paste0(OUTPUT.FOLDER,"/",FILENUMS.OUT.FILENAME)
fwrite(filenums.df,file = filenums.out.filename)

print("Output new .csv files for expression matrix, UMAP layout, and clusters")
print("Finish 01_Downsample_By_Metadata_And_Cluster.R")

  
