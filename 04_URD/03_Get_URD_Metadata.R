#Altered from GetURDMetadata.R written by Corey Williams, University of Virginia
#2021
#Adds optional columns to URD metdata file in addition to age/stage (e.g. Tissue and SampleType)

print("Start Get_URD_Metadata.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(data.table)
library(dplyr)

print("libraries loaded")

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
STAGE.TAG <- "stage" #column name in metadata for developmental stage
METADATA.FILENAME <- "/metadata.csv" #Must include "stage" column for URD
FILENUMS.FILENAME <- "filenums.csv" 
SUBSAMPLE.IDS.FILENAME <- "Subsample_IDs.csv"
CONCATTRANSFORMED.FILENAME <- "/Concat_Transformed_Subsampled.csv" #.csv file for pipeline clustering
CLUSTERS.FILENAME <- "/clusters_Subsampled.csv"
URD.META.FILENAME <- "metadata_URD.csv" #.txt file For URD - unused if no subsampling
#METADATA.TYPE <- "SampleType" #Optional: can add additional metadata for plotting

print("got input parameters, loading input files")

## Read in necessary functions =======================================================
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

#' Concat_Transformed file input
#'
#' Reads Concat_Transformed file for use in Zunder lab pipeline
#' @param input.folder directory containing concat_transformed file
#' @param concat.transformed.filename concat_transformed filename, defaults to
#' Concat_Transformed.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.concat.transformed <- function(input.folder,
                                    concat.transformed.filename = "/Concat_Tranformed.csv"){
  concat.transformed <- read.csv(paste0(input.folder,concat.transformed.filename))
  colnames(concat.transformed) <- gsub("-", "_", colnames(concat.transformed))
  colnames(concat.transformed) <- gsub(".", "_", colnames(concat.transformed), fixed = TRUE)
  return(concat.transformed)
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

## Read input files ===============================================================================
metadata <- read.metadata(INPUT.FOLDER,METADATA.FILENAME)
concat.transformed <- read.concat.transformed(INPUT.FOLDER,CONCATTRANSFORMED.FILENAME)
filenums <- read.table(FILENUMS.FILENAME, header=TRUE, sep=",", check.names=FALSE)[,]
subsample.ids <- as.numeric(unlist(fread(SUBSAMPLE.IDS.FILENAME, stringsAsFactors=FALSE)))
cluster.ids <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)

print("read input files")

## Construct "file_name" column for urd.md (to match original)=====================================
filenames_key <- as.character(metadata[,"file_name"]) # Get list of file names
names(filenames_key) <- 1:length(metadata[,"file_name"]) # Add index to list of stage IDs
filenames.in <- filenames_key[filenums] # Obtain file names for all cells
filenames.out <- filenames.in[subsample.ids]  # Subset file numbers for subsampled cells
filenames.out <- as.data.frame(filenames.out)
colnames(filenames.out) <- "file_name"

print("file names generated for urd.md")

## Construct "stage" column for urd.md ============================================================
stages_key <- as.character(metadata[,STAGE.TAG]) # Get list of stage IDs
names(stages_key) <- 1:length(metadata[,STAGE.TAG]) # Add index to list of stage IDs
stages.in <- stages_key[filenums] # Obtain stage.ids for all cells
stages.out <- stages.in[subsample.ids] # Subset stage.ids for subsampled cells
stages.out <- as.data.frame(stages.out) 
colnames(stages.out) <- STAGE.TAG

print("stage.ids generated for urd.md")

## Get list of file numbers for roots and tips ==============================================
filenums.out <- filenums[subsample.ids]  # Subset file numbers for subsampled cells
filenums.out <- as.data.frame(filenums.out)
colnames(filenums.out) <- "file_number"

print("file numbers imported")

## Construct "Roots" column for urd.md ============================================================
root.file.ids <- which(metadata$Root_clusters != "") # Get ids for files containing root cells
root.cluster.ids <- as.numeric(na.omit(unique(metadata$Root_clusters))) # Get root cluster ids
# NOTE: Can use multiple clusters BUT WILL BE A SINGLE ROOT IN URD ANALYSIS!
roots <- vector(mode = "logical",length = length(subsample.ids)) #Inititalize logical vector
for (i in 1:length(subsample.ids)) {
  if ((cluster.ids[i] %in% root.cluster.ids) & (filenums.out[i,] %in% root.file.ids))  {
    roots[[i]] = TRUE
  } else {
    roots[[i]] = FALSE
  }
}
roots.out <- as.data.frame(roots)
colnames(roots.out) <- "roots"

print("root cells identified for urd.md")

## Construct "Tips" column for urd.md =============================================================
tip.file.ids <- which(metadata$Tip_clusters != "") # Get ids for files containing root cells
tip.cluster.ids <- unique(metadata$Tip_clusters) # Get list of ids for all tip clusters
tip.cluster.ids <- tip.cluster.ids[tip.cluster.ids != ""] # Remove empty elements from list
tip.cluster.ids <- unlist(strsplit(tip.cluster.ids, ",")) # Split cluster ids into a vector
tips <- vector() # Initialize list for tip.ids
`%!in%` <- Negate(`%in%`) # Create not-in function to exclude 

for (i in 1:length(subsample.ids)) {
  if (filenums.out[i,] %!in% tip.file.ids) {
    tips[[i]] = 0
  } else if (cluster.ids[i] %!in% tip.cluster.ids) {
    tips[[i]] = 0
  } else if ((filenums.out[i,] %in% tip.file.ids) & (cluster.ids[i] %in% tip.cluster.ids)) {
    tips[[i]] = match(cluster.ids[i], tip.cluster.ids)
  } 
  else {
    cat("Error obtaining tip for sample", paste0(i))
  }
}
print("The following tips were generated for URD analysis:")
print(sort(unique(tips))) #Sanity check that all tips were sampled
tips.out <- as.data.frame(tips)
colnames(tips.out) <- "tips"

print("tips cells identified for urd.md")

## Generate urd.md output =========================================================================
urd.md.in <- c(filenames.out, stages.out, roots.out, tips.out) # Make list of dataframes to combine
urd.md <- do.call(cbind, urd.md.in) # Combine dataframs into single data frame
fwrite(urd.md, file = URD.META.FILENAME)

print("Finish Get_URD_Metadata.R")








