#Amy Van Deusen, University of Virginia
#1 July, 2021
#Convert outputs of Leiden clustering for URD analysis
#1) Rename expression_matrix.csv to Concat_Transformed.csv
#2) Rename umap_xy_R2.csv to umap_layout.csv
#3) Convert R2_assigns.csv to clusters.csv

print("Start 01_URD_Prep.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(data.table)
library(dplyr)

## Input parameters ===============================================================================

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
EXPRESSION.MATRIX.FILENAME <- "/expression_matrix_analysis.csv" #Standard name for Pipeline2.0
FILE.NUMBERS.FILENAME <- "/filenums.csv" #Standard name for Pipeline2.0
UMAP.LAYOUT.FILENAME <- "/umap_xy_R1.csv"  #Need to change for round of clustering
CLUSTERS.FILENAME <- "/cluster_R1_assigns.csv" #Need to change for round of clustering

## Output parameters ==============================================================================
CONCATTRANSFORMED.OUT.FILENAME <- "./Concat_Transformed.csv" #Standard name for Pipeline1.0
UMAP.LAYOUT.OUT.FILENAME <- "./UMAP_layout.csv" #Standard name for Pipeline1.0
CLUSTERS.OUT.FILENAME <- "./clusters.csv" #Standard name for Pipeline1.0

print("Finished reading input and output parameters, reading files")

## Read in necessary functions ======================================================
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



## Read needed files ==============================================================================
expression.matrix <- fread(paste0(INPUT.FOLDER, EXPRESSION.MATRIX.FILENAME), stringsAsFactors = F)
file.nums <- fread(paste0(INPUT.FOLDER, FILE.NUMBERS.FILENAME), stringsAsFactors = F)
layout.in <- read.layout(INPUT.FOLDER,UMAP.LAYOUT.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME)

print("Finished reading needed files")

## Combine expression_matrix_analysis.csv and filenums.csv to generate new Concat_Transformed.csv ==========================================

concat.transformed <- cbind(expression.matrix, file.nums)
colnames(concat.transformed) <- gsub("\\.", "_", colnames(concat.transformed))
colnames(concat.transformed) [41] <- "File"
fwrite(concat.transformed,file = CONCATTRANSFORMED.OUT.FILENAME,row.names = FALSE)

print("Finished renaming expression_matrix.csv to Concat_Transformed.csv")

## Rename umap_xy_RX.csv to umap_layout.csv ========================================================
#Rename columns and save output file
colnames(layout.in) <- "umap_x"
colnames(layout.in) [2] <- "umap_y"
fwrite(layout.in,file = UMAP.LAYOUT.OUT.FILENAME)

print("Finished renaming umap_xy_RX.csv to UMAP_layout.csv")

## Convert cluster_RX_assigns.csv to clusters.csv ==================================================
#Convert column to row
clusters.out <- t(clusters.in)
#Save output file
fwrite(clusters.out,file = CLUSTERS.OUT.FILENAME,col.names = FALSE,sep = ",")

print("Finished converting cluster_RX_assigns.csv to clusters.csv")

print("End 00_URD_Prep.R")
