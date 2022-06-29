#Corey Williams, University of Virginia
#17 May, 2019
#Subsample .fcs file and output .csv file containing all data points asinh transformed

print("Start 02_FCS_to_CSV.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(flowCore)
library(data.table)

print("Libraries loaded")

## Input parameters ===============================================================================

METADATA_FILENAME <- "metadata.csv"
PANEL_FILENAME <- "panel.csv"
SUBSAMPLE_IDXS_FILENAME <- "subsample_idxs.csv" #unused if no subsampling
EXPRS_MAT_ANALYSIS_FILENAME <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER_FILENAME <- "expression_matrix_other.csv"
FILENUMS_FILENAME <- "filenums.csv"
ASINH_FACTOR <- NULL #set to NULL if set in panel
SUBSAMPLE_BOOL <- FALSE
NUM_SUBSAMPLES <- NULL #change this value if you want subsampling. Number of cells per file

# Get path to fcs files
args = commandArgs(trailingOnly=TRUE) # Get command line arguments
if (length(args) == 0) {
  INPUT_FOLDER <- getwd() # If there's no command line argument, set to current directory
} else {
  INPUT_FOLDER <- args[1] # Use filepath specified by command line argument
}
# Get path to where you want to save csv files
if (length(args) < 2) {
  OUTPUT_FOLDER <- INPUT_FOLDER # If there's no command line argument, make this the same as FCS_FILEPATH
} else {
  OUTPUT_FOLDER <- args[2] # Use filepath specified by command line argument
}

print("Got input parameters, loading fcs files")

## Read metadata ==================================================================================
metadata_path <- paste(INPUT_FOLDER, METADATA_FILENAME, sep="/")
metadata <- read.table(metadata_path, header=TRUE, sep=",", check.names=FALSE)

## Read fcs files =================================================================================
fs <- read.flowSet(metadata$File_Name, path=INPUT_FOLDER, transformation=FALSE,
                   truncate_max_range=FALSE)

print("Loaded fcs files")

## Read panel =====================================================================================
panel_path <- paste(INPUT_FOLDER, PANEL_FILENAME, sep="/")
panel <- read.table(panel_path, header=TRUE, sep=",", check.names=FALSE)

## Pull out markers for analysis ==================================================================
analysis_markers <- panel$Analysis

## Get things set up for subsampling (if this option's selected) ==================================
if (SUBSAMPLE_BOOL){
  # Get subsample id's for each file
  subsample_list <- lapply(seq_along(fs), function(x){sample(1:nrow(fs[[x]]), NUM_SUBSAMPLES)})
  subsample_idxs <- do.call(cbind, subsample_list)
  # Add filenames for clarity
  colnames(subsample_idxs)  <- metadata$File_Name
  # Write to file
  fwrite(subsample_idxs, file=SUBSAMPLE_IDXS_FILENAME, col.names=FALSE, sep=",")
  print("Saved subsample IDX file")
  # Define exprs retrieve function to include subsampling
  retrieve_exprs <- function(x){
    exprs(fs[[x]])[subsample_idxs[, x], ]
  }
} else {
  # Define exprs retrieve function to NOT include subsampling
  retrieve_exprs <- function(x){
    exprs(fs[[x]])
  }
}
  
## Concatenate exprs arrays (with or without subsampling as defined above =========================
## using the function retrieve_exprs() ).  Also, store filenums array.
data_mat <- c()
filenums <- c()
for (i in seq_along(fs)) {
  fs_i_mat <- retrieve_exprs(i)
  data_mat <- rbind(data_mat, fs_i_mat)
  filenums <- c(filenums, rep(i, nrow(fs_i_mat)))
}
colnames(data_mat) <- panel$Fixed_Param

## Perform transforms on expression matrix ========================================================
# Divide by asinh factors (should be 1 for linear transforms)
data_mat <- t(apply(data_mat, 1, function(x) {x/panel$Factor}))
# Perform asinh transform
asinh_channels <- which(panel$Transform=="Asinh")
data_mat[,asinh_channels] <- asinh(data_mat[,asinh_channels])

## Save concatenated exprs arrays (one for analysis, then one for everything else) ================
fwrite(data_mat[, which(panel$Analysis==1)],
       file=EXPRS_MAT_ANALYSIS_FILENAME, row.names=FALSE)
fwrite(data_mat[, which(panel$Analysis!=1)],
       file=EXPRS_MAT_OTHER_FILENAME, row.names=FALSE)
print("Expression matrix files outputted")

fwrite(matrix(filenums, ncol=1), file=FILENUMS_FILENAME, row.names=FALSE)
print("Filenums reference file outputted")

print("Finish 02_FCS_to_CSV.R")
