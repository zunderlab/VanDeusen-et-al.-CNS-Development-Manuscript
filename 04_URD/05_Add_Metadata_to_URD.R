#Amy Van Deusen, University of Virginia
#2021
## Add sample metadata (e.g. clusters, tissue types, etc.) to existing URD Object

print("Start Add_Metadata_to_URD.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(dplyr)
library(stringr)
library(rapport)

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
OVERWRITE.URD.OBJECT <- TRUE #TRUE saves URD_Object.new and updates existing URD_Object, FALSE = saves new object only 
METADATA.CSV <- "metadata.csv"
SUBSAMPLE.IDS.CSV <- "Subsample_IDs.csv"
ADD.FILE.NAMES <- TRUE #TRUE reimports file names using subsample IDs
FILENUMS.CSV <- "filenums_Subsampled.csv" #Must specify .csv if ADD.FILE.NAMES = TRUE, otherwise = NULL
ADD.STAGE.INFO <- TRUE #TRUE reimports stage info using file names 
ADD.CLUSTER.INFO <- TRUE #TRUE adds cluster info, requires a .csv file (below)
CLUSTER.IDS.CSV <- "clusters_Subsampled.csv" #Must specify .csv if ADD.CLUSTER.INFO = TRUE, otherwise = NULL
ADD.TISSUE.INFO <- TRUE #TRUE adds tissue info using filenames

# Read in necessary files ====================================================================
sample.md.filename <- paste0(INPUT.FOLDER,"/",METADATA.CSV)
sample.md <- read.csv(sample.md.filename, header = TRUE)
subsample.ids.filename <- paste0(INPUT.FOLDER,"/",SUBSAMPLE.IDS.CSV)
subsample.ids <- as.numeric(unlist(read.csv(subsample.ids.filename, header=FALSE,stringsAsFactors=FALSE)))
if (ADD.FILE.NAMES == TRUE) {
  filenums.filename <- paste0(INPUT.FOLDER,"/",FILENUMS.CSV)
  filenums <- read.table(filenums.filename, header=TRUE, sep=",", check.names=FALSE)[,]
}
if (ADD.CLUSTER.INFO == TRUE) {
  cluster.id.filename <- paste0(INPUT.FOLDER, "/", CLUSTER.IDS.CSV)
  cluster.ids <- read.csv(cluster.id.filename, header = FALSE)
}

## Load R Data ===============================================================================
# Load URD.RData
URD.RDATA.IN <- "URD_Final.RData" #Change filenames as relevant, section should run without issue
load(paste0(INPUT.FOLDER,"/",URD.RDATA.IN)) #.Rdata from URD run must have completed through Biased Walks (or beyond) to be used in this script
URD_Object.new <- URD_Object # Make a copy of existing URD_Object (NOTE: can replace with URD_Object.tree)

# Add file names to URD_Object.new@meta$file_name ===========================================
if (ADD.FILE.NAMES == TRUE) {
  filenames_key <- as.character(sample.md[,"file_name"]) # Get list of file names
  names(filenames_key) <- 1:length(sample.md[,"file_name"]) # Add index to list of file names
  filenames.out <- filenames_key[filenums] # Obtain file names for all cells
  filenames.df <- as.data.frame(filenames.out) # Convert filenames to a dataframe (needed for next line to work)
  filenames.vec <- filenames.df[,1] # Convert filenames to a character vector
  URD_Object.new@meta$file_name <- filenames.vec # Add filenames to URD_Object.new
  urd_meta$file_name <- filenames.vec # Add filenames to urd_meta (to match URD_Object)
}

# Add stage information to URD_Object.new@meta$stage ========================================
if (ADD.STAGE.INFO == TRUE) {
  urd_meta_stage <- urd_meta %>%
    mutate(stage = case_when(
      str_detect(file_name, "E11") ~ "E11",
      str_detect(file_name, "E12") ~ "E12",
      str_detect(file_name, "E13") ~ "E13",
      str_detect(file_name, "E14") ~ "E14",
      str_detect(file_name, "E15") ~ "E15",
      str_detect(file_name, "E16") ~ "E16",
      str_detect(file_name, "E17") ~ "E17",
      str_detect(file_name, "E18") ~ "E18",
      str_detect(file_name, "P0") ~ "P0",
      str_detect(file_name, "P1") ~ "P1",
      str_detect(file_name, "P2") ~ "P2",
      str_detect(file_name, "P3") ~ "P3",
      str_detect(file_name, "P4") ~ "P4"))
  urd_meta_stage[, 'stage'] <- as.factor(urd_meta_stage[, 'stage']) # Make stages a factor
  URD_Object.new@meta <- urd_meta_stage # Add stages to existing URD_Object
  urd_meta <- urd_meta_stage # Add stages to urd_meta (to match URD_Object)
}
  
# Add cluster information to URD_Object.new@meta$cluster ===================================
if (ADD.CLUSTER.INFO == TRUE) {
  colnames(cluster.ids) <- "cluster" # Rename column
  cluster.ids[, 'cluster'] <- as.factor(cluster.ids[, 'cluster']) # Make stages into factor levels
  cluster.ids.vec <- cluster.ids[,1] # Convert cluster column into vector
  URD_Object.new@meta$cluster <- cluster.ids.vec # Add clusters to the existing URD_Object
  urd_meta$cluster <- cluster.ids.vec # Add clusters to urd_meta (to match URD_Object)
}

#Add tissue information to URD_Object.
if (ADD.TISSUE.INFO == TRUE) {
  urd_meta_tissue <- urd_meta %>%
    mutate(tissue = case_when(
      str_detect(file_name, "Brain") ~ "Brain",
      str_detect(file_name, "Cortex") ~ "Cortex",
      str_detect(file_name, "Diencephalon") ~ "Diencephalon",
      str_detect(file_name, "Midbrain") ~ "Midbrain",
      str_detect(file_name, "Hindbrain") ~ "Hindbrain"
    ))
  urd_meta_tissue[,'tissue'] <- as.factor(urd_meta_tissue[,'tissue']) # Convert tissues into factor levels
  urd_meta_tissue$tissue <- factor(urd_meta_tissue$tissue,  levels = c("Brain","Cortex","Diencephalon","Midbrain","Hindbrain"))
  URD_Object.new@meta$tissue <- urd_meta_tissue$tissue # Add new data to existing URD_Object
  urd_meta$tissue <- urd_meta_tissue$tissue # Add tissues to urd_meta (to match URD_Object)
}

## Save R workspace image =========================================================================
if (OVERWRITE.URD.OBJECT == TRUE) {
  URD_Object <- URD_Object.new
  save.image("URD_Final.RData")
  print("Existing URD Object Updated. Finished Repair_URD_Metadata.R")
} else {
  print("New URD Object Created, Existing URD Object Not Updated. Finished Repair_URD_Metadata.R")
}

print("Finish Add_Metadata_to_URD.R")