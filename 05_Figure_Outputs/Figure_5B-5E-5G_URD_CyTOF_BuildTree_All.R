#Corey Williams, edited by Austin Keeler and Ashley Hirt, University of Virginia
#21 Feb, 2020, edited August 21st, 2020
#Run URD on CyTOF data, cholmod error fixed
#
#To use this code, run each section (denoted by = signs) separately 
#following instructions at top of section
#
#This process and QC is outlined in 
#https://github.com/farrellja/URD/blob/master/Analyses/QuickStart/URD-QuickStart-AxialMesoderm.md
#
#Other parameter troubleshooting can be found on Github through some googling
#
#Most run parameters are left as written in tutorial, so not optimized in our hands. 
#Be ready to troubleshoot
#

print("Start URD_CyTOF_BuildTree.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(umap)
library(URD)
library(dplyr)
library(scales)
library(tibble)
library(data.table)
library(ggplot2)
library(stringr)
library(rapport)
library(viridisLite)

print("Libraries loaded")


## Input parameters ===============================================================================
# Initial parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
URD.RDATA.IN <- "/URD_Final.RData" #Change filenames as relevant, section should run without issue

# Load URD.RData
load(paste0(INPUT.FOLDER, URD.RDATA.IN)) #.Rdata from URD run must have completed through Biased Walks (or beyond) to be used in this script


# Plot parameters to input after load URD Object
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
PLOT.BY.CLUSTER <- TRUE #TRUE exports dendrogram colored by cluster
CLUSTER.IDS.CSV <- "clusters_Subsampled.csv"
CLUSTER.COLOR.CSV <- "colors_R2.csv" # Specify .csv file containing column "color" with colors
BUILD.NEW.TREE <- FALSE #TRUE builds new tree, FALSE uses existing URD_Object.tree (i.e. for redos, faster)
OVERWRITE.RDATA <- FALSE #TRUE overwrites URD_Object.tree in .RData, FALSE keeps original tree
CELLS.PER.PSEUDOTIME.BIN <- 5 #Typically values between 5 and 40
BINS.PER.PSEUDOTIME.WINDOW <- 8 #Typically values between 5 and 20
TIP.NAMES <- c("Pax6pos Neuronal Progenitors","Pax6pos Glial Progenitors", "Pax6pos Tbr2 Neuronal Progenitors",
               "Mixed Interneuron Progenitors", "Interneuron Progenitors","Mixed Olig2 Neuronal Progenitors",
               "Ctip2 Interneuron Precursors", "Ctip2 Early Neurons", "VCAM Glial Precursors", 
               "Olig2 Glial Precursors","Olig2 Neuronal Precursors", "Astroglial Precursors",
               "Sox10 OPCs","DCX OPCs") # If not already in URD_Object
OUTPUT.BASENAME <- "URD_Tree_ColoredBy_" # Base name for plot outputs (e.g. "URD_Tree_ColoredBy_")
OUTPUT.DEVICE <- "pdf"

# Optional plot parameters
PLOT.BY.SEGMENT <- FALSE #TRUE exports dendrogram colored by segment
PLOT.BY.STAGE <- FALSE #TRUE exports dendrogram colored by stage
HIGHLIGHT.STAGE <- FALSE #TRUE exports individual dendrograms highlighting each stage 
HIGHLIGHT.CLUSTER <- TRUE #TRUE exports individual dendrograms highlighting each cluster
PLOT.BY.TISSUE <- FALSE #TRUE exports dendrogram colored by tissue (can specific colors below)
HIGHLIGHT.TISSUE <- FALSE #TRUE exports individual plots highlighting each tissue
STAGE.COLORS <- NULL #NULL uses default plasma color palette, otherwise specify colors as c("red","blue","yellow",...)
SEGMENT.COLORS <- NULL #NULL uses default rainbow color palette, otherwise specify colors as c("red","blue","yellow",...)
STAGE.COLORS <- c("E11"="#0D0887FF",
                  "E12"= "#3B049AFF",
                  "E13"= "#5D01A6FF",
                  "E14" = "#7E03A8FF",
                  "E15" = "#9C179EFF",
                  "E16" = "#B52F8CFF",
                  "E17" = "#CC4678FF",
                  "E18" = "#DE5F65FF",
                  "P0" = "#ED7953FF",
                  "P1" = "#F89441FF",
                  "P2" = "#FDB32FFF",
                  "P3" = "#FBD424FF",
                  "P4" = "#F0F921FF")

TISSUE.COLORS <- c("Brain" = "#F242F5",
                   "Cortex" = "#4E42F5",
                   "Diencephalon" = "#18C92D",
                   "Midbrain" = "#F0F013",
                   "Hindbrain" = "#F01313")

cluster.colors <- read.csv(CLUSTER.COLOR.CSV)
CLUSTER.COLORS <- cluster.colors$color
names(CLUSTER.COLORS) <- as.character(cluster.colors$cluster)
 

## Build tree =====================================================================================
if (BUILD.NEW.TREE == TRUE) {
  #Load the cells used for each tip into the URD object
  URD_Object.tree <- loadTipCells(URD_Object, "tip.clusters")
  #Build the tree
  num.tips <- max(URD_Object@group.ids[rownames(final.cells@group.ids), "tip.clusters"])
  URD_Object.tree <- buildTree(URD_Object.tree, pseudotime = "pseudotime", 
                               tips.use=1:num.tips, divergence.method = "preference", 
                               cells.per.pseudotime.bin = CELLS.PER.PSEUDOTIME.BIN, 
                               bins.per.pseudotime.window = BINS.PER.PSEUDOTIME.WINDOW, 
                               save.all.breakpoint.info = T, p.thresh=0.001)
  #Name tips
  if (is.null(TIP.NAMES) == TRUE) {
    tip.names <- sapply(1:num.tips,function(x){paste("Tip",x)}) 
  } else {
    tip.names <- TIP.NAMES
  }
  # Name the segments based on our previous determination of the identity of tips 1 and 2.
  URD_Object.tree <- nameSegments(URD_Object.tree, segments=as.character(1:num.tips), 
                                  segment.names = tip.names)
}

# Output dendrogram colored by segment =============================================================
if (PLOT.BY.SEGMENT == TRUE) {
  p.segment <- plotTree(URD_Object.tree, "segment", title="URD Tree Segment") # Tree colored by segment
  if (!is.null(SEGMENT.COLORS)) {
    segment.colors <- SEGMENT.COLORS
  } else {
    n.segments <- length(unique(URD_Object.tree@group.ids$segment)) # Uses number of stages from URD object
    segment.colors <- rainbow(n.segments, alpha = 1) # Calculates discrete rainbow palette
  }
  p.segment.final <- p.segment + scale_color_manual(values = segment.colors) # Add custom color scale to ssegment plot
  p.segment.filename <- paste0(OUTPUT.BASENAME, "Segment", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.segment.filename, plot = p.segment.final, device = OUTPUT.DEVICE) # Save plot output
}


# Calculate percent metadata per segment =======================================================
# Get all cell IDs and corresponding metadata 
all_segments <- URD_Object.tree@tree$cells.in.segment # Get list of lists with cells in each segment
all_tissues <- data.frame(Cell_ID = as.character(URD_Object.tree@group.ids$init, stringsAsFactors=FALSE)) # Get all cell IDs
all_tissues$Tissue <- URD_Object.tree@meta$tissue # Get tissue for all cells
unique_segments <- unique(URD_Object.tree@tree$segments) # Get unique segments (sometimes a number is skipped!)

# Initialize lists to collect counts
total.counts <- vector(mode = "list", length = length(unique_segments))
brain.counts <- vector(mode = "list", length = length(unique_segments))
cortex.counts <- vector(mode = "list", length = length(unique_segments))
dien.counts <- vector(mode = "list", length = length(unique_segments))
mid.counts <- vector(mode = "list", length = length(unique_segments))
hind.counts <- vector(mode = "list", length = length(unique_segments))

# Get cell counts for each tissue in each segment
for (i in unique_segments) {
  cells_segment <- as.vector(unlist(URD_Object.tree@tree$cells.in.segment[i])) # Isolate cells in segment
  tissues_segments <- all_tissues[match(cells_segment, all_tissues$Cell_ID),] # Obtain tissue for each cell
  total.counts[[i]] <- nrow(tissues_segments) # Count all cells in segment
  brain.counts[[i]] <- length(which(tissues_segments$Tissue == "Brain")) # Count brain cells in segment
  cortex.counts[[i]] <- length(which(tissues_segments$Tissue == "Cortex")) # Count cortex cells in segment
  dien.counts[[i]] <- length(which(tissues_segments$Tissue == "Diencephalon")) # Count dien cells in segment
  mid.counts[[i]] <- length(which(tissues_segments$Tissue == "Midbrain")) # Count midbrain cells in segment
  hind.counts[[i]] <- length(which(tissues_segments$Tissue == "Hindbrain")) # Count hindbrain cells in segment
}

# Combine counts into single data.frame
total.counts.out <- as.data.frame(unlist(total.counts))
brain.counts.out <- as.data.frame(unlist(brain.counts))
cortex.counts.out <- as.data.frame(unlist(cortex.counts))
dien.counts.out <- as.data.frame(unlist(dien.counts))
mid.counts.out <- as.data.frame(unlist(mid.counts))
hind.counts.out <- as.data.frame(unlist(hind.counts))
counts.out <- cbind(total.counts.out,brain.counts.out,cortex.counts.out,
                                 dien.counts.out,mid.counts.out,hind.counts.out)
colnames(counts.out) <- c("Total","Brain","Cortex","Diencephalon","Midbrain","Hindbrain")

# Generate pie charts showing tissue percentages per segment
for (i in 1:length(unique_segments)) {
  pie_df <- as.data.frame(t(counts.out[i,2:6]))
  colnames(pie_df) <- "cell_count"
  pie_df <- tibble::rownames_to_column(pie_df, "Tissue")
  pie_df <- pie_df %>% mutate(prop = cell_count / counts.out[i,1]) %>%
    mutate(Tissue_Percentage = paste0(Tissue," (", percent(prop),")")) %>%
    mutate(Tissue_Percentage = factor(x=Tissue_Percentage,levels=Tissue_Percentage))
  PIE.COLORS <- c("#F242F5","#4E42F5","#18C92D","#F0F013","#F01313")
  names(PIE.COLORS) <- pie_df$Tissue_Percentage
  p <- ggplot(pie_df, aes(x="", y=cell_count, fill=Tissue_Percentage)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) + 
    scale_fill_manual(values = PIE.COLORS) +
    theme_void()
  pie_filename <- paste0(INPUT.FOLDER,"/","PieChart_Segment_",i,".png")
  ggsave(pie_filename, height = 7, width = 7)
}

#Output .csv file with counts
counts.out <- tibble::rownames_to_column(counts.out, "Segment")
counts_filename <- paste0(INPUT.FOLDER,"/Tissue_Counts_Per_Segment.csv")
fwrite(counts.out, file = counts_filename)

# Output dendrogram colored by stage ===============================================================
if (PLOT.BY.STAGE == TRUE) {
  p.stage <- plotTree(URD_Object.tree, "stage", title="URD Stage", cell.alpha=0.25, cell.size=0.5) # Tree colored by stage
  cells.included <- rownames(p.stage$layers[[1]]$data) # Get ordered list of cells in p.stage plot
  cells.inc.IDs <- substring(cells.included, 2) # Removes first character "V" from strings
  cell.IDs <- as.vector(as.numeric(cells.inc.IDs)) # Convert cell IDs to vector
  stages.all <- URD_Object@meta$stage # Get stages for all cells from URD metadata
  stages.included <- stages.all[cell.IDs] # Subset stage data for included cells
  p.stage.final <- p.stage # Make a copy of p.stage plot to preserve original
  p.stage.final$layers[[1]]$data$expression <- stages.included # Add stage information to p.stage plot
  if (!is.null(STAGE.COLORS)) {
    stage.colors <- STAGE.COLORS # Use custom stage colors specified in parameters at top
  } else {
    n.stages <- length(unique(URD_Object.tree@meta$stage)) # Uses number of stages from URD object
    stage.colors <- plasma(n.stages, alpha = 1, begin = 0, end = 1, direction = 1) # Calculates discrete plasma palette
  }
  p.stage.final <- p.stage.final + scale_color_manual(values = stage.colors) # Add custom color scale to stage plot #NOT WORKING!
  p.stage.filename <- paste0(OUTPUT.BASENAME, "Stage", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.stage.filename, plot = p.stage.final, device = OUTPUT.DEVICE) # Save plot output
  #(Optional) Generate plots highlighting each stage
  if (HIGHLIGHT.STAGE == TRUE) { #NOT FUNCTIONAL YET!!!
    for (i in unique(URD_Object.tree@meta$stage)) {
      p.stage.final.hl <- plotTreeHighlight(URD_Object.tree, "stage", i) # Plot only cells in selected stages in red
      p.stage.filename.hl <- paste0(OUTPUT.BASENAME, "Stage_", i, ".", OUTPUT.DEVICE) # Generate output filename
      ggsave(p.stage.filename.hl, plot = p.stage.final.hl, device = OUTPUT.DEVICE) # Save plot output
    }
  }
}

## Output dendrogram colored by cluster ============================================================
if (PLOT.BY.CLUSTER == TRUE) {
  p.cluster <- plotTree(URD_Object.tree, label.type="meta", "cluster", title="Cluster") # Tree colored by cluster
  cells.included <- rownames(p.cluster$layers[[1]]$data) # Get ordered list of cells included in dendrogram
  cells.inc.IDs <- substring(cells.included, 2) # Removes first character "V" from strings
  cell.IDs <- as.vector(as.numeric(cells.inc.IDs)) # Convert cell IDs to vector
  clusters.all <- URD_Object@meta$cluster # Get cluster for all cells from URD metadata
  clusters.included <- clusters.all[cell.IDs] # Subset cluster data for included cells
  p.cluster.final <- p.cluster # Make a copy of p.cluster
  p.cluster.final$layers[[1]]$data$expression <- clusters.included # Add cluster data to plot layer 1 expression
  cluster.colors <- CLUSTER.COLORS
  #color.vec <- sort(unique(URD_Object.tree@meta$cluster)) # Get cluster numbers included in plot
  #cluster.colors <- cluster.colors.all[color.vec] # Get only colors needed for clusters in plot
  p.cluster.final <- p.cluster.final + scale_color_manual(values = cluster.colors) # Add custom color scale to stage plot 
  p.cluster.filename <- paste0(OUTPUT.BASENAME, "Cluster", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.cluster.filename, plot = p.cluster.final, device = OUTPUT.DEVICE) # Save plot output
}


if (HIGHLIGHT.CLUSTER == TRUE) {
  cluster.nums <- unique(URD_Object.tree@meta$cluster)
  for (i in cluster.nums) {
  p.cluster <- plotTree(URD_Object.tree, label.type="meta", "cluster", title="Cluster", cell.alpha = 1) # Tree colored by cluster
  cells.included <- rownames(p.cluster$layers[[1]]$data) # Get ordered list of cells included in dendrogram
  cells.inc.IDs <- substring(cells.included, 2) # Removes first character "V" from strings
  cell.IDs <- as.vector(as.numeric(cells.inc.IDs)) # Convert cell IDs to vector
  clusters.all <- URD_Object@meta$cluster # Get cluster for all cells from URD metadata
  clusters.included <- clusters.all[cell.IDs] # Subset cluster data for included cells
  p.cluster.final <- p.cluster # Make a copy of p.cluster
  p.cluster.final$layers[[1]]$data$expression <- clusters.included # Add cluster data to plot layer 1 expression
  p.cluster.i <- p.cluster.final
  p.cluster.i$layers[[1]]$data <- subset(p.cluster.i$layers[[1]]$data, expression == i) #Subset for cluster
  colors.i <- CLUSTER.COLORS # make copy of custom cluster colors
  colors.i[which(colors.i!=i)] <- "#ececec"# change all colors to background grey
  colors.i <- replace(colors.i,i,'red') # Indicate highlight color for cluster i
  p.cluster.i <- p.cluster.i + scale_color_manual(values = colors.i)
  p.cluster.i.filename <- paste0(OUTPUT.BASENAME, "Cluster_",i, ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.cluster.i.filename, plot = p.cluster.i, device = OUTPUT.DEVICE) # Save plot output
  }
}


#Should work but doesn't
#if (HIGHLIGHT.CLUSTER == TRUE) {
  #cluster.nums <- unique(URD_Object.tree@meta$cluster)
  #for (i in cluster.nums) {
    #p.cluster.highlight <- plotTreeHighlight(URD_Object.tree, label.type = "meta", label.name="cluster", label.value = i,
                                             #color = "red", bg.color = "#CCCCCC", highlight.alpha = 1, highlight.size = 10)
    #p.cluster.highlight.filename <- paste0(OUTPUT.BASENAME, "Cluster_",i,".", OUTPUT.DEVICE) # Generate output filename
    #ggsave(p.cluster.highlight.filename, plot = p.cluster.highlight, device = OUTPUT.DEVICE) # Save plot output
  #}
  
#}


## Output dendrogram colored by tissue ============================================================
if (PLOT.BY.TISSUE == TRUE) {
  p.tissue <- plotTree(URD_Object.tree, label.type = "meta", label = "tissue", title="Tissue") # Tree colored by tissue
  cells.included <- rownames(p.tissue$layers[[1]]$data) # Get ordered list of cells in p.tissue plot
  cells.inc.IDs <- substring(cells.included, 2) # Removes first character "V" from strings
  cell.IDs <- as.vector(as.numeric(cells.inc.IDs)) # Convert cell IDs to vector
  tissue.all <- URD_Object@meta$tissue # Get tissue info for all cells from URD metadata
  tissue.included <- tissue.all[cell.IDs] # Subset tissue info for included cells
  p.tissue.final <- p.tissue # Make a copy of p.tissue plot to preserve original
  p.tissue.final$layers[[1]]$data$expression <- tissue.included # Add tissue info to p.tissue plot
  if (!is.null(TISSUE.COLORS)) {
    tissue.colors <- TISSUE.COLORS # Use custom tissue colors specified in parameters at top
  } else {
    n.tissue <- length(unique(URD_Object.tree@meta$tissue)) # Uses number of tissue from URD object
    tissue.colors <- rainbow(n.tissue, alpha = 1, begin = 0, end = 1, direction = 1) # Calculates discrete rainbow palette
  }
  p.tissue.final <- p.tissue.final + scale_color_manual(values = tissue.colors) # Add custom color scale to stage plot #NOT WORKING!
  p.tissue.filename <- paste0(OUTPUT.BASENAME, "Tissue", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.tissue.filename, plot = p.tissue.final, device = OUTPUT.DEVICE) # Save plot output
  
  if (HIGHLIGHT.TISSUE == TRUE) {
    for (i in unique(URD_Object.tree@meta$tissue)) {
      p.tissue <- plotTree(URD_Object.tree, label.type = "meta", label = "tissue", title="Tissue") # Tree colored by tissue
      cells.included <- rownames(p.tissue$layers[[1]]$data) # Get ordered list of cells in p.tissue plot
      cells.inc.IDs <- substring(cells.included, 2) # Removes first character "V" from strings
      cell.IDs <- as.vector(as.numeric(cells.inc.IDs)) # Convert cell IDs to vector
      tissue.all <- URD_Object@meta$tissue # Get tissue info for all cells from URD metadata
      tissue.included <- tissue.all[cell.IDs] # Subset tissue info for included cells
      p.tissue.final <- p.tissue # Make a copy of p.tissue plot to preserve original
      p.tissue.final$layers[[1]]$data$expression <- tissue.included # Add tissue info to p.tissue plot
      p.tissue.i <- p.tissue.final # Make copy of p.tissue.final plot
      p.tissue.i$layers[[1]]$data <- subset(p.tissue.i$layers[[1]]$data, expression == i) #Subset for tissue
      p.tissue.i <- p.tissue.i + scale_color_manual(values = )
      p.tissue.i.filename <- paste0(OUTPUT.BASENAME, "Tissue_",i, ".", OUTPUT.DEVICE) # Generate output filename
      ggsave(p.tissue.i.filename, plot = p.tissue.i, device = OUTPUT.DEVICE) # Save plot output
    }
  }
}

## Save R workspace image =========================================================================
if (OVERWRITE.RDATA == TRUE) {
  save.image("URD_Final.RData")
  #reset working directory
  setwd("..")
  print("Finished URD_CyTOF_BuildTree.R")
} else {
  print("Finished URD_CyTOF_BuildTree.R")
}

print("Finish URD_CyTOF_BuildTree.R")