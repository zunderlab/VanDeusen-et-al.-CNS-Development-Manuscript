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
library(data.table)
library(ggplot2)
library(stringr)
library(rapport)
library(viridisLite)

print("Libraries loaded")


## Input parameters ===============================================================================
# Load R data
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
URD.RDATA.IN <- "/URD_Final.RData" #Change filenames as relevant, section should run without issue
load(paste0(INPUT.FOLDER, URD.RDATA.IN)) #.Rdata from URD run must have completed through Biased Walks (or beyond) to be used in this script

# Initial parameters
PLOT.BY.CLUSTER <- FALSE #TRUE exports dendrogram colored by cluster
CLUSTER.IDS.CSV <- "clusters_Subsampled.csv"
CLUSTER.COLOR.CSV <- "colors_R2.csv" # NULL = default colors, otherqise specify .csv file containing column "color" with colors
BUILD.NEW.TREE <- FALSE #TRUE builds new tree, FALSE uses existing URD_Object.tree (i.e. for redos, faster)
OVERWRITE.RDATA <- FALSE #TRUE overwrites URD_Object.tree in .RData, FALSE keeps original tree
CELLS.PER.PSEUDOTIME.BIN <- 5 #Typically values between 5 and 40
BINS.PER.PSEUDOTIME.WINDOW <- 8 #Typically values between 5 and 20
TIP.NAMES <- c("Pax6 Glial Progenitors", "OPCs", "Astroglial Precursors", 
               "Olig2 Glial Precursors",
               "Ctip2 Interneurons", "Tbr1 Neurons", "Interneurons", 
               "Tbr1Ctip2 Neurons", "Ly6C Endothelial Cells", "Fibroblasts",
               "A2B5 Endothelial Cells", "Nonneural Cells") # If not already in URD_Object
OUTPUT.BASENAME <- "URD_Tree_ColoredBy_" # Base name for plot outputs (e.g. "URD_Tree_ColoredBy_")
OUTPUT.DEVICE <- "pdf"

# Optional plot parameters
PLOT.BY.MARKER <- TRUE #TRUE exports dendrograms colored by each marker
PLOT.BY.SEGMENT <- FALSE #TRUE exports dendrogram colored by segment
PLOT.BY.STAGE <- FALSE #TRUE exports dendrogram colored by stage
HIGHLIGHT.STAGE <- FALSE #TRUE exports individual dendrograms highlighting each stage 
HIGHLIGHT.CLUSTER <- FALSE #TRUE exports individual dendrograms highlighting each cluster
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
                   "Cortex" = "#4E42F5")

# Load necessary files before URD Object (to prevent errors)
if (PLOT.BY.CLUSTER == TRUE) {
  cluster_id_filename <- paste0(INPUT.FOLDER, "/", CLUSTER.IDS.CSV)
  cluster_ids <- read.csv(cluster_id_filename)
  if (!is.null(CLUSTER.COLOR.CSV)) {
    clusters.df <- as.data.frame(read.csv(CLUSTER.COLOR.CSV)) # Read in cluster colors from csv
    cluster.colors.ll <- as.vector(clusters.df$color) # Create vector with cluster colors
  } else {
    n.clusters <- length(unique(URD_Object.tree@meta$cluster))
    cluster.colors.all <- rainbow(n.clusters, alpha = 1) # Apply default rainbow palette to dendrogram
  }
}

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

# Output dendrograms colored by each marker
if (PLOT.BY.MARKER == TRUE) {
  OUTPUT.BASENAME <- "URD_Tree_by_"
  p.marker.TuJ1 <- plotTree(URD_Object.tree, "TuJ1_89Y", title="TuJ1")
  #p.marker.TuJ1 <- p.marker.TuJ1 + scale_color_viridis() # Apply stadnard iridis color scale
  p.marker.TuJ1 <- p.marker.TuJ1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.TuJ1.filename <- paste0(OUTPUT.BASENAME,"TuJ1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.TuJ1.filename, plot = p.marker.TuJ1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Olig2 <- plotTree(URD_Object.tree, "Olig2_113In", title="Olig2")
  #p.marker.Olig2 <- p.marker.Olig2 + scale_color_viridis() # Apply viridis color scale
  p.marker.Olig2 <- p.marker.Olig2 +scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Olig2.filename <- paste0(OUTPUT.BASENAME,"Olig2", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Olig2.filename, plot = p.marker.Olig2, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.A2B5 <- plotTree(URD_Object.tree, "A2B5_115In", title="A2B5")
  #p.marker.A2B5 <- p.marker.A2B5 + scale_color_viridis() # Apply viridis color scale
  p.marker.A2B5 <- p.marker.A2B5 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.A2B5.filename <- paste0(OUTPUT.BASENAME,"A2B5", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.A2B5.filename, plot = p.marker.A2B5, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.PDGFRb <- plotTree(URD_Object.tree, "PDGFRb_141Pr", title="PDGFRb")
  #p.marker.PDGFRb <- p.marker.PDGFRb + scale_color_viridis() # Apply viridis color scale
  p.marker.PDGFRb <- p.marker.PDGFRb + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.PDGFRb.filename <- paste0(OUTPUT.BASENAME,"PDGFRb", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.PDGFRb.filename, plot = p.marker.PDGFRb, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.VCAM1 <- plotTree(URD_Object.tree, "VCAM1_142Nd", title="VCAM1")
  #p.marker.VCAM1 <- p.marker.VCAM1 + scale_color_viridis() # Apply viridis color scale
  p.marker.VCAM1 <- p.marker.VCAM1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.VCAM1.filename <- paste0(OUTPUT.BASENAME,"VCAM1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.VCAM1.filename, plot = p.marker.VCAM1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.PECAM1 <- plotTree(URD_Object.tree, "PECAM1_143Nd", title="PECAM1")
  #p.marker.PECAM1 <- p.marker.PECAM1 + scale_color_viridis() # Apply viridis color scale
  p.marker.PECAM1 <- p.marker.PECAM1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.PECAM1.filename <- paste0(OUTPUT.BASENAME,"PECAM1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.PECAM1.filename, plot = p.marker.PECAM1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Nestin <- plotTree(URD_Object.tree, "Nestin_144Nd", title="Nestin")
  #p.marker.Nestin <- p.marker.Nestin + scale_color_viridis() # Apply viridis color scale
  p.marker.Nestin <- p.marker.Nestin + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Nestin.filename <- paste0(OUTPUT.BASENAME,"Nestin", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Nestin.filename, plot = p.marker.Nestin, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Sox1 <- plotTree(URD_Object.tree, "Sox1_145Nd", title="Sox1")
  #p.marker.Sox1 <- p.marker.Sox1 + scale_color_viridis() # Apply viridis color scale
  p.marker.Sox1 <- p.marker.Sox1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Sox1.filename <- paste0(OUTPUT.BASENAME,"Sox1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Sox1.filename, plot = p.marker.Sox1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Tbr2 <- plotTree(URD_Object.tree, "Tbr2_146Nd", title="Tbr2")
  #p.marker.Tbr2 <- p.marker.Tbr2 + scale_color_viridis() # Apply viridis color scale
  p.marker.Tbr2 <- p.marker.Tbr2 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Tbr2.filename <- paste0(OUTPUT.BASENAME,"Tbr2", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Tbr2.filename, plot = p.marker.Tbr2, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.PDGFRa <- plotTree(URD_Object.tree, "PDGFRa_147Sm", title="PDGFRa")
  #p.marker.PDGFRa <- p.marker.PDGFRa + scale_color_viridis() # Apply viridis color scale
  p.marker.PDGFRa <- p.marker.PDGFRa + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.PDGFRa.filename <- paste0(OUTPUT.BASENAME,"PDGFRa", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.PDGFRa.filename, plot = p.marker.PDGFRa, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Prominin <- plotTree(URD_Object.tree, "Prominin_148Nd", title="Prominin")
  #p.marker.Prominin <- p.marker.Prominin + scale_color_viridis() # Apply viridis color scale
  p.marker.Prominin <- p.marker.Prominin + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Prominin.filename <- paste0(OUTPUT.BASENAME,"Prominin", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Prominin.filename, plot = p.marker.Prominin, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.CD45 <- plotTree(URD_Object.tree, "CD45_149Sm", title="CD45")
  #p.marker.CD45 <- p.marker.CD45 + scale_color_viridis() # Apply viridis color scale
  p.marker.CD45 <- p.marker.CD45 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.CD45.filename <- paste0(OUTPUT.BASENAME,"CD45", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.CD45.filename, plot = p.marker.CD45, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.NeuN <- plotTree(URD_Object.tree, "NeuN_150Nd", title="NeuN")
  #p.marker.NeuN <- p.marker.NeuN + scale_color_viridis() # Apply viridis color scale
  p.marker.NeuN <- p.marker.NeuN + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.NeuN.filename <- paste0(OUTPUT.BASENAME,"NeuN", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.NeuN.filename, plot = p.marker.NeuN, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Sox10 <- plotTree(URD_Object.tree, "Sox10_151Eu", title="Sox10")
  #p.marker.Sox10 <- p.marker.Sox10 + scale_color_viridis() # Apply viridis color scale
  p.marker.Sox10 <- p.marker.Sox10 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Sox10.filename <- paste0(OUTPUT.BASENAME,"Sox10", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Sox10.filename, plot = p.marker.Sox10, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Ki67 <- plotTree(URD_Object.tree, "Ki67_152Sm", title="Ki67")
  #p.marker.Ki67 <- p.marker.Ki67 + scale_color_viridis() # Apply viridis color scale
  p.marker.Ki67 <- p.marker.Ki67 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Ki67.filename <- paste0(OUTPUT.BASENAME,"Ki67", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Ki67.filename, plot = p.marker.Ki67, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.OligoO4 <- plotTree(URD_Object.tree, "OligoO4_153Eu", title="OligoO4")
  #p.marker.OligoO4 <- p.marker.OligoO4 + scale_color_viridis() # Apply viridis color scale
  p.marker.OligoO4 <- p.marker.OligoO4 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.OligoO4.filename <- paste0(OUTPUT.BASENAME,"OligoO4", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.OligoO4.filename, plot = p.marker.OligoO4, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Pax6 <- plotTree(URD_Object.tree, "Pax6_154Sm", title="Pax6")
  #p.marker.Pax6 <- p.marker.Pax6 + scale_color_viridis() # Apply viridis color scale
  p.marker.Pax6 <- p.marker.Pax6 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Pax6.filename <- paste0(OUTPUT.BASENAME,"Pax6", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Pax6.filename, plot = p.marker.Pax6, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.MCAM <- plotTree(URD_Object.tree, "MCAM_155Gd", title="MCAM")
  p.marker.MCAM <- p.marker.MCAM + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  #p.marker.MCAM <- p.marker.MCAM + scale_color_viridis() # Apply viridis color scale
  p.marker.MCAM.filename <- paste0(OUTPUT.BASENAME,"MCAM", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.MCAM.filename, plot = p.marker.MCAM, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.NeuroD1 <- plotTree(URD_Object.tree, "NeuroD1_156Gd", title="NeuroD1")
  #p.marker.NeuroD1 <- p.marker.NeuroD1 + scale_color_viridis() # Apply viridis color scale
  p.marker.NeuroD1 <- p.marker.NeuroD1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.NeuroD1.filename <- paste0(OUTPUT.BASENAME,"NeuroD1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.NeuroD1.filename, plot = p.marker.NeuroD1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.CD11b <- plotTree(URD_Object.tree, "CD11b_158Gd", title="CD11b")
  #p.marker.CD11b <- p.marker.CD11b + scale_color_viridis() # Apply viridis color scale
  p.marker.CD11b <- p.marker.CD11b + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.CD11b.filename <- paste0(OUTPUT.BASENAME,"CD11b", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.CD11b.filename, plot = p.marker.CD11b, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.CD24 <- plotTree(URD_Object.tree, "CD24_159Tb", title="CD24")
  #p.marker.CD24 <- p.marker.CD24 + scale_color_viridis() # Apply viridis color scale
  p.marker.CD24 <- p.marker.CD24 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.CD24.filename <- paste0(OUTPUT.BASENAME,"CD24", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.CD24.filename, plot = p.marker.CD24, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Sox2 <- plotTree(URD_Object.tree, "Sox2_160Gd", title="Sox2")
  #p.marker.Sox2 <- p.marker.Sox2 + scale_color_viridis() # Apply viridis color scale
  p.marker.Sox2 <- p.marker.Sox2 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Sox2.filename <- paste0(OUTPUT.BASENAME,"Sox2", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Sox2.filename, plot = p.marker.Sox2, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.GAD65 <- plotTree(URD_Object.tree, "GAD65_162Dy", title="GAD65")
  #p.marker.GAD65 <- p.marker.GAD65 + scale_color_viridis() # Apply viridis color scale
  p.marker.GAD65 <- p.marker.GAD65 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.GAD65.filename <- paste0(OUTPUT.BASENAME,"GAD65", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.GAD65.filename, plot = p.marker.GAD65, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.DCX <- plotTree(URD_Object.tree, "Dcx_163Dy", title="DCX")
  #p.marker.DCX <- p.marker.DCX + scale_color_viridis() # Apply viridis color scale
  p.marker.DCX <- p.marker.DCX + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.DCX.filename <- paste0(OUTPUT.BASENAME,"DCX", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.DCX.filename, plot = p.marker.DCX, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.MAP2 <- plotTree(URD_Object.tree, "MAP2_164Dy", title="MAP2")
  #p.marker.MAP2 <- p.marker.MAP2 + scale_color_viridis() # Apply viridis color scale
  p.marker.MAP2 <- p.marker.MAP2 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.MAP2.filename <- paste0(OUTPUT.BASENAME,"MAP2", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.MAP2.filename, plot = p.marker.MAP2, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.ALDH1A1 <- plotTree(URD_Object.tree, "ALDH1A1_165Ho", title="ALDH1A1")
  #p.marker.ALDH1A1 <- p.marker.ALDH1A1 + scale_color_viridis() # Apply viridis color scale
  p.marker.ALDH1A1 <- p.marker.ALDH1A1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.ALDH1A1.filename <- paste0(OUTPUT.BASENAME,"ALDH1A1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.ALDH1A1.filename, plot = p.marker.ALDH1A1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.GLAST <- plotTree(URD_Object.tree, "GLAST_167Er", title="GLAST")
  #p.marker.GLAST <- p.marker.GLAST + scale_color_viridis() # Apply viridis color scale
  p.marker.GLAST <- p.marker.GLAST + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.GLAST.filename <- paste0(OUTPUT.BASENAME,"GLAST", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.GLAST.filename, plot = p.marker.GLAST, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Ctip2 <- plotTree(URD_Object.tree, "Ctip2_168Er", title="Ctip2")
  #p.marker.Ctip2 <- p.marker.Ctip2 + scale_color_viridis() # Apply viridis color scale
  p.marker.Ctip2 <- p.marker.Ctip2 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Ctip2.filename <- paste0(OUTPUT.BASENAME,"Ctip2", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Ctip2.filename, plot = p.marker.Ctip2, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.GFAP <- plotTree(URD_Object.tree, "GFAP_169Tm", title="GFAP")
  #p.marker.GFAP <- p.marker.GFAP + scale_color_viridis() # Apply viridis color scale
  p.marker.GFAP <- p.marker.GFAP + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.GFAP.filename <- paste0(OUTPUT.BASENAME,"GFAP", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.GFAP.filename, plot = p.marker.GFAP, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Cux1 <- plotTree(URD_Object.tree, "Cux1_170Er", title="Cux1")
  #p.marker.Cux1 <- p.marker.Cux1 + scale_color_viridis() # Apply viridis color scale
  p.marker.Cux1 <- p.marker.Cux1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Cux1.filename <- paste0(OUTPUT.BASENAME,"Cux1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Cux1.filename, plot = p.marker.Cux1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Tbr1 <- plotTree(URD_Object.tree, "Tbr1_171Yb", title="Tbr1")
  #p.marker.Tbr1 <- p.marker.Tbr1 + scale_color_viridis() # Apply viridis color scale
  p.marker.Tbr1 <- p.marker.Tbr1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Tbr1.filename <- paste0(OUTPUT.BASENAME,"Tbr1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Tbr1.filename, plot = p.marker.Tbr1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.BLBP <- plotTree(URD_Object.tree, "BLBP_172Yb", title="BLBP")
  #p.marker.BLBP <- p.marker.BLBP + scale_color_viridis() # Apply viridis color scale
  p.marker.BLBP <- p.marker.BLBP + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.BLBP.filename <- paste0(OUTPUT.BASENAME,"BLBP", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.BLBP.filename, plot = p.marker.BLBP, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.TrkB <- plotTree(URD_Object.tree, "TrkB_174Yb", title="TrkB")
  #p.marker.TrkB <- p.marker.TrkB + scale_color_viridis() # Apply viridis color scale
  p.marker.TrkB <- p.marker.TrkB + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.TrkB.filename <- paste0(OUTPUT.BASENAME,"TrkB", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.TrkB.filename, plot = p.marker.TrkB, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.F480 <- plotTree(URD_Object.tree, "F4_80_175Lu", title="F4/80")
  #p.marker.F480 <- p.marker.F480 + scale_color_viridis() # Apply viridis color scale
  p.marker.F480 <- p.marker.F480 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.F480.filename <- paste0(OUTPUT.BASENAME,"F4_80", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.F480.filename, plot = p.marker.F480, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.P75NTR <- plotTree(URD_Object.tree, "p75NTR_176Yb", title="P75NTR")
  #p.marker.P75NTR <- p.marker.P75NTR + scale_color_viridis() # Apply viridis color scale
  p.marker.P75NTR <- p.marker.P75NTR + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.P75NTR.filename <- paste0(OUTPUT.BASENAME,"P75NTR", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.P75NTR.filename, plot = p.marker.P75NTR, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.ClCasp3 <- plotTree(URD_Object.tree, "cl_Casp3_173Yb", title="ClCasp3")
  #p.marker.ClCasp3 <- p.marker.ClCasp3 + scale_color_viridis() # Apply viridis color scale
  p.marker.ClCasp3<- p.marker.ClCasp3 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.ClCasp3.filename <- paste0(OUTPUT.BASENAME,"ClCasp3", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.ClCasp3.filename, plot = p.marker.ClCasp3, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.SSEA1 <- plotTree(URD_Object.tree, "SSEA_1_157Gd", title="SSEA-1")
  #p.marker.SSEA1 <- p.marker.SSEA1 + scale_color_viridis() # Apply viridis color scale
  p.marker.SSEA1 <- p.marker.SSEA1 + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.SSEA1.filename <- paste0(OUTPUT.BASENAME,"SSEA1", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.SSEA1.filename, plot = p.marker.SSEA1, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Ncad <- plotTree(URD_Object.tree, "N_Cadherin_161Dy", title="Ncad")
  #p.marker.Ncad <- p.marker.Ncad + scale_color_viridis() # Apply viridis color scale
  p.marker.Ncad <- p.marker.Ncad + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Ncad.filename <- paste0(OUTPUT.BASENAME,"Ncad", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Ncad.filename, plot = p.marker.Ncad, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.Ly6C <- plotTree(URD_Object.tree, "Ly_6C_166Er", title="Ly6C")
  #p.marker.Ly6C <- p.marker.Ly6C + scale_color_viridis() # Apply viridis color scale
  p.marker.Ly6C <- p.marker.Ly6C + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.Ly6C.filename <- paste0(OUTPUT.BASENAME,"Ly6C", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.Ly6C.filename, plot = p.marker.Ly6C, device = OUTPUT.DEVICE) # Save plot output
  
  p.marker.PSANCAM <- plotTree(URD_Object.tree, "PSA_NCAM_139La", title="PSA-NCAM")
  #p.marker.PSANCAM <- p.marker.PSANCAM + scale_color_viridis() # Apply viridis color scale
  p.marker.PSANCAM <- p.marker.PSANCAM + scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1))
  p.marker.PSANCAM.filename <- paste0(OUTPUT.BASENAME,"PSANCAM", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.marker.PSANCAM.filename, plot = p.marker.PSANCAM, device = OUTPUT.DEVICE) # Save plot output
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
  p.segment.filename <- paste0(OUTPUT.BASENAME, "Segment", ".", "svg") # Generate output filename
  ggsave(p.segment.filename, plot = p.segment.final, device = OUTPUT.DEVICE) # Save plot output
}

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
  p.cluster.final <- p.cluster.final + scale_color_manual(values = cluster.colors) # Add custom color scale to stage plot 
  p.cluster.filename <- paste0(OUTPUT.BASENAME, "Cluster", ".", OUTPUT.DEVICE) # Generate output filename
  ggsave(p.cluster.filename, plot = p.cluster.final, device = OUTPUT.DEVICE) # Save plot output
}

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
  print("Finish URD_CyTOF_BuildTree.R")
} else {
  print("Finish URD_CyTOF_BuildTree.R")
}

