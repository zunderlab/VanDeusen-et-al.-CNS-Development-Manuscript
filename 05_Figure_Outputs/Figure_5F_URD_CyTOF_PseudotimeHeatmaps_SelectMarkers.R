#Corey Williams, University of Virginia
#2019
#Make marker expression vs pseudotime plot for each tip
#Warning: Sometimes timepoint heatmap will be shifted 1-2 pixels out of alignment with expression

print("Start PseudotimeHeatmaps_SelectMarkers.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(dplyr)
library(reshape2)
library(ggfortify)
library(egg)
library(raster)
library(Matrix)
library(URD)

print("Libraries loaded")

## Load R Data ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
URD.RDATA.IN <- "/URD_Final.RData" #Change filenames as relevant, section should run without issue
load(paste0(INPUT.FOLDER, URD.RDATA.IN)) #.Rdata from URD run must have completed through Biased Walks (or beyond) to be used in this script

print("URD Data loaded")

## Input parameters ===============================================================================
# Initial parameters (must restate in case overwritten by URD.RData!!!)
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
PLOT_BASENAME <- "Pseudotime_Heatmap_Tip" #Must include "Tip" for automated output naming
OUTPUT_DEVICE <- "pdf"
PSEUDOTIME_BINS <- 3000 #number of columns in heatmap
NORMALIZE_GLOBALLY <- TRUE
NORMALIZE_LOCALLY <- FALSE
PLOT.PARAMS.CSV <- "Dot_Plot_Parameters.csv" # Only needed if ADD.CLUSTER.BAR = TRUE; NULL = default colors, otherwise .csv file with column "color"
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 5
LEGEND_WIDTH_FACTOR <- 50
TEXT_SIZE <- 10
MARKER_ORDER <- c("Pax6_154Sm","Sox2_160Gd","Nestin_144Nd","N_Cadherin_161Dy","CD24_159Tb",
                  "Tbr2_146Nd","NeuroD1_156Gd","Tbr1_171Yb","Cux1_170Er","Ctip2_168Er",
                  "Dcx_163Dy","TuJ1_89Y","MAP2_164Dy","NeuN_150Nd","GAD65_162Dy",
                  "A2B5_115In", "BLBP_172Yb","GLAST_167Er","GFAP_169Tm","Olig2_113In",
                  "PDGFRa_147Sm", "Sox10_151Eu","OligoO4_153Eu","TrkB_174Yb","VCAM1_142Nd",
                  "Ki67_152Sm","cl_Casp3_173Yb")
#NULL = heirarchical clustering or markers; to reorder, use a vector like c(1,3,2)  urd_markers (line 51). 

print("User inputs loaded")

# Load necessary files
plot.params.filename <- paste0(INPUT.FOLDER,"/",PLOT.PARAMS.CSV)
plot.params <- read.csv(plot.params.filename, header = TRUE)

# Set color scales ===============================================================================================
#For stages (sets color to plasma by default)
n.stages <- length(unique(URD_Object.tree@meta$stage)) # Uses number of stages from URD object
stage.colors <- plasma(n.stages, alpha = 1, begin = 0, end = 1, direction = 1) # Calculates discrete plasma palette
names(stage.colors) <- unique(URD_Object.tree@meta$stage) # Assign names so stage & color always match

# For 
n.segments <- length(URD_Object.tree@tree$segments)
segment.colors <- rainbow(n.segments, alpha = 1)
names(segment.colors) <- URD_Object.tree@tree$segments

print("Color scales defined")

## Establish initial parameters from URD object ==============================================
object <- URD_Object.tree
urd_markers <- clustering.vars
START_SEGMENT <- length(unique(URD_Object.tree@group.ids$segment))
END_SEGMENT <- 1 #Initializes generation of plots at Tip #1
if (is.null(START_SEGMENT) | is.null(END_SEGMENT)){
  stop("Specify START_SEGMENT and END_SEGMENT")
}

## Get cells on trajectory ========================================================================
for (i in tip.names) {
  END_SEGMENT <- match(i,tip.names)
  #get order of segments from end to start
  segments_trajectory <- END_SEGMENT
  this_segment <- as.character(segments_trajectory)
  while (this_segment %in% object@tree$segment.joins$child){
    this_segment <- dplyr::filter(object@tree$segment.joins,child == this_segment)$parent
    segments_trajectory <- c(segments_trajectory,as.integer(this_segment))
  }
  #Get cells in each segment
  cells_trajectory <- unlist(lapply(segments_trajectory,function(this_segment){
    object@tree$cells.in.segment[[toString(this_segment)]]
  }))
  
  ## Make plotting data frame =======================================================================
  #get pseudotime for each cell
  pseudotime_trajectory <- object@pseudotime[cells_trajectory,"pseudotime"]
  #get pseudotime quantiles
  pseudotime_quantiles <- quantile(pseudotime_trajectory,seq(0,1,length.out=PSEUDOTIME_BINS))
  #get pseudotime bin info
  pseudotime_bins <- cut(pseudotime_trajectory,unique(pseudotime_quantiles[1:(PSEUDOTIME_BINS)]),
                         include.lowest=TRUE)
  #get cells for each bin
  cell_ids_in_bin <- lapply(sort(unique(pseudotime_bins)),
                            function(this_bin){which(pseudotime_bins==this_bin)})
  #get expression info for each bin
  exprs_bin <- lapply(urd_markers,function(this_marker){
    data.frame(sort(unique(pseudotime_bins)),
               sapply(cell_ids_in_bin,function(this_bin){
                 mean(object@logupx.data[this_marker,cells_trajectory[this_bin]])}),
               rep(this_marker,times = length(unique(pseudotime_bins))))
  })
  heatmap_df <- do.call(rbind,exprs_bin)
  colnames(heatmap_df) <- c("pseudotime","expression","marker")
  #normalize expression of markers (need for hierarchical clustering or for plotting)
  if (NORMALIZE_GLOBALLY == TRUE) {
    for(this_marker in urd_markers){
      these_bins <- which(heatmap_df$marker == this_marker)
      heatmap_df$expression[these_bins] <- heatmap_df$expression[these_bins]/
        max(object@logupx.data[this_marker,])
    }  
  } else if (NORMALIZE_LOCALLY == TRUE){
    for(this_marker in urd_markers){
      these_bins <- which(heatmap_df$marker == this_marker)
      heatmap_df$expression[these_bins] <- heatmap_df$expression[these_bins]/
        max(heatmap_df$expression[these_bins])
    }
  }
  #normalize markers for plotting
  if (NORMALIZE_LOCALLY == TRUE){
    heatmap_df <- heatmap_normalized
  }
  #do hierarchical clustering whenever marker order is not manually specified
  if (is.null(MARKER_ORDER) == TRUE){
    #do hierarchical clustering to get order of markers
    hclust_df <- dcast(heatmap_df,marker ~ pseudotime,value.var = "expression")
    hclust_dist <- dist(hclust_df,method = "euclidean")
    hclust_out <- hclust(hclust_dist,method = "ward.D")
    #reorder levels of plotting df based on hclust
    heatmap_df$marker <- factor(heatmap_df$marker,levels = urd_markers[hclust_out$order])
  } else {
    #reorder levels of plotting df based on hclust
    heatmap_df$marker <- factor(heatmap_df$marker,ordered=TRUE, levels = rev(MARKER_ORDER))
  }
  
  # Remove NA values from heatmap dataframe (necessary if select only certain markers)
  heatmap_df <- na.omit(heatmap_df)
  
  #get stage for each cell
  all_stages <- data.frame(Cell_ID = as.character(object@group.ids$init, stringsAsFactors=FALSE))
  all_stages$stage <- object@meta$stage
  stages_trajectory <- all_stages[match(cells_trajectory, all_stages$Cell_ID),]
  stages_trajectory$pseudotime <- pseudotime_trajectory
  stages_trajectory$timepoint <- 1

  ## Generate plots ===========================================================================================
  # Make expression heatmap
  exprs_plot <- ggplot(heatmap_df) +
    geom_tile(aes(x = pseudotime,y = marker,fill = expression)) +
    scale_fill_viridis_c() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = TEXT_SIZE,color = "black",hjust = 1),
          axis.title.y = element_blank(),axis.title.x = element_blank(),
          axis.ticks = element_blank(),legend.title = element_text(size = TEXT_SIZE),
          legend.key.width = unit(PLOT_WIDTH/LEGEND_WIDTH_FACTOR,"in"),plot.margin = margin())
  
  # Make timepoint heatmap bar
  timepoint_plot <- ggplot(stages_trajectory, aes(factor(pseudotime), timepoint, fill = stage, color = stage)) +
    geom_tile() +
    scale_color_manual(values = stage.colors) +
    labs(y = "Age") + labs(x = "Pseudotime") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.text = element_blank(),
          axis.title.x = element_text(size = 8),axis.title.y = element_text(size=8,angle=0),
          axis.ticks = element_blank(),legend.position = "none",plot.margin = margin())
  
  # Output heatmap plot with stage (timepoint) bar
  PLOT_FILENAME <- sub("Tip", i, PLOT_BASENAME)
  PLOT_FILENAME <- gsub(" ", "_", PLOT_FILENAME)
  PLOT_FILENAME <- paste0(PLOT_FILENAME, ".", OUTPUT_DEVICE)
  ggsave(PLOT_FILENAME, egg::ggarrange(exprs_plot,timepoint_plot,ncol = 1,heights = c(length(urd_markers)-1,1)), device = OUTPUT_DEVICE)
  print(paste0("Heatmap plot output for ", i))
  
}

print("Finish PseudotimeHeatmaps_SelectMarkers.R")