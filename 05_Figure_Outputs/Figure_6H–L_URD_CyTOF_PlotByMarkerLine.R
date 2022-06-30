#Corey Williams, University of Virginia
#20 Aug, 2020
#Make marker expression vs pseudotime plots
#INSTRUCTIONS: In panel file, add column, default name 'URD' and mark for inclusion in plot

print("Start PlotByMarkerLine.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(dplyr)
library(ggfortify)

INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
URD_IN <- "/URD_Final.RData"
load(paste0(INPUT.FOLDER,URD_IN))

#Functions for Zunder lab pipeline
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

#' Panel file input
#'
#' Reads panel file and produces format friendly to scripts in Zunder lab pipeline
#' @param input.folder directory containing metadata file
#' @param panel.filename panel filename, defaults to panel.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.panel <- function(input.folder,panel.filename = "/panel.csv"){
  panel <- read.csv(paste0(input.folder,panel.filename))
  panel$Antigen <- gsub("-", "_", panel$Antigen)
  panel$Antigen <- gsub("\\.", "_", panel$Antigen)
  return(panel)
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

#' Clustering replicate file input
#'
#' Reads clustering replicates for consensus clustering and stability analysis
#' @param input.folder directory containing layout file
#' @param clustering.rep.filename clustering replicate filename, defaults to
#' ClusterStabilityCheck.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.clustering.rep <- function(input.folder,clustering.rep.filename = "/ClusterStabilityCheck.csv")
{
  #Add 1 since clustering done in Python begins indexing at zero
  clustering.rep <- read.csv(paste0(input.folder,clustering.rep.filename),header = FALSE) + 1
  return(clustering.rep)
}

#' Cluster stability file input
#'
#' Reads cluster stability values
#' @param input.folder directory containing layout file
#' @param cluster.stability.filename cluster stability filename, defaults to
#' Final_Cluster_stability.csv as outputted by generation script
#' @importFrom utils read.csv
#' @export
read.cluster.stability <- function(input.folder,
                                   cluster.stability.filename = "/Final_Cluster_Stability.csv"){
  cluster.stability <- as.vector(read.csv(paste0(input.folder,cluster.stability.filename),
                                          header = FALSE)[[1]])
  return(cluster.stability)
}
#' Pull out markers for transformation
#'
#' Finds markers in panel file marked for either clustering or plotting
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.markers <- function(panel){
  transform.markers <- as.character(panel$Metal[panel$Clustering == 1 | panel$Plotting == 1])
  return(transform.markers)
}

#' Get marker names in more legibile format for Concat_Transformed file
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.annotate <- function(panel){
  transform.markers.annotate <- as.character(paste0(panel$Antigen[panel$Clustering == 1 |
                                                                    panel$Plotting == 1],
                                                    "_",panel$Metal[panel$Clustering == 1 |
                                                                      panel$Plotting == 1]))
  return(transform.markers.annotate)
}

#' Get clustering variables in legible format
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.clustering.annotate <- function(panel){
  clustering.markers.annotate <- as.character(paste0(panel$Antigen[panel$Clustering == 1],"_",
                                                     panel$Metal[panel$Clustering == 1]))
  return(clustering.markers.annotate)
}

#' Get panel variables in legible format
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.plotting.annotate <- function(panel){
  plotting.markers.annotate <- as.character(paste0(panel$Antigen[panel$Plotting == 1],"_",
                                                   panel$Metal[panel$Plotting == 1]))
  return(plotting.markers.annotate)
}

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
PANEL_IN <- "panel_URD.csv"
URD_COLNAME <- "URD" #column name in panel for plotting markers
END_SEGMENT <-  10 #tip to generate plot for
OUTPUT_DEVICE <- "png"
PLOT_FILENAME <- "Markers_Pseudotime_TipNumber"

MARKERS.TO.PLOT <- c("Sox2_160Gd","Nestin_144Nd","A2B5_115In", "Prominin_148Nd","PECAM1_143Nd","Ly_6C_166Er",
                     "Pax6_154Sm","p75NTR_176Yb","Tbr1_171Yb","BLBP_172Yb","PDGFRa_147Sm","Tbr2_146Nd",
                     "NeuroD1_156Gd","Ctip2_168Er","TuJ1_89Y","MAP2_164Dy","NeuN_150Nd","GAD65_162Dy","ALDH1A1_165Ho",
                     "GLAST_167Er","GFAP_169Tm","Olig2_113In", "OligoO4_153Eu","Sox10_151Eu")
MARKER.COLORS <- c("#A90F7B","#A1E0FF","#00B0B2","#FFED00","#793F0D","#FF881A","#E19CFF",
                   "#11AB19","#92EA78","#1B9AFF","#7077FF","#DD10A0","#BE26FF","#FFB300",
                   "#0039FF","#36E044","#EC610F","#5EF8FF","#FC1212","#FF85A2","#FF5656",
                   "#FFCC75","#8405E6","#C42323")
#For OPCs (Tip 2, Cluster 10)
MARKERS.TO.PLOT <- c("Sox2_160Gd","p75NTR_176Yb","PDGFRa_147Sm","Olig2_113In", "OligoO4_153Eu","Sox10_151Eu")
MARKER.COLORS <- c("#A90F7B","#11AB19","#7077FF","#FFCC75","#8405E6","#C42323")
#For Glia (Tip 3, Cluster 11)
MARKERS.TO.PLOT <- c("Sox2_160Gd","A2B5_115In","Olig2_113In","BLBP_172Yb","GLAST_167Er","GFAP_169Tm")
MARKER.COLORS <- c("#A90F7B","#00B0B2","#FFCC75","#1B9AFF","#FF85A2","#FF5656")
#For Ctip2 Interneurons (Tip 5, Cluster 19)
MARKERS.TO.PLOT <- c("Sox2_160Gd","Tbr2_146Nd","Ctip2_168Er","TuJ1_89Y","MAP2_164Dy","GAD65_162Dy")
MARKER.COLORS <- c("#A90F7B","#DD10A0","#FFB300","#0039FF","#36E044","#5EF8FF")
#For Tbr1 Neurons (Tip 6, Cluster 21)
MARKERS.TO.PLOT <- c("Pax6_154Sm","Tbr1_171Yb","Tbr2_146Nd","NeuroD1_156Gd","TuJ1_89Y","NeuN_150Nd")
MARKER.COLORS <- c("#E19CFF","#92EA78","#DD10A0","#BE26FF","#0039FF","#EC610F")
#For Ly6C Endothelial Cells (Tip 9, Cluster 32)
MARKERS.TO.PLOT <- c("Sox2_160Gd","Nestin_144Nd","A2B5_115In", "Prominin_148Nd","PECAM1_143Nd","Ly_6C_166Er")
MARKER.COLORS <- c("#A90F7B","#A1E0FF","#00B0B2","#FFED00","#793F0D","#FF881A")
#For P75 Nonneural Cells (Fibroblasts; Tip 10, Cluster 33)
MARKERS.TO.PLOT <- c("Nestin_144Nd","Prominin_148Nd","p75NTR_176Yb","PDGFRa_147Sm","ALDH1A1_165Ho","GLAST_167Er")
MARKER.COLORS <- c("#A1E0FF","#FFED00", "#11AB19","#7077FF","#FC1212","#FF85A2")


## Load URD workspace =============================================================================
input.folder.temp <- INPUT.FOLDER
output.folder.temp <- OUTPUT.FOLDER
INPUT.FOLDER <- input.folder.temp
OUTPUT.FOLDER <- output.folder.temp
object <- URD_Object.tree

START_SEGMENT <-  length(unique(URD_Object.tree@group.ids$segment))

## Get cells on trajectory ========================================================================
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

## Get markers for plotting =======================================================================
# Assign marker colors


if (is.null(MARKERS.TO.PLOT)) {
  urd_markers <- clustering.vars
  names(MARKER.COLORS) <- clustering.vars
} else {
  urd_markers <- MARKERS.TO.PLOT
  names(MARKER.COLORS) <- MARKERS.TO.PLOT
  #MARKER.COLORS <- MARKER.COLORS[names(MARKER.COLORS) %in% urd_markers]
}


## Calculate splines ==============================================================================
#get pseudotime for each cell
pseudotime_trajectory <- object@pseudotime[cells_trajectory,"pseudotime"]
#get stage for each cell
all_stages <- data.frame(Cell_ID = as.character(object@group.ids$init, stringsAsFactors=FALSE))
all_stages$stage <- object@meta$stage
stages_trajectory <- all_stages[match(cells_trajectory, all_stages$Cell_ID),]
stages_trajectory$pseudotime <- pseudotime_trajectory
stages_trajectory$timepoint <- 1

#get expression data for selected markers
exprs_trajectory <- t(as.matrix(object@logupx.data[urd_markers,cells_trajectory]))
#get splines for each marker
splines_trajectory <- lapply(urd_markers,function(this_marker){
  smooth.spline(pseudotime_trajectory,exprs_trajectory[,this_marker],df=5)
})

## Plot ===========================================================================================
#Make output directory
output.dir <- OUTPUT.FOLDER

#Generate df for plotting splines
splines_markers <- lapply(seq_along(urd_markers),function(this_marker){data.frame(splines_trajectory[[this_marker]]$x,splines_trajectory[[this_marker]]$y,rep(urd_markers[this_marker],length(splines_trajectory[[this_marker]]$x)))})
plotting_df_1 <- do.call(rbind,splines_markers)
colnames(plotting_df_1) <- c("pseudotime","expression","marker")
# Get minimum and maximum pseudotime values to scale timepoint bar for x-axis
time.min <- min(stages_trajectory$pseudotime)
time.max <- max(stages_trajectory$pseudotime)


#output all spline curves on one plot
spline_plot <- ggplot(plotting_df_1, aes(x=pseudotime,y=expression,color=marker)) + 
  geom_line(size = 2) +
  scale_x_continuous(limits=c(time.min,time.max), expand=c(0,0)) +
  scale_color_manual(values = MARKER.COLORS) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(), axis.ticks.x = element_blank())

n.stages <- length(unique(URD_Object.tree@meta$stage)) # Uses number of stages from URD object
stage.colors <- plasma(n.stages, alpha = 1, begin = 0, end = 1, direction = 1) # Calculates discrete plasma palette
names(stage.colors) <- unique(URD_Object.tree@meta$stage) # Assign names so stage & color always match


timepoint_plot <- ggplot(stages_trajectory, aes(factor(pseudotime), timepoint, fill = stage, color = stage)) +
  geom_tile() +
  scale_color_manual(values = stage.colors) +
  labs(y = "Age") + labs(x = "Pseudotime") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_blank(),
        axis.title.x = element_text(size = 8),axis.title.y = element_text(size=8,angle=0),
        axis.ticks = element_blank(),legend.position = "none",plot.margin = margin())

output_filename <- paste0(gsub("Number",END_SEGMENT,PLOT_FILENAME),".",OUTPUT_DEVICE)
ggsave(output_filename, egg::ggarrange(spline_plot,timepoint_plot,ncol = 1,heights = c(15,1)), device = "png")

print("Finish PlotByMarkerLine.R")
