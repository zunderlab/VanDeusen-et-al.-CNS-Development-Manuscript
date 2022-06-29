#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start UMAP_R1_Layout_ColoredBy_Select_R2_Clusters_and_Metadata.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(ggfortify)
library(data.table)

print("libraries loaded")

## Input parameters ===============================================================================
INPUT_FOLDER <- getwd()
OUTPUT_FOLDER <- INPUT_FOLDER
CLUSTER_ROUND <- 1
LAYOUT_BASENAME <- "umap_xy_R1_shifted.csv"
CLUSTERS_BASENAME <- "cluster_R2_assigns.csv"
OUTPUT_FILENAME <- "UMAP_R1_layout_R2_clusters_Tips.png"
COLORS_FILENAME <- "colors_R2.csv"  #.csv with column 'color' containing hex codes or colors for R 
FILENUMS_FILENAME <- "filenums.csv"
GRAY.COLOR <- "#bdbdbd"
CLUSTERS.GRAY <- c(1:6,9:14,17,18,23:26,28,31:33,35:37,40:63) #Which metadata conditions to gray. Use string or integer to index (e.g. 1,2,3,5)
  #TipClusters = c(1:6,9:14,17,18,23:26,28,31:33,35:37,40:63)
GRAY.FILES <- c(1:65) #Which metadata conditions to gray. Use string or integer to index
  #Brain = c(5:112), Cortex = c(1:4,7:14,17:23,26:34,37:44,47:55,58:65,68:75,78:85,88:93,96:104,107:112),
  #Dien = c(1:6,9:16,19:25,29:36,39:46,51:57,62:67,72:77,81:87,90:95,99:106,109:112),
  #Mid = c(1:11,15:21,24:31,35:41,45:53,56:63,66:73,76:82,86:91,94:101,105:110),
  #Hind = c(1:8,12:18,22:28,32:38,42:50,54:61,64:71,74:80,83:89,92:98,102:108,111,112)
  #E11 = c(3:112), E15 = c(1:23,35:112), P0 = c(1:65,76:112), P4 = c(1:104), P0-P4 = c(1:65)
OUTPUT_DEVICE <- "png" # "png", "pdf", "jpeg", etc
POINT_SIZE <- 1


print("Input parameters loaded, reading needed files")

## Read needed files ==============================================================================
layout_filename <- LAYOUT_BASENAME
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
clusters_filename <- CLUSTERS_BASENAME
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

filenums <- fread(paste0(INPUT_FOLDER, "/", FILENUMS_FILENAME), stringsAsFactors = F)

print("Needed files read, prepping data to plot")

## Prep dataframe for plotting ====================================================================
plot_df <- as.data.frame(cbind(layout_in, clusters_in))
plot_df <- cbind(plot_df, filenums)
colnames(plot_df) <- c("umap_x", "umap_y", "Cluster", "File")

## Randomize order for plotting ===================================================================
# to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

print("Data ready to plot, plotting")

## Read in custom color palette ===================================================================
colors <- data.frame(read.csv(paste0(COLORS_FILENAME)))
CUSTOM.PALETTE <- as.vector(colors$color)
names(CUSTOM.PALETTE) <- as.character(1:length(CUSTOM.PALETTE))
CUSTOM.PALETTE[CLUSTERS.GRAY] <- GRAY.COLOR

## Save plots colored by each marker ==============================================================
#set output folder
setwd(OUTPUT_FOLDER)
#loop through variables to plot
ggsave(OUTPUT_FILENAME,plot = ggplot() + 
         geom_point(data = dplyr::filter(plot_df,File %in% GRAY.FILES),aes(x=umap_x,y=umap_y),
                    color = GRAY.COLOR, size = POINT_SIZE) + 
         geom_point(data = dplyr::filter(plot_df,!(File %in% GRAY.FILES)),
                    aes(x=umap_x,y=umap_y,color=factor(Cluster)), size = POINT_SIZE) +
         scale_color_manual(values = CUSTOM.PALETTE) +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black")) +
         guides(colour = guide_legend(override.aes = list(shape=15, size=8))),height = 21,width = 21)

print("File outputted")
print("Start UMAP_R1_Layout_ColoredBy_Select_R2_Clusters_and_Metadata.R")
