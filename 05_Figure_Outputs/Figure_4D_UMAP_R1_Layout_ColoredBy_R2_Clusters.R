#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start UMAP_R1_Layout_ColoredBy_R2_Clusters.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(ggfortify)

print("libraries loaded")

## Input parameters ===============================================================================
CLUSTER_ROUND <- 1
LAYOUT_BASENAME <- "umap_xy_R1_shifted.csv"
CLUSTERS_BASENAME <- "cluster_R2_assigns.csv"
OUTPUT_FILENAME <- "UMAP_R1_layout_R2_clusters_shifted.png"
COLORS_FILENAME <- "colors_R2.csv"  #.csv with column 'color' containing hex codes or colors for R 
OUTPUT_DEVICE <- "png" # "png", "pdf", "jpeg", etc
POINT_SIZE <- 0.000001

print("Input parameters loaded, reading needed files")

## Read needed files ==============================================================================
layout_filename <- LAYOUT_BASENAME
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
clusters_filename <- CLUSTERS_BASENAME
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

print("Needed files read, prepping data to plot")

## Prep dataframe for plotting ====================================================================
plot_df <- as.data.frame(cbind(layout_in, clusters_in))
colnames(plot_df) <- c("umap_x", "umap_y", "Cluster")

## Randomize order for plotting ===================================================================
# to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

print("Data ready to plot, plotting")

## Read in custom color palette ===================================================================
colors <- data.frame(read.csv(paste0(COLORS_FILENAME)))
CUSTOM.PALETTE <- as.vector(colors$color)

## Save plots colored by each marker ==============================================================
p <- ggplot(plot_df, aes(x=umap_x, y=umap_y, color=factor(Cluster))) + scale_color_manual(values = CUSTOM.PALETTE) +
    geom_point(size = POINT_SIZE) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    guides(colour = guide_legend(override.aes = list(shape=15, size=8)))
    #^^should work for changing size/shape of legend elements... might have to tweak size per preference
       
output_filename <- OUTPUT_FILENAME
ggsave(output_filename, p,  device = OUTPUT_DEVICE, height = 21,width = 21)

print("Finish UMAP_R1_Layout_ColoredBy_R2_Clusters.R")
