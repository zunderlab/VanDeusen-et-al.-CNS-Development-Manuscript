#Corey Williams, University of Virginia
#15 Jul, 2019
#Generate UMAP colored by clusters

print("Start UMAP_R2_Layout_ColoredBy_R2_Clusters.R")

rm(list = ls(all = TRUE))
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ggfortify)

print("libraries loaded")

## Input parameters ===============================================================================
CLUSTER_ROUND <- 2
LAYOUT_BASENAME <- "umap_xy_RX.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
COLORS_FILENAME <- "colors_R2_groups.csv"  #.csv with column 'color' containing hex codes or colors for R 
OUPUT_BASENAME <- sub("RX_assigns.csv", paste0("R", CLUSTER_ROUND, "_xy_Custom_Group_#."), CLUSTERS_BASENAME)
OUTPUT_DEVICE <- "png" # "png", "pdf", "jpeg", etc
POINT_SIZE <- 0.000001

print("Input parameters loaded, reading needed files")

## Read needed files ==============================================================================
layout_filename <- sub("RX", paste0("R", CLUSTER_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

print("Needed files read, prepping data to plot")

## Prep dataframe for plotting ====================================================================
plot_df <- as.data.frame(cbind(layout_in, clusters_in))
colnames(plot_df) <- c("Group", "umap_x", "umap_y", "Cluster")

## Randomize order for plotting ===================================================================
# to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

print("Data ready to plot, plotting")

## Read in custom color palette ===================================================================
colors <- data.frame(read.csv(paste0(COLORS_FILENAME)))

print("Custom color palette loaded")

## Save plots colored by cluster==============================================================
groups <- sort(unique(plot_df$Group))

for (g in groups) {
  GROUP.PALETTE <- filter(colors, group == g)
  GROUP.PALETTE <- as.vector(GROUP.PALETTE$color)
  plot_df_group <- filter(plot_df, Group == g)
  p <- ggplot(plot_df_group, aes(x=umap_x, y=umap_y, color=factor(Cluster))) + geom_point(size = POINT_SIZE) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    guides(colour = guide_legend(override.aes = list(shape=15, size=8)))
  #^^should work for changing size/shape of legend elements... might have to tweak size per preference
  
  print("Custom colored plot generated")
  
  
  p <- p + scale_color_manual(values = GROUP.PALETTE)
  output_filename <- paste0(OUPUT_BASENAME, OUTPUT_DEVICE)
  output_filename <- gsub("#", g, output_filename)
  #output_filename <- "Panel_2E_Test.png"
  ggsave(output_filename, p,  device = OUTPUT_DEVICE, height = 21,width = 21)
  
  print("File outputted")
  
  }

print("Finish UMAP_R2_Layout_ColoredBy_R2_Clusters.R")
