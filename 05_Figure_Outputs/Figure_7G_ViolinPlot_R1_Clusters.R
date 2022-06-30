#Corey Williams, University of Virginia
#Corey Williams, University of Virginia
#15 Jul, 2019
#Generate violin plot colored by cluster


print("Start Violin_Plots_R1_Clusters.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(ggfortify)
library(ggstance)
library(ggpubr)
library(forcats)
library(dplyr)

print("Libraries loaded")

## Input parameters
INPUT_FOLDER <- getwd()
CLUSTER_ROUND <- 1
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
PANEL_FILENAME <- "panel.csv" # Only needed if want to pull in markers old way
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
MARKER_ORDER_FILENAME <- "marker_order_select.csv"
PLOT_PARAMETERS_FILENAME <- "Dot_Plot_Parameters.csv" #has to be individually crafted with labels and desired colors
OUTPUT_FILENAME <- "ViolinPlots_Select_R1_Clusters.pdf"

VIOLIN_HEIGHT_FACTOR <- 5

print("Input parameters loaded")

## Read needed files
exprs_mat <- cbind(read.table(EXPRS_MAT_INFILE, header=TRUE, sep=",",
                              check.names=FALSE),
                   read.table(EXPRS_MAT_OTHER, header=TRUE, sep=",",
                              check.names=FALSE))
panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

plot_parameters <- fread(paste0(INPUT_FOLDER, "/", PLOT_PARAMETERS_FILENAME), stringsAsFactors = F)

marker_order <- data.frame(read.csv(paste0(MARKER_ORDER_FILENAME)))

print("Needed files successfully loaded")

## Read in custom marker order, cluster order, and color palette
plot_vars_ordered <- as.vector(marker_order$order) # Marker order
VIOLIN.ORDER <- as.vector(plot_parameters$Cluster) # Cluster order
CUSTOM.PALETTE <- as.vector(plot_parameters$Color) # Custom colors (in order)

print("Custom color palette loaded")

## Prep dataframe for plotting
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars_ordered],clusters_in))
colnames(plot_df)[ncol(plot_df)] <- "Cluster"

print("Dataframe prepared with standard order")

#Apply custom cluster order to dataframe
plot_df_ordered <- plot_df %>%
    mutate(Cluster=factor(Cluster,levels=VIOLIN.ORDER))

print("Custom order applied to dataframe, ready for plotting")

#Make list of violin plots by cluster
plist = sapply(plot_vars_ordered, function(marker_plot) {
  if (marker_plot == plot_vars_ordered[1]){
    ggplot(plot_df_ordered, aes(x = plot_df_ordered[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            legend.position = "none",axis.text.x = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
  else {
    ggplot(plot_df_ordered, aes(x = plot_df_ordered[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            legend.position = "none",axis.text = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
}, simplify=FALSE)

print("List of violin plots by cluster generated")

#save Violin plots
ggsave(OUTPUT_FILENAME,
       annotate_figure(ggarrange(plotlist = plist,ncol=length(plist)),
                       left = text_grob("Cluster",rot=90)),
       height=max(clusters_in)/VIOLIN_HEIGHT_FACTOR,width=length(plot_vars_ordered))

print("Finish Violin_Plots_R1_Clusters.R")
