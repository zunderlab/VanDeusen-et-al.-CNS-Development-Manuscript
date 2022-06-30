#Amy Van Deusen and Eli Zunder
#2022
# Generate heatmap of lag correlation values

print("Start CorrelationValuesHeatmap.R")

rm(list = ls(all = TRUE))
.libPaths("/R/4.1.1")

library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dendextend)
library(ggdendro)
library(gridExtra)

# Input parameters
INPUT_FOLDER <- getwd()
OUTPUT_FOLDER <- INPUT_FOLDER
LAG_VALUES_FILENAME <- "protein_vs_rna_cross_corr_lags_ab_names.csv"
OUTPUT_BASENAME <- "CorrelationValues_Heatmap"
OUTPUT_DEVICE <- "pdf"

TOP_MARGIN <- 15 # for heatmap
BOTTOM_MARGIN <- 5 # for heatmap
RIGHT_MARGIN <- 75 # for dendrogram

# Read in lag values
lag_filename <- paste0(INPUT_FOLDER,"/",LAG_VALUES_FILENAME)
lags_in <- read.csv(lag_filename)

# Combine rows by marker, to be used for clustering/ordering
markers_all <- lags_in[,"Marker"]
markers_unique <- unique(lags_in[,"Marker"])
vals_only <- lags_in[,4:ncol(lags_in)]
comb_vals <- c()
for(m in markers_unique) {
  comb <- colMeans(vals_only[which(markers_all==m),])
  comb_vals <- rbind(comb_vals, comb)
}
rownames(comb_vals) <- markers_unique

# Cluster based on combined rows (pooled by marker)
hc <- hclust(dist(comb_vals))
hd <- as.dendrogram(hc)

# References for leave swapping.  Run with this line commented out after final configuration is determined.
labels(hd) <- paste0(1:34 , "_",labels(hd)) 


# Set up function for outputting PDFs
heat_with_dend <- function(basename, in_dend, in_data, file_out) {
  # Reorder heatmap dataset by dendrogram
  reordered_vals <- c()
  for (m in markers_unique[order.dendrogram(in_dend)]) {
    m_vals <- in_data[which(markers_all==m),]
    reordered_vals <- rbind(reordered_vals,m_vals)
  }
  row.names(reordered_vals) <- NULL
  
  # Reshape data for plotting
  plot_df <- melt(reordered_vals)
  plot_df$Sample <- factor(plot_df$Sample,levels=rev(unique(plot_df$Sample)))
  
  # Plot heatmap -- note that margins above should be adjusted to line up with dendrogram!
  plot_heat <- ggplot(plot_df, aes(variable, y = Sample, fill=value)) +
    geom_tile() + theme(axis.text = element_text(size = 3), axis.title.y = element_blank()) +
    scale_fill_gradientn(colours=c("#0014FF","#FFFFFF","#FF0000"),values=c(0,0.5,1)) +
    theme(plot.margin = unit(c(TOP_MARGIN,0,BOTTOM_MARGIN,0), "pt")) +
    scale_y_discrete(labels=rev(rep(c("telencephalon", "diencephalon", "mesencephalon", "rhombencephalon"),length(markers_unique))))
  
  # Need this to have it start at the top
  dend_data <- dendro_data(rev(in_dend))
  
  # Plot dendrogram -- not that margin above should be adjusted to fit text
  plot_dend <- ggplot() + theme_void() +
    geom_segment(data = segment(dend_data), aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) +
    geom_text(data = label(dend_data), aes_string(x = "x", y = "y", label = "label"), hjust = 0, angle = 0) +
    #scale_y_continuous(expand=c(0.2, 0)) +
    #scale_x_continuous(expand=c(0.1, 1)) +
    coord_flip(clip="off") +
    scale_y_reverse() +
    scale_x_continuous(limits=c(0,34)) +
    theme(plot.margin = unit(c(0,RIGHT_MARGIN,0,0), "pt"))
  
  # Create ggplot object with heatmap and dendrogram side-by-side
  require(gridExtra)
  plot_grid <- grid.arrange(plot_dend, plot_heat, ncol=2, widths = c(1,2))
  
  # Output to R Plot or PDF
  if(file_out==TRUE) {
    outfile <- paste0(OUTPUT_BASENAME,"_",basename,".",OUTPUT_DEVICE)
    ggsave(outfile, plot_grid, device = OUTPUT_DEVICE)
  } else {
    plot_grid
  }
}

heat_with_dend("plot_original", hd, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd,hc$order[16:25],hc$order[26:34])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[11:15],hc$order[9:10])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[1:8],hc$order[c(11:15,9:10,26:34)])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[19:25],hc$order[16:18])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[1],hc$order[2:5])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[2],hc$order[3])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

hd_mod <- flip_leaves(hd_mod,hc$order[3:2],hc$order[4:5])
heat_with_dend("plot_mod", hd_mod, lags_in, file_out=FALSE)

heat_with_dend("plot_mod", hd_mod, lags_in, file_out=TRUE)

print("Start CorrelationValuesHeatmap.R")
