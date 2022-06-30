#Corey Williams, University of Virginia
#15 July 2019
#Plot colored by expression of markers

print("Start Heatmap_Age_versus_Marker.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))


library(ggfortify)
library(ggstance)
library(ggpubr)
library(forcats)
library(dplyr)
library(reshape2)

print("Libraries loaded")

## Input parameters ============================================================
INPUT_FOLDER <- getwd()
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
PANEL_FILENAME <- "panel.csv"
METADATA_FILENAME <- "metadata.csv"
FILENUMS_FILENAME <- "filenums.csv"
MARKER_ORDER_FILENAME <- "marker_order.csv" # NULL = hierarchically cluster markers (NOT FUNCTIONAL YET!)
METADATA_TO_PLOT <- "Age"
METADATA_ORDER <- c("E11", "E12", "E13", "E14","E15","E16","E17","E18","P0","P1","P2","P3","P4")
OUTPUT_FILENAME <- "Heatmap_Age_Marker.pdf"

print("Input parameters loaded")

## Read needed files ===========================================================
exprs_mat <- cbind(read.table(EXPRS_MAT_INFILE, header=TRUE, sep=",",
                              check.names=FALSE),
                   read.table(EXPRS_MAT_OTHER, header=TRUE, sep=",",
                              check.names=FALSE))
panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

marker_order <- data.frame(read.csv(paste0(MARKER_ORDER_FILENAME)))

metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)

filenums <- read.table(FILENUMS_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)[,]

print("Needed files successfully loaded")

## Get metadata and information for dataframe =================================
# Make metadata key, which will be used to convert metadata and file_nums
metadata_key <- as.character(metadata[,METADATA_TO_PLOT])
names(metadata_key) <- 1:length(metadata[,METADATA_TO_PLOT])

# Get metadata order (important for plotting with discrete and especially continuous scales!)
# If no order is specified by METADATA_ORDER, then calculate alphabetical order.
if (is.null(METADATA_ORDER)) {
  unique_sorted <- sort(unique(as.character(metadata[,METADATA_TO_PLOT])))
  metadata_order <- 1:length(unique_sorted)
  names(metadata_order) <- sort(unique(as.character(metadata[,METADATA_TO_PLOT])))
} else {
  metadata_order <- 1:length(METADATA_ORDER)
  names(metadata_order) <- METADATA_ORDER
}

# Make filenum key, which will be used to convert File_Num to plotting level
filenum_key <- metadata_order[metadata_key]
names(filenum_key) <- 1:length(metadata[,METADATA_TO_PLOT])

# Convert filenames for each cell event to metadata values
md_converted_for_plot <- filenum_key[filenums]

print("Metadata loaded")

## Prep dataframe with metadata and expression values ==========================
# Read in custom marker order 
if (!is.null(MARKER_ORDER_FILENAME)) {
  plot_vars_ordered <- rev(as.vector(marker_order$order)) # Marker order
} else {
  ## ADD MARKER NAMES FROM OTHER SOURCE (FOR HIERARCHICAL CLUSTERING)
  print("Hierarchical clustering of markers not functional. Please add marker_order.csv file.")
}

print("Custom marker order loaded")

# Combine data into single plotting dataframe
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars_ordered],md_converted_for_plot))
names(plot_df)[names(plot_df) == 'md_converted_for_plot'] <- 'Age'
print("Dataframe prepared with standard order")

## Generate heatmap dataframe with mean expression values ======================
# Initiate plotting list
plot_list <- vector(mode = "list", length = length(plot_vars_ordered))

# Get mean values for each age for each marker
for (i in 1:length(METADATA_ORDER)) {
  plot_df_meta <- filter(plot_df, Age == i)
  plot_list[[i]] <- colMeans(plot_df_meta)
}

plot_list_df <- t(do.call("cbind", plot_list))
pheatmap_df <- as.matrix(plot_list_df[,1:40])
pheatmap_df <- t(pheatmap_df)
colnames(pheatmap_df) <- METADATA_ORDER

low_colors <- viridis(500, begin=0, end=0.2)
mid_colors <- viridis(500, begin=0.2, end=0.95)
high_colors <- viridis(300, begin = 0.95, end = 1)
custom_colors <- c(low_colors,mid_colors,high_colors)

plot <- pheatmap(mat = pheatmap_df, filename = OUTPUT_FILENAME, 
                 cluster_cols = FALSE, scale = "row", 
                 color = custom_colors, border_color = NA)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot, OUTPUT_FILENAME)

# Transform dataframe into format for ggplot heatmap
plot_df_final <- reshape::melt(plot_list_df, id.vars = 'Age')
colnames(plot_df_final) <- c("Age","Marker","Mean_Expression")

print("Custom order applied to dataframe, ready for plotting")

## Generate heatmap ============================================================
if (!is.null(MARKER_ORDER_FILENAME)) {
  
  ggplot(plot_df_final, aes(x = Age, y = Marker, fill = Mean_Expression)) +
    geom_tile() +
    coord_fixed()
  
} else {
  ## ADD HEATMAP WITH HIERARCHICAL CLUSTERING HERE
  print("Hierarchical clustering of markers not functional. Please add marker_order.csv file.")
}

ggsave("Heatmap_Age_Marker.pdf")

print("Finish Heatmap_Age_versus_Marker.R")
