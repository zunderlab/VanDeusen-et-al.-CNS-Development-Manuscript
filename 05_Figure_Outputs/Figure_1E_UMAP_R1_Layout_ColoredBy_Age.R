# Created by Austin Keeler & Ashley Hirt, University of Virginia
# 15 July 2020
# Plot a single UMAP colored by age (or other metadata input)

print("Start UMAP_R1_Layout_ColoredBy_Age.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ggplot2)
library(data.table)
library(viridis)

# ONLY WORKS IF SAME # CELLS IN UMAP.LAYOUT AS CONCAT.TRANSFORMED.

# Input parameters ----
CLUSTER_ROUND <- 1
METADATA_TO_PLOT <- "Age" # make sure this matches your metadata column name exactly

#METADATA_ORDER <- NULL
METADATA_ORDER <- c("E11", "E12", "E13", "E14","E15","E16","E17","E18","P0","P1","P2","P3","P4") # if NULL, then metadata will be sorted alphabetically
#^^ This needs to match with the values in metadata.csv in the METADATA_TO_PLOT column

METADATA_FILENAME <- "metadata.csv"
FILENUMS_FILENAME <- "filenums.csv"
LAYOUT_BASENAME <- "umap_xy_RX_shifted.csv"

POINT_SIZE <- 0.00001
OUTPUT_FILENAME <- paste0("Colored_by_", METADATA_TO_PLOT, "_shifted.png")
LABEL_SIZE <- 5
METADATA_SCALE <- "continuous" # "continuous" or "discrete"
CONTINUOUS_COLORSCALE <- "plasma" # "plasma", "viridis", etc.
OUTPUT_FORMAT <- "png" # "png" or "pdf"
OUTPUT_HEIGHT <- 21
OUTPUT_WIDTH <- 21
# ^^ this is just a default, you can change--supply a full string, eg "My_File_Name.png")

# Read files ----
# read metadata
metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)

# read in filenum labels
filenums <- read.table(FILENUMS_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)[,]

# Read UMAP layout, has two columns, x and y, where each row is a cell
layout_filename <- sub("RX", paste0("R", CLUSTER_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)

# Make metadata key, which will be used to convert metadata and file_nums
metadata_key <- as.character(metadata[,METADATA_TO_PLOT])
names(metadata_key) <- 1:length(metadata[,METADATA_TO_PLOT])

# Get metadata order (important for plotting with discrete and especially continuous scales!)
# If no order is specified by METADATA_ORDER, then calculate alphabetical order.
# Otherwise, use order specified by METADATA_ORDER
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
# Merge with XY coordinates
plot_df <- cbind(layout_in, md_converted_for_plot)
# Assign column names
colnames(plot_df) <- c("umap_x", "umap_y", "md_level")

# Randomize order for plotting, to avoid cells at the end (often
# the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order
plot_df$md_level <- factor(plot_df$md_level,
                           levels=1:length(metadata_order),
                           labels=names(metadata_order))

# save plot
p <- ggplot(data=plot_df, aes(x=umap_x, y=umap_y,
                              group=md_level, color=md_level)) +
    geom_point(size=POINT_SIZE) +
    labs(color=METADATA_TO_PLOT) +
    guides(color = guide_legend(override.aes = list(size=LABEL_SIZE))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #legend.key.size = unit(3,"line")
          axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill = NA))

if (METADATA_SCALE == "discrete") {
    p <- p + scale_color_discrete()
}
if (METADATA_SCALE == "continuous") {
    p <- p + scale_color_viridis(discrete = TRUE,
                                 option = CONTINUOUS_COLORSCALE)
}
ggsave(OUTPUT_FILENAME, p, device = OUTPUT_FORMAT,
       height = OUTPUT_HEIGHT, width = OUTPUT_WIDTH)

print("Finish UMAP_R1_Layout_ColoredBy_Age.R")
