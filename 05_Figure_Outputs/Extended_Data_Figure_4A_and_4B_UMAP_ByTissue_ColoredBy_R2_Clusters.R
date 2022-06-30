#Corey Williams, University of Virginia
#15 Jul, 2019
#Plot colored by expression of markers

print("Start UMAP_ByTissue_ColoredBy_R2_Clusters.R")

rm(list = ls(all = TRUE))
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ggfortify)

print("libraries loaded")

## Input parameters ===============================================================================
METADATA_TO_PLOT <- "Tissue" # make sure this matches your metadata column name exactly

#METADATA_ORDER <- NULL
METADATA_ORDER <- c("Brain","Cortex","Diencephalon","Midbrain","Hindbrain") # if NULL, then metadata will be sorted alphabetically
#^^ This needs to match with the values in metadata.csv in the METADATA_TO_PLOT column

METADATA_FILENAME <- "metadata.csv"
FILENUMS_FILENAME <- "filenums.csv"

CLUSTER_ROUND <- 1
LAYOUT_BASENAME <- "umap_xy_RX.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
OUPUT_BASENAME <- sub("RX_assigns.csv",
                      paste0("R", CLUSTER_ROUND,
                             "_xy_by_",
                             METADATA_TO_PLOT,
                             "_MX."),
                      CLUSTERS_BASENAME)
OUTPUT_DEVICE <- "png" # "png", "pdf", "jpeg", etc
POINT_SIZE <- 0.001
PLOT_WTH <- 21 #1000 for png, 10 for pdf
PLOT_HGT <- 21 #1000 for png, 10 for pdf
BASE_PLOT_WTH <- 2000 #1000 for png, 10 for pdf
BASE_PLOT_HGT <- 2000 #1000 for png, 10 for pdf

print("Input parameters loaded, reading needed files")

## Read needed files ==============================================================================
# read metadata
metadata <- read.table(METADATA_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)
# read in filenum labels
filenums <- read.table(FILENUMS_FILENAME, header=TRUE, sep=",",
                       check.names=FALSE)[,]

layout_filename <- sub("RX", paste0("R", CLUSTER_ROUND), LAYOUT_BASENAME)
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

print("Needed files read, prepping data to plot")

## Prep dataframe for plotting ====================================================================
full_df <- as.data.frame(cbind(layout_in, clusters_in))
colnames(full_df) <- c("umap_x", "umap_y", "Cluster")

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
full_df <- cbind(full_df, md_converted_for_plot)

# Randomize order for plotting, to avoid cells at the end (often
# the last file concatenated) going on top of everything else
set.seed(42)
full_df <- full_df[sample(nrow(full_df)),]
row.names(full_df) <- NULL #remove old indices, so everything has a new order

print("Data ready to plot, plotting")

## Save plots colored by each marker ==============================================================
for (m in metadata_order) {
    meta_df <- full_df[which(full_df$md_converted_for_plot==m),]    
    p <- ggplot() +
        geom_point(data = full_df,
                   aes(x=umap_x, y=umap_y),
                   size=POINT_SIZE, color="grey") +
        geom_point(data = meta_df,
                   aes(x=umap_x, y=umap_y, color=factor(Cluster)),
                   size=POINT_SIZE) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        guides(colour = guide_legend(override.aes = list(shape=15, size=8)))
    #^^should work for changing size/shape of legend elements... might have to tweak size per preference

    output_filename <- paste0(sub("MX", names(metadata_order)[m],
                                  OUPUT_BASENAME), OUTPUT_DEVICE)
    ggsave(output_filename, p, device=OUTPUT_DEVICE,
           height=PLOT_HGT, width=PLOT_WTH)
    print(paste0(output_filename, " outputted"))
}

# For each group, pull out clusters and XY coordinates, then plot

# Try and calculate the best point size --may want to try and find
# a better function that gives good plots for all datasets.
point.size <- function(n) {
    #pt.size <- 0.5+10/(nrow(um_R1)^0.5) # min=0.5, max=100
    #pt.size <- 1 + 20*(0.999^(n*10))
    if (n <= 70) {
        pt.size <- 25
    } else if (n <= 100000000) {
        pt.size <- 20/log10(0.1*n) - 2.5  # meets at 23.66589
    } else {
        pt.size <- 2.5/log10(0.1*n) # meets at 0.3571429
    }
    return(pt.size)
}
point_size <- point.size(nrow(full_df))
output_filename <- paste0(sub("MX", "All", OUPUT_BASENAME), OUTPUT_DEVICE)
pr <- ceiling(sqrt(length(metadata_order)))
adj.wth <- BASE_PLOT_WTH*pr/2
adj.hgt <- BASE_PLOT_HGT*pr/2
# Write out the plot as PNG or PDF
if (OUTPUT_DEVICE == "png") {
    png(output_filename, width=adj.wth, height=adj.hgt)
} else if (OUTPUT_DEVICE == "pdf") {
    pdf(output_filename, width=adj.wth, height=adj.hgt)
}
par(mfrow=c(pr,pr))
#par(mar=c(0,0,0,0))

for (m in metadata_order) {
    meta_df <- full_df[which(full_df$md_converted_for_plot==m),]
    
    plot(x=full_df[,"umap_x"], y=full_df[,"umap_y"], pch=".",
         cex=point_size, col="grey", main=names(metadata_order)[m],
         axes=FALSE, ylab="", xlab="", cex.main=5)
    points(x=meta_df[,"umap_x"], y=meta_df[,"umap_y"], pch=".",
           cex=point_size, col="red")
    box()
}
dev.off()

print("Finish UMAP_ByTissue_ColoredBy_R2_Clusters.R")
