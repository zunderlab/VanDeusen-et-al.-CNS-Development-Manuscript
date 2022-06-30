#Amy Van Deusen, University of Virginia
#2022
#Visualize 3D UMAP layout with plotly
#See https://plotly.com/r/3d-scatter-plots/#basic-3d-scatter-plot 
#This will allow for .html and other outputs to be interactive for sample metadata

print("Start View_3D_UMAP.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

## Load libraries ==============================================================
library(data.table)
library(plotly)
library(umap)
library(dplyr)
library(RColorBrewer) #Optional: depends on how color plot
library(viridis) #Optional: depends on how color plot
library(htmlwidgets) #Optional: necessary if output HTML file

## Input parameters ============================================================
INPUT.FOLDER <- getwd()
OUTOUT.FOLDER <- INPUT.FOLDER
UMAP.3D.LAYOUT <- "umap_xyz_3D.csv"
EXPRESSION.FILENAME <- "Concat_Transformed_by_Cluster.csv"
CLUSTERS.FILENAME <- "Clusters_Subsampled.csv"
METADATA.FILENAME <- "metadata.csv"
FILENUMS.FILENAME <- "filenums_Subsampled.csv"
SUBSAMPLE.IDS.FILENAME <- "Subsample_ID_by_Cluster.csv" #NULL if no subsample IDs
METADATA.TO.INCLUDE <- c("Age", "Tissue", "SampleType") #list of sample metadata (column IDs) to include in plot dataframe
PLOT.COLOR.FACTOR <- NULL #NULL = Color plot by "Cluster" as default, can also choose metadata (column IDs) 
CLUSTER.COLORS <- "colors_R1.csv" # .csv file or NULL will used default colors
METADATA.COLORS <- NULL # Can add c("color1", "color2", etc.) or palette (e.g. viridis(n))
OUTPUT.BASENAME <- "UMAP_3D_EqualClusters"
OUTPUT.TYPE <- "png" # Choices are "svg" (preferred), "pdf", "png", etc.
SAVE.HTML <- FALSE # TRUE = output .html file, FALSE = do not output file
HTML.FILENAME <- "UMAP_3D_EqualClusters_10K.html"

## Read necessary files ========================================================
expression.filename <- paste0(INPUT.FOLDER, "/", EXPRESSION.FILENAME)
exprs.in <- fread(expression.filename, header=TRUE, sep=",", check.names=FALSE)
layout.filename <- paste0(INPUT.FOLDER, "/", UMAP.3D.LAYOUT)
layout.in <- read.table(layout.filename, header=TRUE, sep=",", check.names=FALSE)
clusters.filename <- paste0(INPUT.FOLDER, "/", CLUSTERS.FILENAME)
clusters.in <- read.table(clusters.filename, header=FALSE, sep=",", check.names=FALSE)
metadata.filename <- paste0(INPUT.FOLDER, "/", METADATA.FILENAME)
metadata.in <- read.table(metadata.filename, header=TRUE, sep=",", check.names=FALSE)
filenums.filename <- paste0(INPUT.FOLDER, "/", FILENUMS.FILENAME)
filenums.in <- read.table(filenums.filename, header=FALSE, sep=",", check.names=FALSE)
if (!is.null(SUBSAMPLE.IDS.FILENAME)) {
  subsamples.filename <- paste0(INPUT.FOLDER, "/", SUBSAMPLE.IDS.FILENAME)
  subsamples.in <- read.table(subsamples.filename, header=FALSE, sep=",", check.names=FALSE)
}

## Prepare plotting dataframe ==================================================
# Combine layout, clusters, and filenums into single dataframe
colnames(filenums.in) <- "File"
plot.df.init <- cbind(filenums.in, layout.in) # Combine file numbers and UMAP layout
colnames(clusters.in) <- "Cluster"
plot.df.init2 <- cbind(plot.df.init, clusters.in) # Add clusters to file numbers and UMAP layout

# Get cell IDs or subsample IDs and add to plotting dataframe
if (!is.null(SUBSAMPLE.IDS.FILENAME)) {
  cell.ids <- subsamples.in 
} else {
  cell.ids <- as.data.frame(1:length(clusters.in$Cluster))
}
colnames(cell.ids) <- "cell.id"
plot.df <- cbind(plot.df.init2, cell.ids)

# Prepare sample metadata (specified in METADATA.TYPE) and add to plotting dataframe
metadata.out <- vector(mode = "list") # Initialize empty list to collect metadata
filenums <- as.vector(filenums.in$File) # Convert file numbers to vector for subsetting below

for (i in METADATA.TO.INCLUDE) {
  metadata_key <- as.character(metadata.in[,i]) # Generate metadata for selected metadata
  names(metadata_key) <- 1:length(metadata.in[,i]) # Generate new index for selected metadata
  metadata.out[[i]] <- metadata_key[filenums] # Add selected metadata to compiled list
}
metadata.final <- as.data.frame(metadata.out) # Convert selected metadata into dataframe
plot.df.final <- cbind(plot.df, metadata.final) # Add metadata to plotting dataframe (quick)
plot.df.exprs <- cbind(plot.df.final, exprs.in) # Generate larger plotting dataframe with expression data


## Generate plotly plot ========================================================
# Colored by clusters (Default)
if (is.null(PLOT.COLOR.FACTOR) || (PLOT.COLOR.FACTOR == "Cluster")) {
  if (!is.null(CLUSTER.COLORS)) {
    colors.filename <- paste0(INPUT.FOLDER, "/", CLUSTER.COLORS) # Find colors.csv
    colors.in <- read.csv(colors.filename) # Read colors.csv
    colors.final <- as.vector(colors.in$color) # Convert colors to vector
    names(colors.final) <- unique(clusters.in)
  } else {
    numclusters <- length(unique(plot.df.final$Cluster)) # Calculate number of clusters to color
    colors.final <- rainbow(numclusters, alpha = 1) # Generate color palette: Can use any RColorBrewer or viridis palette
  }
  p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
               x = ~umap_x, y = ~umap_y, z = ~umap_z, 
               marker = list(size = 0.3),
               alpha = 1, color = ~Cluster, colors = colors.final,
               hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                 '</br> Cluster ID: ', Cluster,
                                                 '</br> Tissue: ', Tissue,
                                                 '</br> Age: ', Age)) %>%
  config(plot_ly(), toImageButtonOptions = list(format = OUTPUT.TYPE, 
                                        filename = OUTPUT.BASENAME, 
                                        height= 1200, width= 1600, scale = 3))

} else if (PLOT.COLOR.FACTOR == "Age") {
  color.palette <- plasma(13, alpha = 1) # Generate color palette: Can use any RColorBrewer or viridis palette
  p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
               x = ~umap_x, y = ~umap_y, z = ~umap_z, 
               marker = list(size = 0.3),
               alpha = 1, color = ~Age, colors = color.palette,
               hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                 '</br> Cluster ID: ', Cluster,
                                                 '</br> Tissue: ', Tissue,
                                                 '</br> Age: ', Age))
  config(plot_ly(), toImageButtonOptions = list(format = OUTPUT.TYPE, 
                                                filename = OUTPUT.BASENAME, 
                                                height= 1200, width= 1600, scale = 3))
} else if (PLOT.COLOR.FACTOR == "Tissue") {
  metadata.col <- match(PLOT.COLOR.FACTOR,names(plot.df.final)) # Figure out column # for metadata
  numfactors <- length(unique(plot.df.final[,metadata.col])) # Calculate number of factors to color
  color.palette <- c("#F242F5", "#4E42F5", "#18C92D","#F01313","#F0F013") # Generate color palette: Can use any RColorBrewer or viridis palette
  p <- plot_ly(plot.df.final, type = 'scatter3d', mode = 'markers', 
               x = ~umap_x, y = ~umap_y, z = ~umap_z, 
               marker = list(size = 0.3),
               alpha = 1, color = ~Tissue, colors = color.palette,
               hoverinfo = 'text', text = ~paste('</br> Cell ID: ', cell.id,
                                                 '</br> Cluster ID: ', Cluster,
                                                 '</br> Tissue: ', Tissue,
                                                 '</br> Age: ', Age))
  config(plot_ly(), toImageButtonOptions = list(format = OUTPUT.TYPE, 
                                                filename = OUTPUT.BASENAME, 
                                                height= 1200, width= 1600, scale = 3))
}  
 

## Visualize 3D plot ===========================================================
# Visualize in R using plotly
p
# Note: Here is where to add/customize what is viewed in the plot 

## Save 3D plot outputs ========================================================
if (isTRUE(SAVE.HTML)) {
  p.bundle <- partial_bundle(p) # Significantly reduces the output file size
  saveWidget(p.bundle, HTML.FILENAME, selfcontained = T)
  # Note: If want single HTML file (e.g. simple output to share), use selfcontained = T
  # If want to generate multiple HTML files for single data library, consider using
  # selfcontained = F to save single "lib" for all .html outputs (way more efficient!)
  # e.g. saveWidget(p.bundle, HTML.FILENAME, selfcontained = F, libdir = "lib")
}

print("Finish View_3D_UMAP.R")
