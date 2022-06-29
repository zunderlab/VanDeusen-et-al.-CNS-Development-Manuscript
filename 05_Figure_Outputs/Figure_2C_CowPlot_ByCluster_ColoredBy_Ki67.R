#Kristen Fread, University of Virginia
#20 Oct 2021
#Generate cowplot with cluster abundances and color according to expression of one marker

print("Start CowPlot_ByCluster_ColoredBy_Marker.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ggfortify)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(viridis)
library(cowplot)

## Input parameters
INPUT.FOLDER <- getwd()
FILES.FILENAME <- "filenums.csv" #to get this file from pipeline 1, run Save_file_nums_pipeline1.R to separate the filesnumbs from concattransformed
METADATA.FILENAME <- "metadata.csv"
CLUSTERS.FILENAME <- "cluster_R2_assigns.csv"
OUTPUT.FILENAME.Cowplot <- "CowPlot_R2_ByTissue_Ki67_4_2.pdf"
CLUSTER.NAMES.ORDER.COLOR <- "Dot_Plot_Parameters_4_2.csv" #has to be individually crafted with labels and desired colors
CONCAT.TRANSFORMED.FILENAME <- "/Concat_Transformed.csv"
OUTPUT.DEVICE <- "pdf" # Choose "png" or "pdf"
METADATA.TYPE <- "Tissue" # Metadata to sort cow plots by on Y axis
METADATA.ORDER <- c("Hindbrain", "Midbrain", "Diencephalon", "Cortex", "Brain") # Reverse order want metadata on y axis
MARKER.TO.PLOT <- "Ki67_152Sm"
SAMPLE.GROUP.ID <- "SampleType" # Metadata value that groups/sorts metadata on x and y axis
SAMPLE.GROUP.ORDER <- c("E11B", "E12B", "E13C", "E13D", "E13M", "E13H", "E14C", "E14D", "E14M", "E14H", "E15C", "E15D", "E15M", "E15H", "E16C", "E16D", "E16M", "E16H", "E17C", "E17D", "E17M", "E17H", "E18C", "E18D", "E18M", "E18H", "P0C", "P0D", "P0M", "P0H", "P1C", "P1D", "P1M", "P1H", "P2C", "P2D", "P2M", "P2H", "P3C", "P3D", "P3M", "P3H", "P4C", "P4D", "P4M", "P4H")


## Necessary import functions ====================================================
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
  md$file_name <- as.character(md$File_Name)
  return(md)
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

## Read in necessary files =================================================
# Read in file numbers
fil_nums <- fread(paste0(INPUT.FOLDER, "/", FILES.FILENAME), stringsAsFactors = F)
colnames(fil_nums) <- "File"
# Read in clusters
clusters.in <- read.csv(paste0(INPUT.FOLDER,"/",CLUSTERS.FILENAME))
colnames(clusters.in) <- "Cluster"
# read metadata
md <- fread(paste0(INPUT.FOLDER, "/", METADATA.FILENAME), stringsAsFactors = F)
# assign file numbering
md <- cbind(md, "File" = 1:nrow(md)) 
#read this in for dot plot ordering
ord_name_color <- fread(paste0(INPUT.FOLDER, "/", CLUSTER.NAMES.ORDER.COLOR), stringsAsFactors = F)
# red in expression data
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCAT.TRANSFORMED.FILENAME), stringsAsFactors = F)

## Prep cluster abundances for plotting data frame ===============================================
plotting_df <- as.data.frame(cbind(fil_nums,clusters.in))
colnames(plotting_df) <- c("File_Num", "Cluster")

#join metadata to the dataframe with file and cluster number
cluster_file_metadata <- plotting_df %>% 
  left_join(md, by="File_Num")

#calculate number of cells per age
totalcellstimepoint <- cluster_file_metadata %>% 
  group_by(Tissue, Age) %>%
  dplyr::summarise(num_cells=n())

total_cells <- sum(totalcellstimepoint$num_cells)

#group matrix by cluster and file then calculate how much of each file in the cluster
clusterz <- plotting_df %>%
  group_by(File_Num, Cluster) %>% 
  dplyr::summarize(num_cells=n())

#calculate the percent breakdown of each cluster contributing to a file/timepoint (on a per file basis, find the makeup of each file)
percent_similar <- clusterz %>%
  group_by(File_Num) %>%
  dplyr::summarize(total_cells_per_file = sum(num_cells)) %>%
  inner_join(clusterz) %>%
  mutate(fraction_cells_in_file = num_cells/total_cells_per_file) %>% 
  mutate(percent_cells_in_file = fraction_cells_in_file*100)

#calculate the normalized fraction of each timepoint within each cluster
cluster_plot_values <- percent_similar %>% 
  group_by(Cluster) %>% 
  dplyr::summarize(total_percents_in_cluster = sum(percent_cells_in_file)) %>% 
  inner_join(percent_similar) %>% 
  mutate(normed_fraction_file_in_cluster = percent_cells_in_file/total_percents_in_cluster) %>% 
  mutate(normed_percent_file_in_cluster = normed_fraction_file_in_cluster*100)

#pull out matrix of just required info for plotting
normed_fractions <- cluster_plot_values %>% 
  select(Cluster, File_Num, normed_fraction_file_in_cluster)

#add the day labels
md.df <- as.data.frame(md) #Convert md to dataframe (needed to run with metadata)
days_fractions <- normed_fractions %>% 
  left_join(md.df[,c("File_Num", "Age", "SampleType", METADATA.TYPE)], by="File_Num")

#sum each day together, this is specific for my plots since they have multiple replicates
summed_fractions <- ddply(days_fractions,.(Cluster,Age,SampleType,Tissue),summarize,sum=sum(normed_fraction_file_in_cluster),number=length(Age))

ordered_plotting <- summed_fractions %>% 
  left_join(ord_name_color, by="Cluster")

## Add marker expression to plotting data frame ===============================================
# Note: Ugly and inefficient (reloads stuff) because lazily copied from another script ;)
METADATA.FILENAME.2 <- paste0("/",METADATA.FILENAME)
metadata <- read.metadata(INPUT.FOLDER,METADATA.FILENAME.2)
concat.transformed <- fread(paste0(INPUT.FOLDER, CONCAT.TRANSFORMED.FILENAME), stringsAsFactors = F)
CLUSTERS.FILENAME.2 <- paste0("/",CLUSTERS.FILENAME)
clusters.in <- read.clusters(INPUT.FOLDER,CLUSTERS.FILENAME.2)
filenums <- read.table(FILES.FILENAME, header=TRUE, sep=",", check.names=FALSE)[,]

#Prepare sample metadata for each group (specified in SAMPLE.GROUP.ID) and clusters for addition to dataframe
metadata_key <- as.character(metadata[,SAMPLE.GROUP.ID])
names(metadata_key) <- 1:length(metadata[,SAMPLE.GROUP.ID])
metadata.in <- metadata_key[filenums]
metadata.in <- as.data.frame(metadata.in)
colnames(metadata.in) <- SAMPLE.GROUP.ID

#Prepare sample metadata (specified in METADATA.TYPE) and clusters for addition to dataframe
clusters.in.named <- as.data.frame(clusters.in)
colnames(clusters.in.named) <- "Cluster"

#Combine expression data, umap layout, metadata, and cluster ids into full dataframe
full_exprs_df <- do.call("cbind", list(concat.transformed, metadata.in, clusters.in.named))

#Isolate expression data for marker to plot 
cols.keep <- c(SAMPLE.GROUP.ID,"Cluster",MARKER.TO.PLOT)
marker_exprs_df <- subset(full_exprs_df, select = cols.keep)

# Get median expression values for MARKER.TO.PLOT for each SAMPLE.GROUP.ID
expression.out <- vector(mode = "list") #Initialize list for median expresion data

for (i in SAMPLE.GROUP.ORDER) {
  sample_exprs_df <- subset(marker_exprs_df, SampleType==i)
  cluster.ids <- unique(sample_exprs_df$Cluster) # Get list of unique clusters
  median_exprs_df <- vector(mode = "list", length = length(cluster.ids))
  for (k in cluster.ids) {
    k_exprs_df <- subset(sample_exprs_df, Cluster == k)
    if (is.na(median(k_exprs_df$Ki67_152Sm))) {
      median_exprs <- k_exprs_df %>%
        mutate(Median = 0)
    } else {
      median_exprs <- k_exprs_df %>%
        mutate(Median = median(k_exprs_df$Ki67_152Sm))
    }
    median_exprs <- median_exprs[1,] # Take only first row (summary)
    median_exprs_df[[k]] <- median_exprs[,-3] # Keep everything but raw data column
  }
  expression.out[[i]] <- do.call(rbind, median_exprs_df)
}
exprs_df <- do.call(rbind, expression.out)
colnames(exprs_df)[colnames(exprs_df) == "Median"] <- MARKER.TO.PLOT # Rname column to marker

ordered_plotting <- ordered_plotting %>% right_join(exprs_df, by=c("SampleType","Cluster"))
ordered_plotting <- subset(ordered_plotting, select = -c(number,Color)) # Remove unecessary columns

# Make the cowplot =============================================================
#plotting trying to order by manually put in numbers
order_for_plotting <- ordered_plotting %>% 
  arrange(Plotting_Order)
order_for_plotting$Tissue <- factor(order_for_plotting$Tissue, levels = METADATA.ORDER)
order_for_plotting$Dot_Plot_Name <- factor(order_for_plotting$Dot_Plot_Name, levels=unique(order_for_plotting$Dot_Plot_Name)[unique(order_for_plotting$Plotting_Order)], ordered=TRUE)

#Remove rows with NA (in case only plotting certain clusters, otherwise does nothing)
order_for_plotting <- order_for_plotting[!is.na(order_for_plotting$Plotting_Order), ]

#inputs
dot.scale = 8
scale.min = 0
scale.max = 1 # Can change to max(order_for_plotting$sum) to reduce scale
scale.by = 'size'

scale.func <- switch(
  EXPR = scale.by,
  'size' = scale_size,
  'radius' = scale_radius,
  stop("'scale.by' must be either 'size' or 'radius'")
)


#Generate plot of all markers colored by particular marker
plot <- ggplot(data = order_for_plotting, mapping = aes_string(x = 'Age', y = 'Tissue')) + #put in matrix containing x (day) and y (cluster)
  geom_point(mapping = aes_string(size = 'sum', color = MARKER.TO.PLOT)) +
  scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.33,0.8,1)) + #scale_color_viridis()
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + #optional scaling method
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +  #sets x and y axis title to be blank
  guides(size = guide_legend(title = 'Percent Expressed')) +  #sets the legend title
  labs(x = 'Timepoint', y = 'Cluster') +
  theme_cowplot() +
  theme(plot.title = element_text(size=3), axis.text.x = element_text(angle = 70, hjust=1, size = 8)) +
  facet_wrap(~ Dot_Plot_Name, ncol = 7)

ggsave(OUTPUT.FILENAME.Cowplot, plot = plot, device = OUTPUT.DEVICE, width = 16, height = 3.3) # Max width and height ~30

print("Finish CowPlot_ByCluster_ColoredBy_Marker.R")
