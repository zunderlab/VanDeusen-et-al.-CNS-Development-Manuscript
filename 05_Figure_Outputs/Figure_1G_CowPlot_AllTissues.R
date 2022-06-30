#Kristen Fread, University of Virginia
#20 October 2021
#Generate cowplot with cluster abundances

print("Start CowPlot_AllTissues.R")

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
FILES.FILENAME <- "filenums.csv" 
METADATA.FILENAME <- "metadata.csv"
CLUSTERS.FILENAME <- "cluster_R2_assigns.csv"
OUTPUT.FILENAME.Totalcells <- "Total_cells_per_timepoint.pdf"
OUTPUT.FILENAME.Cowplot <- "CowPlot_All.pdf"
CLUSTER.NAMES.ORDER.COLOR <- "Dot_Plot_Parameters.csv" #has to be individually crafted with labels and desired colors
OUTPUT.DEVICE <- "pdf" # Choose "png" or "pdf"

## Read needed files
fil_nums <- fread(paste0(INPUT.FOLDER, "/", FILES.FILENAME), stringsAsFactors = F)
colnames(fil_nums) <- "File"
clusters.in <- read.csv(paste0(INPUT.FOLDER,"/",CLUSTERS.FILENAME))
colnames(clusters.in) <- "Cluster"
# read metadata
md <- fread(paste0(INPUT.FOLDER, "/", METADATA.FILENAME), stringsAsFactors = F)
# assign file numbering
md <- cbind(md, "File" = 1:nrow(md)) 
#read this in for dot plot ordering
ord_name_color <- fread(paste0(INPUT.FOLDER, "/", CLUSTER.NAMES.ORDER.COLOR), stringsAsFactors = F)


## Prep dataframe for plotting
plotting_df <- as.data.frame(cbind(fil_nums,clusters.in))
colnames(plotting_df) <- c("File_Num", "Cluster")

#join metadata to the dataframe with file and cluster number
cluster_file_metadata <- plotting_df %>% 
  left_join(md, by="File_Num")

#calculate number of cells per age
totalcellstimepoint <- cluster_file_metadata %>% 
  group_by(Age) %>% 
  dplyr::summarise(num_cells=n())

total_cells <- sum(totalcellstimepoint$num_cells)

####### Plot the total cells per timepoint ####### =========================================
p <- ggplot(totalcellstimepoint, aes(Age, num_cells, fill = Age)) + geom_col() + 
  scale_fill_viridis(discrete=TRUE, name = "Timepoint", guide="none") +
  theme_bw() +
  ggtitle("Total cells per timepoint") +
  labs(x = "Timepoint", y = "Number of cells") +
  theme(plot.title = element_text(size=30), axis.title = element_text(size = 30), 
        axis.text = element_text(size = 22), axis.text.x = element_text(angle = 70, hjust=1)) 

# save plot
ggsave(OUTPUT.FILENAME.Totalcells, plot = p, device = OUTPUT.DEVICE)

##Calculate values for Dot Plot ============================================================
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
days_fractions <- normed_fractions %>% 
  left_join(md[,c("File_Num", "Age")], by="File_Num")

#sum each day together, this is specific for my plots since they have multiple replicates
summed_fractions <- ddply(days_fractions,.(Cluster,Age),summarize,sum=sum(normed_fraction_file_in_cluster),number=length(Age))

ordered_plotting <- summed_fractions %>% 
  left_join(ord_name_color, by="Cluster")

#make the dotplot
#inputs
dot.scale = 8
scale.min = 0
scale.max = 1
scale.by = 'size'

scale.func <- switch(
  EXPR = scale.by,
  'size' = scale_size,
  'radius' = scale_radius,
  stop("'scale.by' must be either 'size' or 'radius'")
)

#plotting trying to order by manually put in numbers
order_for_plotting <- ordered_plotting %>% 
  arrange(Plotting_Order)
order_for_plotting$Dot_Plot_Name <- factor(order_for_plotting$Dot_Plot_Name, levels=rev
                                           (unique(order_for_plotting$Dot_Plot_Name)[unique(order_for_plotting$Plotting_Order)]), ordered=TRUE)
coloring <- rev(unique(order_for_plotting$Color))

order_for_plotting <- order_for_plotting[!is.na(order_for_plotting$Age), ]

plot <- ggplot(data = order_for_plotting, mapping = aes_string(x = 'Age', y = 'Dot_Plot_Name')) + #put in matrix containing x (day) and y (cluster)
  geom_point(mapping = aes_string(size = 'sum', color = 'Dot_Plot_Name')) +
  scale_color_manual(values = coloring) +
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + #optional scaling method
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +  #sets x and y axis title to be blank
  guides(size = guide_legend(title = 'Percent Expressed')) +  #sets the legend title
  labs(x = 'Timepoint', y = 'Cluster') +
  theme_cowplot() +
  theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 70, hjust=1))

ggsave(OUTPUT.FILENAME.Cowplot, plot = plot, device = OUTPUT.DEVICE, width = 28.5, height = 14) # Max width and height ~30

print("Finish CowPlot_AllTissues.R")
