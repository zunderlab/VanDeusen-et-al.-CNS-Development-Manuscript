## Get cell counts for each sample and perform statistical analyses

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggthemes)
library(reshape2)

# Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
SAMPLE.CELLS.FILENAME <- "Total_Cells_ByClass.csv" # Can use this instead of FILENUMS.FILENAME
METADATA.FILENAME <- "metadata_master.csv" 
OUTPUT.BASENAME <- "CellClasses_Plot_Cluster.pdf"
OUTPUT.DEVICE <- "pdf"
COLOR.PALETTE <- c("Brain"="#F242F5", "Cortex"="#4E42F5","Diencephalon"="#18C92D",
                   "Midbrain"="#FFC316","Hindbrain"="#FF4E41")
INDIVIDUAL.SCALES <- TRUE # TRUE = each plot individually scaled for Y axis; FALSE = all plot have same Y scale

# Read in necessary files
cell_counts <- read.csv(paste0(INPUT.FOLDER,"/",SAMPLE.CELLS.FILENAME))
metadata <- read.csv(METADATA.FILENAME)
metadata.in <- metadata[,5:7]

# Generate initial plot data frame
cell_counts_init <- cbind(cell_counts,metadata.in)

# Calculate percentages for each class
cell_counts_pct <- cell_counts_init %>%
  mutate(NSC = (NSC/Total_Cells)*100) %>% 
  mutate(NPCs = (NPCs/Total_Cells)*100) %>%
  mutate(Neurons = (Neurons/Total_Cells)*100) %>%
  mutate(Interneurons = (Interneurons/Total_Cells)*100) %>%
  mutate(RGC_GP = (RGC_GP/Total_Cells)*100) %>%
  mutate(Olig2_NPCs = (Olig2_NPCs/Total_Cells)*100) %>%
  mutate(Olig2_GPs = (Olig2_GPs/Total_Cells)*100) %>%
  mutate(OPCs = (OPCs/Total_Cells)*100) %>%
  mutate(Nonneural = (Nonneural/Total_Cells)*100) %>%
  mutate(Endo = (Endo/Total_Cells)*100) %>%
  mutate(Mural = (Mural/Total_Cells)*100) %>%
  mutate(Microglia = (Microglia/Total_Cells)*100)
  
# Generate dataframes for plots comparing relative abundances of each class for each tissue
tissues <- unique(cell_counts_pct$Tissue)
cluster_classes <- colnames(cell_counts_pct[,2:13]) # This line will need to be edited based on number of classes!
sample_ages <- unique(cell_counts_pct$Age)

# Fill plotting dataframe with summary data and ranges for each metadata type
for (i in cluster_classes) {
  cols_keep <- c("Age","Tissue",i)
  cluster_pct <- subset(cell_counts_pct, select=cols_keep) # Isolate each cluster
  cluster_replicates <- melt(cluster_pct)
  colnames(cluster_replicates) <- c("Age","Tissue","Sample_Class","Abundance_Pct")
  
  # Set y axis scaling
  if (isTRUE(INDIVIDUAL.SCALES)) {
    max_value <- max(cluster_replicates$Abundance_Pct) + 0.02
  } else {
    max_value <- max(cell_counts_pct[,2:13]) # This line will need to be edited based on number of classes!
  }
  
  # Set up basic plots
  p <- ggplot(cluster_replicates, aes(x=Age, y=Abundance_Pct, colour = Tissue, group = Tissue))
  age_labels <- unique(cluster_replicates$Age)
  plot <- p + geom_point(aes(colour=Tissue,y=Abundance_Pct), shape=1, size=5, stroke = 3) +
    scale_x_discrete(labels = age_labels) +
    xlab("Age") + ylab("% Abundance") +
    ylim(0,max_value) +
    stat_smooth(aes(fill = Tissue), size = 3) +
    scale_colour_manual(name = "Tissue", values = COLOR.PALETTE) +
    scale_fill_manual(values = COLOR.PALETTE) +
    theme_hc()
  
  # Export plot
  OUTPUT.FILENAME <- gsub("Plot", "ScatterPlot", OUTPUT.BASENAME)
  OUTPUT.FILENAME <- gsub("Cluster", i, OUTPUT.FILENAME)
  ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)
}



