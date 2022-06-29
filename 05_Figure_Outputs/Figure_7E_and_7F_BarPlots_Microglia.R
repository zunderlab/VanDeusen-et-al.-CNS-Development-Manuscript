#Amy Van Deusen, University of Virginia
#2022
#Get cell counts for each sample and perform statistical analyses
#Output bar plots

print("Start BarPlot_Microglia.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata_cluster.csv"
TOTAL.CELLS.FILENAME <- "Total_Cells.csv" # Can use this instead of FILENUMS.FILENAME
SAMPLE.CELLS.FILENAME <- "Total_PhagocyticCells.csv"
FILENUMS.FILENAME <- NULL # Can use this instead of SAMPLE.CELLS.FILENAME
OUTPUT.FILENAME <- "PhagocyticCells_BarPlot.pdf"
OUTPUT.DEVICE <- "pdf"
COLOR.PALETTE <- c("#F242F5", "#4E42F5", "#18C92D","#FFC316","#FF4E41")

# Read in metadata
metadata <- read.csv(paste0(INPUT.FOLDER,"/",METADATA.FILENAME))

if (is.null(FILENUMS.FILENAME)) {
  sample_counts_filename <- paste0(INPUT.FOLDER,"/",SAMPLE.CELLS.FILENAME)
  file_counts_df <- read.csv(sample_counts_filename)
  file_counts_df <- file_counts_df[,-1] 
} else {
  # Read in file numbers for sampled cells
  filenums <- read.csv(FILENUMS.FILENAME)
  colnames(filenums) <- "File"
  file_counts <- table(filenums$File) # Extract counts for each file
  file_counts_df <- as.data.frame(file_counts) # Convert to dataframe
  file_counts_df <- file_counts_df[,-1]
}
#Generate dataframe containing cell counts and metadata
cell_counts <- cbind(metadata, file_counts_df) # Combine counts and metadata
colnames(cell_counts)[7] <- "Cell_Count"
sample_types <- unique(cell_counts$SampleType)


#Calculate relative cell abundances (instead of comparing raw counts)
cell_totals <- as.data.frame(read.csv(paste0(INPUT.FOLDER,"/",TOTAL.CELLS.FILENAME))) # Input total cell numbers from each sample type 
cell_totals <- cell_totals[,-1]
cell_counts <- cbind(cell_counts, cell_totals) # Add total cells coulmn to df
colnames(cell_counts)[8] <- "Total_Cells"

cell_counts <- cell_counts %>%
  mutate(Abundance = Cell_Count/Total_Cells) %>% 
  mutate(Abundance_Pct = Abundance*100)


# Initialize plotting dataframe containing columns for sample range values
plot_df <- data.frame(SampleType = character(), 
                      Age = character(),
                      Tissue = character(),
                      Mean = double(),
                      Mean_Pct = double(),
                      Min = double(),
                      Min_Pct = double(),
                      Max = double(),
                      Max_Pct = double(),
                      SD = double(),
                      SD_Pct = double()) # Initialze dataframe 

# Fill plotting dataframe with summary data and ranges for each metadata type
for (i in sample_types) {
  iteration <- match(i, sample_types) # Iteration = new row number for plot_df
  filtered_cc <- filter(cell_counts, SampleType == i) # Isolate each sample type
  sample_age <- unique(filtered_cc$Age) # Obtain value for age
  sample_tissue <- unique(filtered_cc$Tissue) # Obtain value for tissue
  sample_mean <- mean(filtered_cc$Abundance)
  sample_mean_pct <- mean(filtered_cc$Abundance_Pct)
  sample_min <- min(filtered_cc$Abundance) # Obtain minimum abundance
  sample_min_pct <- min(filtered_cc$Abundance_Pct) # Obtain minimum abundance percentage
  sample_max <- max(filtered_cc$Abundance) # Obtain maximum abundance
  sample_max_pct <- max(filtered_cc$Abundance_Pct) # Obtain maximum abundance percentage
  sample_sd <- sd(filtered_cc$Abundance) # Calculate standard deviation
  sample_sd_pct <- sd(filtered_cc$Abundance_Pct) # Calculate standard deviation as percentage
  sample_summary <- c(i, sample_age, sample_tissue, sample_mean, sample_mean_pct, 
                      sample_min, sample_min_pct, sample_max, sample_max_pct, sample_sd, sample_sd_pct)
  plot_df[nrow(plot_df)+1,] <- sample_summary # Add summary data to ploting dataframe
}

# Convert values in plotting dataframe to numeric for ggplot
plot_df$Mean <- as.numeric(plot_df$Mean) # Converts character to numeric
plot_df$Mean_Pct <- as.numeric(plot_df$Mean_Pct) # Converts character to numeric
plot_df$Min <- as.numeric(plot_df$Min) # Converts character to numeric
plot_df$Min_Pct <- as.numeric(plot_df$Min_Pct) # Converts character to numeric
plot_df$Max <- as.numeric(plot_df$Max) # Converts character to numeric
plot_df$Max_Pct <- as.numeric(plot_df$Max_Pct) # Converts character to numeric
plot_df$SD <- as.numeric(plot_df$SD) # Converts character to numeric
plot_df$SD_Pct <- as.numeric(plot_df$SD_Pct) # Converts character to numeric
plot_df$Tissue <- factor(plot_df$Tissue, levels = c('Brain', 'Cortex', 'Diencephalon', 'Midbrain', 'Hindbrain')) # Converts character to a factor

# Generate range plots with various aesthetics
p <- ggplot(plot_df, aes(x=Age, y=Mean_Pct, colour = Tissue, group = Tissue))
plot <- p + geom_bar(aes(fill = Tissue), colour="black", stat = "identity",size = 0, width = 1, position = position_dodge2(width=0.8, preserve ='single')) +
  geom_errorbar(aes(ymin = Min_Pct,ymax = Max_Pct),colour = "black", position = position_dodge2(width=0.8, preserve='single'), width=1, size=0.2) +
  xlab("Age") + ylab("% Abundance") +
  scale_fill_manual(values = COLOR.PALETTE) + 
  scale_y_continuous(breaks=c(0,0.05,0.10,0.15)) +
  theme_minimal() + facet_grid(cols = vars(Tissue),  scales = "free_x", space = "free_x")
plot

ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)

#Copy plot_df to new dataframe
plot_df_interacting <- plot_df


## Process data for Noninteracting Cells =================================================

SAMPLE.CELLS.FILENAME <- "Total_NoninteractingMicroglia.csv"
OUTPUT.FILENAME <- "NoninteractingMicroglia_BarPlot.pdf"

if (is.null(FILENUMS.FILENAME)) {
  sample_counts_filename <- paste0(INPUT.FOLDER,"/",SAMPLE.CELLS.FILENAME)
  file_counts_df <- read.csv(sample_counts_filename)
  file_counts_df <- file_counts_df[,-1] 
} else {
  # Read in file numbers for sampled cells
  filenums <- read.csv(FILENUMS.FILENAME)
  colnames(filenums) <- "File"
  file_counts <- table(filenums$File) # Extract counts for each file
  file_counts_df <- as.data.frame(file_counts) # Convert to dataframe
  file_counts_df <- file_counts_df[,-1]
}
#Generate dataframe containing cell counts and metadata
cell_counts <- cbind(metadata, file_counts_df) # Combine counts and metadata
colnames(cell_counts)[7] <- "Cell_Count"
sample_types <- unique(cell_counts$SampleType)

#Calculate relative cell abundances (instead of comparing raw counts)
cell_totals <- as.data.frame(read.csv(paste0(INPUT.FOLDER,"/",TOTAL.CELLS.FILENAME))) # Input total cell numbers from each sample type 
cell_totals <- cell_totals[,-1]
cell_counts <- cbind(cell_counts, cell_totals) # Add total cells coulmn to df
colnames(cell_counts)[8] <- "Total_Cells"

cell_counts <- cell_counts %>%
  mutate(Abundance = Cell_Count/Total_Cells) %>% 
  mutate(Abundance_Pct = Abundance*100)

# Initialize plotting dataframe containing columns for sample range values
plot_df <- data.frame(SampleType = character(), 
                      Age = character(),
                      Tissue = character(),
                      Mean = double(),
                      Mean_Pct = double(),
                      Min = double(),
                      Min_Pct = double(),
                      Max = double(),
                      Max_Pct = double(),
                      SD = double(),
                      SD_Pct = double()) # Initialze dataframe 

# Fill plotting dataframe with summary data and ranges for each metadata type
for (i in sample_types) {
  iteration <- match(i, sample_types) # Iteration = new row number for plot_df
  filtered_cc <- filter(cell_counts, SampleType == i) # Isolate each sample type
  sample_age <- unique(filtered_cc$Age) # Obtain value for age
  sample_tissue <- unique(filtered_cc$Tissue) # Obtain value for tissue
  sample_mean <- mean(filtered_cc$Abundance)
  sample_mean_pct <- mean(filtered_cc$Abundance_Pct)
  sample_min <- min(filtered_cc$Abundance) # Obtain minimum abundance
  sample_min_pct <- min(filtered_cc$Abundance_Pct) # Obtain minimum abundance percentage
  sample_max <- max(filtered_cc$Abundance) # Obtain maximum abundance
  sample_max_pct <- max(filtered_cc$Abundance_Pct) # Obtain maximum abundance percentage
  sample_sd <- sd(filtered_cc$Abundance) # Calculate standard deviation
  sample_sd_pct <- sd(filtered_cc$Abundance_Pct) # Calculate standard deviation as percentage
  sample_summary <- c(i, sample_age, sample_tissue, sample_mean, sample_mean_pct, 
                      sample_min, sample_min_pct, sample_max, sample_max_pct, sample_sd, sample_sd_pct)
  plot_df[nrow(plot_df)+1,] <- sample_summary # Add summary data to ploting dataframe
}

# Convert values in plotting dataframe to numeric for ggplot
plot_df$Mean <- as.numeric(plot_df$Mean) # Converts character to numeric
plot_df$Mean_Pct <- as.numeric(plot_df$Mean_Pct) # Converts character to numeric
plot_df$Min <- as.numeric(plot_df$Min) # Converts character to numeric
plot_df$Min_Pct <- as.numeric(plot_df$Min_Pct) # Converts character to numeric
plot_df$Max <- as.numeric(plot_df$Max) # Converts character to numeric
plot_df$Max_Pct <- as.numeric(plot_df$Max_Pct) # Converts character to numeric
plot_df$SD <- as.numeric(plot_df$SD) # Converts character to numeric
plot_df$SD_Pct <- as.numeric(plot_df$SD_Pct) # Converts character to numeric
plot_df$Tissue <- factor(plot_df$Tissue, levels = c('Brain', 'Cortex', 'Diencephalon', 'Midbrain', 'Hindbrain')) # Converts character to a factor

# Generate range plots with various aesthetics
p <- ggplot(plot_df, aes(x=Age, y=Mean_Pct, colour = Tissue, group = Tissue))
plot <- p + geom_bar(aes(fill = Tissue), colour="black", stat = "identity",size = 0, width = 1, position = position_dodge2(width=0.8, preserve ='single')) +
  geom_errorbar(aes(ymin = Min_Pct,ymax = Max_Pct),colour = "black", position = position_dodge2(width=0.8, preserve='single'), width=1, size=0.2) +
  xlab("Age") + ylab("% Abundance") +
  scale_fill_manual(values = COLOR.PALETTE) +
  theme_minimal() + facet_grid(cols = vars(Tissue),  scales = "free_x", space = "free_x")
plot

ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)

#Copy plot_df to new dataframe
plot_df_noninteracting <- plot_df

## Combine plot dataframes to compare relative abundances of interacting and noninteracting microglia
# Add columns indicating cell type
plot_df_interacting$CellType <- "Interacting"
plot_df_noninteracting$CellType <- "Noninteracting"

# Combine plot data frames
plot_df_all <- rbind(plot_df_interacting, plot_df_noninteracting)
plot_df_brain <- filter(plot_df_all, Tissue == "Brain")
plot_df_cortex <- filter(plot_df_all, Tissue == "Cortex")
plot_df_dien <- filter(plot_df_all, Tissue == "Diencephalon")
plot_df_mid <- filter(plot_df_all, Tissue == "Midbrain")
plot_df_hind <- filter(plot_df_all, Tissue == "Hindbrain")

# Export paired line graph comparing flow and CyTOF doublets for each sample type
plot <- ggpaired(plot_df_all, x = 'CellType', y = 'Mean_Pct',
                 color = 'CellType', line.color = "gray", line.size = 0.4,
                 palette = "npg") + 
  stat_compare_means(paired = TRUE) +
  xlab("Method") + ylab("Doublets (% Total Cells)")
plot
OUTPUT.FILENAME <- "MicrogliaComparison_PairedBoxPlot.pdf"
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)

print("Finish BarPlot_Microglia.R")
