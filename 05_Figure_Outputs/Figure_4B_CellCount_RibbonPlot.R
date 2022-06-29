#Amy Van Deusen, University of Virginia
#2021
#Get cell counts for each sample and perform statistical analyses
#Output line graph with ribbons indicating range

print("Start CellCount_RibbonPlot.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(dplyr)
library(ggplot2)
library(tidyverse)

# Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
FILENUMS.FILENAME <- "filenums.csv"
METADATA.FILENAME <- "metadata.csv"
TOTAL.CELLS.FILENAME <- "Total_Cells.csv"
OUTPUT.FILENAME <- "CellCounts_RibbonPlot.pdf"
OUTPUT.DEVICE <- "pdf"
COLOR.PALETTE <- c("#F242F5", "#4E42F5", "#18C92D","#FFC316", "#FF4E41")

# Read in metadata
metadata <- read.csv(paste0(INPUT.FOLDER,"/",METADATA.FILENAME))

# Read in file numbers
filenums <- read.csv(FILENUMS.FILENAME)
colnames(filenums) <- "File"
file_counts <- table(filenums$File) # Extract counts for each file
file_counts_df <- as.data.frame(file_counts) # Convert to dataframe
file_counts_df <- file_counts_df[,-1]

#Generate dataframe containing cell counts and metadata
cell_counts <- cbind(metadata, file_counts_df) # Combine counts and metadata
colnames(cell_counts)[7] <- "Cell_Count" # Rename column containing cell counts
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
                      Max_Pct = double()) # Initialze dataframe 

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
  sample_summary <- c(i, sample_age, sample_tissue, sample_mean, sample_mean_pct, 
                      sample_min, sample_min_pct, sample_max, sample_max_pct)
  plot_df[nrow(plot_df)+1,] <- sample_summary # Add summary data to ploting dataframe
}

# Convert values in plotting dataframe to numeric for ggplot
plot_df$Mean <- as.numeric(plot_df$Mean) # Converts character to numeric
plot_df$Mean_Pct <- as.numeric(plot_df$Mean_Pct) # Converts character to numeric
plot_df$Min <- as.numeric(plot_df$Min) # Converts character to numeric
plot_df$Min_Pct <- as.numeric(plot_df$Min_Pct) # Converts character to numeric
plot_df$Max <- as.numeric(plot_df$Max) # Converts character to numeric
plot_df$Max_Pct <- as.numeric(plot_df$Max_Pct) # Converts character to numeric
plot_df$Tissue <- as.factor(plot_df$Tissue) # Converts character to a factor

# Generate range plots with various aesthetics
p <- ggplot(plot_df, aes(x=Age, y=Mean_Pct, colour = Tissue, group = Tissue))

# As a ribbon plot
max_value <- max(plot_df$Max_Pct) + 0.02
age_labels <- unique(plot_df$Age)
plot <- p + geom_line(aes(colour=factor(Tissue)), size = 1) + 
  scale_x_discrete(labels = age_labels) +
  xlab("Age") + ylab("% Abundance") +
  ylim(0,max_value) +
  geom_ribbon(aes(ymin=Min_Pct, ymax=Max_Pct, fill = factor(Tissue)), 
              alpha = 0.25, colour = NA) +
  scale_colour_manual(name = "Tissue", values = COLOR.PALETTE) +
  scale_fill_manual(values = COLOR.PALETTE) +
  theme_minimal() 

# Export plot
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)

print("Finish CellCount_RibbonPlot.R")
