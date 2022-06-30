#Amy Van Deusen, University of Virginia
#2022
#Get cell counts for each sample and perform statistical analyses
#Output bar plots

print("Start BarPlot_CellDoublets.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Input parameters
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata.csv"
TOTAL.CELLS.FILENAME <- "Ungated_Cells.csv"
SAMPLE.CELLS.FILENAME <- "Total_BCDistxWidth.csv" # Can use this instead of FILENUMS.FILENAME
FILENUMS.FILENAME <- NULL # Can use this instead of SAMPLE.CELLS.FILENAME
OUTPUT.FILENAME <- "BCDistxWidth_BarPlot_ByTissue.pdf"
OUTPUT.DEVICE <- "pdf"
COLOR.PALETTE <- c("#F242F5", "#4E42F5", "#18C92D", "#FFC316", "#FF4E41")

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
colnames(cell_counts)[8] <- "Cell_Count" # Rename column containing cell counts
sample_types <- unique(cell_counts$SampleType)

#Calculate relative cell abundances (instead of comparing raw counts)
cell_totals <- as.data.frame(read.csv(paste0(INPUT.FOLDER,"/",TOTAL.CELLS.FILENAME))) # Input total cell numbers from each sample type 
cell_totals <- cell_totals[,-1]
cell_counts <- cbind(cell_counts, cell_totals) # Add total cells coulmn to df
colnames(cell_counts)[9] <- "Total_Cells"
cell_counts <- cell_counts %>%
  mutate(Abundance = Cell_Count/Total_Cells) %>% 
  mutate(Abundance_Pct = Abundance*100)

# Initialize plotting dataframe containing columns for sample range values
plot_df <- data.frame(SampleType = character(),
                      Tissue_Method = character(),
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
  tissue_method <- unique(filtered_cc$Tissue_Method) # Obtain value for tissue and method type (for plotting)
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
  sample_summary <- c(i, tissue_method, sample_age, sample_tissue, sample_mean, sample_mean_pct, 
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
plot_df$Tissue_Method <- as.factor(plot_df$Tissue_Method) # Converts character to a factor


#Copy to CyTOF data frame
cell_counts_CyTOF <- cell_counts # Make a copy of CyTOF cell counts
plot_df_CyTOF <- plot_df
plot_df_CyTOF <- plot_df_CyTOF %>%
  mutate(method.type = "CyTOF")

# Output plot
p <- ggplot(plot_df, aes(x=Age, y=Mean_Pct, colour = Tissue, group = Tissue))
plot <- p + geom_bar(aes(fill = Tissue), colour="black", stat = "identity",size = 0, width = 1, position = position_dodge2(width=0.8, preserve ='single')) +
  geom_errorbar(aes(ymin = Min_Pct,ymax = Max_Pct),colour = "black", position = position_dodge2(width=0.8, preserve='single'), width=1, size=0.2) +
  xlab("Age") + ylab("% Abundance") +
  scale_fill_manual(values = COLOR.PALETTE) +
  theme_minimal() + facet_grid(cols = vars(Tissue), scales = "free_x", space = "free_x")
plot
# Export plot
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)


## Generate Flow doublte dataframe ==================================================
METADATA.FILENAME <- "metadata_Flow.csv"
TOTAL.CELLS.FILENAME <- "Total_Cells_Flow.csv"
SAMPLE.CELLS.FILENAME <- "Doublets_Flow.csv" # Can use this instead of FILENUMS.FILENAME
FILENUMS.FILENAME <- NULL # Can use this instead of SAMPLE.CELLS.FILENAME
OUTPUT.FILENAME <- "FSCxDRAQ7_BarPlot_ByTissue.pdf"

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
colnames(cell_counts)[8] <- "Cell_Count" # Rename column containing cell counts
sample_types <- unique(cell_counts$SampleType)

#Calculate relative cell abundances (instead of comparing raw counts)
cell_totals <- as.data.frame(read.csv(paste0(INPUT.FOLDER,"/",TOTAL.CELLS.FILENAME))) # Input total cell numbers from each sample type 
cell_totals <- cell_totals[,-1]
cell_counts <- cbind(cell_counts, cell_totals) # Add total cells coulmn to df
colnames(cell_counts)[9] <- "Total_Cells"
cell_counts <- cell_counts %>%
  mutate(Abundance = Cell_Count/Total_Cells) %>% 
  mutate(Abundance_Pct = Abundance*100)

# Initialize plotting dataframe containing columns for sample range values
plot_df <- data.frame(SampleType = character(), 
                      Tissue_Method = character(),
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
  tissue_method <- unique(filtered_cc$Tissue_Method) # Obtain value for tissue and method type (for plotting)
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
  sample_summary <- c(i, tissue_method, sample_age, sample_tissue, sample_mean, sample_mean_pct, 
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
plot_df$Tissue <-  factor(plot_df$Tissue, levels = c('Brain', 'Cortex', 'Diencephalon', 'Midbrain', 'Hindbrain')) # Converts character to a factor
plot_df$Tissue_Method <- as.factor(plot_df$Tissue_Method) # Converts character to a factor

# Copy to CyTOF data frame
cell_counts_Flow <- cell_counts # Make a copy of flow cell counts
plot_df_Flow <- plot_df
plot_df_Flow <- plot_df_Flow %>%
  mutate(method.type = "Flow")

# Output plot
p <- ggplot(plot_df, aes(x=Age, y=Mean_Pct, colour = Tissue, group = Tissue))
plot <- p + geom_bar(aes(fill = Tissue), colour="black", stat = "identity",size = 0, width = 1, position = position_dodge2(width=0.8, preserve ='single')) +
  geom_errorbar(aes(ymin = Min_Pct,ymax = Max_Pct),colour = "black", position = position_dodge2(width=0.8, preserve='single'), width=1, size=0.2) +
  xlab("Age") + ylab("% Abundance") +
  scale_fill_manual(values = COLOR.PALETTE) +
  theme_minimal() + facet_grid(cols = vars(Tissue), scales = "free_x", space = "free_x")
plot
# Export plot
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)


## Combine flow and CyTOF data frames ============================================
plot_df_all <- rbind(plot_df_Flow, plot_df_CyTOF)
COLOR.PALETTE <- c("#F242F5", "#FFCCF2","#4E42F5", "#99C2FF","#18C92D","#ADEBAD","#FFC316","#FFEB99","#FF4E41", "#FFAD99")
names(COLOR.PALETTE) <- c('Brain_CyTOF', 'Brain_Flow', 'Cortex_CyTOF', 'Cortex_Flow', 'Dien_CyTOF', 'Dien_Flow',
                          'Mid_CyTOF', 'Mid_Flow', 'Hind_CyTOF', 'Hind_Flow')

# Export bar graph comparing flow and CyTOF doublets, grouped by tissue
p <- ggplot(plot_df_all, aes(x=Age, y=Mean_Pct, colour = Tissue_Method, group = Age, dodge=factor(Tissue)))
plot <- p + geom_bar(aes(fill = Tissue_Method), colour="black", stat = "identity",
                     size = 0, width = 1, position = position_dodge2(width=2, preserve = c('single'))) +
  geom_errorbar(aes(ymin = Mean_Pct-SD_Pct,ymax = Mean_Pct+SD_Pct),colour = "black", 
                position = position_dodge2(width=2, preserve=c('single')), width=1, size=0.2) +
  xlab("Age") + ylab("% Abundance") + 
  scale_fill_manual(values = COLOR.PALETTE) +
  scale_y_continuous(limits = c(0,100), breaks=c(0,20,40,60,80,100)) +
  theme_minimal() + facet_grid(cols = vars(Tissue), scales = "free_x", space = "free_x")
OUTPUT.FILENAME <- "DoubletsComparison_BarPlot.pdf"
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)

# Export paired box plotcomparing flow and CyTOF doublets for each sample type
plot <- ggpaired(plot_df_all, x = 'method.type', y = 'Mean_Pct',
              color = 'method.type', line.color = "gray", line.size = 0.4,
              palette = "npg", ylim = c(0,100)) + 
  stat_compare_means(paired = TRUE) +
  xlab("Method") + ylab("Doublets (% Total Cells)")
OUTPUT.FILENAME <- "DoubletsComparison_PairedBoxPlot.pdf"
ggsave(OUTPUT.FILENAME, plot = plot, device = OUTPUT.DEVICE)


## Perform paired t-test for each sample type ============================================================
# Generate dataframes with mean CyTOF and flow values
ttest_init <- cell_counts_CyTOF[3:10] # Copy CyTOF dataframe as initial dataframe
colnames(ttest_init)[5:8] <- c("CyTOF_Doublets", "CyTOF_Total_Cells", 
                              "CyTOF_Abundance", "CyTOF_Abundance_Pct") # Rename columns with  CyTOF values
ttest_df <- merge(ttest_init,cell_counts_Flow[3:10]) # Add flow counts
colnames(ttest_df)[9:12] <- c("Flow_Doublets", "Flow_Total_Cells", 
                             "Flow_Abundance", "Flow_Abundance_Pct") # Rename columns with flow values

# Used paired t test to calculate p value of two methods to analyze doublets in each sample type
n.sample.types <- length(unique(ttest_df$SampleType))
sample.pvalues <- vector(mode = "list", length = n.sample.types) # Initiate list to collect p values
for (i in unique(ttest_df$SampleType)) {
  subset_df <- filter(ttest_df, SampleType == i)
  if (length(unique(subset_df$File_Num)) > 1) {
    doublets.cytof <- subset_df$CyTOF_Abundance
    doublets.flow <- subset_df$Flow_Abundance
    sample.pvalues[[i]]<- t.test(doublets.cytof, doublets.flow, paired = FALSE, alternative = 'two.sided')$p.value
    print(paste0("P value for ", i, ":", t.test(doublets.cytof, doublets.flow, paired = TRUE, alternative = 'two.sided')$p.value))
  } else {
    sample.pvalues[[i]]<- NA
    print(paste0("Not enough values to perform paired t test for ", i))
  }
}

sample.pvalues.unlist <- unlist(sample.pvalues)

print(sample.pvalues.unlist)

print("Finish BarPlot_CellDoublets.R")
