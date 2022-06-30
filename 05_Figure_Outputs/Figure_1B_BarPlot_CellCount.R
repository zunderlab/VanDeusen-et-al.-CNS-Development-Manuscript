#Zunder Lab, University of Virginia
#Amy Van Deusen
#25 September, 2021
#Generate plots from Cell Count Data

print("Start CellCountBarPlot.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(ggplot2)

## INPUT PARAMETERS ===============================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
COUNT.FILENAME <- "CellCounts.csv"
COLOR.FILENAME <- "Colors.csv"
OUTPUT.FILENAME <- "/CellCounts.pdf"
OUTPUT.DEVICE <- "pdf" # "png", "pdf", "jpeg", etc

## READ IN NECESSARY FILES ========================================================
counts <- read.csv(COUNT.FILENAME)
colors <- data.frame(read.csv(paste0(COLOR.FILENAME)))
colors.in <- as.vector(colors$color)
output_filename <- paste0(INPUT.FOLDER, OUTPUT.FILENAME)

## GENERATE AND OUTPUT STACKED BAR PLOT ================================================================
p <- ggplot(counts, aes(fill=Tissue, y=Count, x=Age)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = colors.in) + theme_classic() + theme(axis.text.x = element_text(face="bold", size=10, vjust=0.5, angle = 45))

ggsave(output_filename, p,  device = OUTPUT.DEVICE)

print("Finish CellCountBarPlot.R")
