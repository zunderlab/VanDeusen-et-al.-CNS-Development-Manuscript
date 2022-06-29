#Amy Van Deusen, University of Virginia
#2021
#Generate bar plots and area plots (outputs both to see which more aesthetically pleasing)

print("Start BarPlot_Clusters_ByRegion.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)

print("Libraries loaded")

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata.csv" # Metadata for combined dataset
ABUNDANCES.FILENAME <- "ClusterAbundance_SampleType.csv" # Should be organized accodring to METADATA.TYPE
PLOT.PARAMETER.FILENAME <- "Dot_Plot_Parameters.csv"
METADATA.ORDER <- c("Hindbrain", "Midbrain", "Diencephalon", "Cortex", "Brain")
OUTPUT.BASENAME <- "Plot_Cluster_Abundances_XX.pdf" #XX Must be used in name for outputs!
INDIVIDUAL.OUTPUTS <- TRUE #TRUE = output each plot individually, FALSE = single file
OUTPUT.TYPE <- "pdf" # choose png or pdf


## Read in necessary files =====================================================
metadata <- read.csv(paste0(INPUT.FOLDER,"/",METADATA.FILENAME))
abundances <- read.table(ABUNDANCES.FILENAME, header=TRUE, sep=",", check.names=FALSE)
plot.parameters <- fread(paste0(INPUT.FOLDER, "/", PLOT.PARAMETER.FILENAME), stringsAsFactors = F)

# Generate plotting dataframe
metadata.in <- metadata[c("Age","Tissue")]
abundances.in <- t(abundances)
abundances.in <- abundances.in[-1,]
plot.df.init <- cbind(metadata.in, abundances.in)
plot.df.all <- melt(plot.df.init, id = c("Age", "Tissue"), value.name = "Abundance")
names(plot.df.all)[3] <- "Cluster"

#Generate a line graph for each cluster =======================================================
# Initialize plots and set custom color palette
barplotlist = list() # Generate empty list in which to collect plots
areaplotlist = list() # Generate empty list in which to collect plots

#Cycle through each cluster and build a line graph plot for each metadata type
for (i in unique(plot.df.all$Tissue)) {
  # Prepare plot dataframe
  brain.ids <- which(plot.df.all$Tissue == "Brain")
  i.ids <- which(plot.df.all$Tissue == i)
  sample.ids <- c(brain.ids,i.ids)
  plot.title <- i
  plot.filename <- sub("XX", i, OUTPUT.BASENAME)
  plot.df <- plot.df.all[c(sample.ids),]
  plot.df <- plot.df[ , c("Age","Cluster","Abundance")]
  plot.df <- transform(plot.df, Cluster = as.numeric(Cluster))
  plot.df.to.order <- plot.df %>% 
    left_join(plot.parameters, by="Cluster")
  plot.df.ordered <- plot.df.to.order %>% 
    arrange(Plotting_Order)
  plot.df.ordered$Dot_Plot_Name <- factor(plot.df.ordered$Dot_Plot_Name, levels=unique(plot.df.ordered$Dot_Plot_Name)[unique(plot.df.ordered$Plotting_Order)], ordered=TRUE)
  custom.colors <- unique(plot.df.ordered$Color)
  
  # Generate bar plots
  barplot <- ggplot(data=plot.df.ordered, aes(x=Age, y=Abundance, fill=Dot_Plot_Name)) + 
    geom_bar(aes(fill=Dot_Plot_Name),position = "fill", stat='identity') +
    scale_fill_manual(values = custom.colors) + 
    theme_grey() + labs(title = plot.title) +
    guides(shape = guide_legend(override.aes = list(size = 0.5)),
           color = guide_legend(override.aes = list(size = 0.1))) +
    guides(fill = guide_legend(ncol=1)) +
    theme_classic() +
    theme(axis.line=element_blank()) +
    theme(legend.key.size = unit(0.1, 'cm'), 
          legend.title = element_text(size=5),
          legend.text = element_text(size=5))
  barplotlist[[i]] = barplot
  
  # Generate area plots
  areaplot <- ggplot(data=plot.df.ordered, aes(x=Age, y=Abundance, group = Dot_Plot_Name, fill=Dot_Plot_Name)) + 
    geom_area(aes(fill=Dot_Plot_Name), position = "fill", stat='identity') +
    scale_fill_manual(values = custom.colors) + 
    theme_grey() + labs(title = plot.title) +
    guides(fill = guide_legend(ncol=1)) + 
    theme_classic() + theme(axis.line=element_blank()) +
    theme(legend.key.size = unit(0.1, 'cm'), 
          legend.title = element_text(size=5),
          legend.text = element_text(size=5))
  areaplotlist[[i]] = areaplot # Can also generate one .pdf with all plots
  
  if (INDIVIDUAL.OUTPUTS == TRUE) {
    barplot.filename <- paste0("Bar_", sub("XX", i, OUTPUT.BASENAME))
    ggsave(barplot.filename, barplot, device = OUTPUT.TYPE)
    areaplot.filename <- paste0("Area_", sub("XX", i, OUTPUT.BASENAME))
    ggsave(areaplot.filename, areaplot, device = OUTPUT.TYPE)
    
  }
  print("Individual plots generated")
}

# Output all plots as single .pdf file
if (INDIVIDUAL.OUTPUTS == FALSE) {
  all.bar.plots <- grid.arrange(grobs=barplotlist, ncol=3)
  all.bar.plots <- marrangeGrob(barplotlist, nrow = 4, ncol = 3)
  all.area.plots <- grid.arrange(grobs=areaplotlist, ncol=3)
  all.area.plots <- marrangeGrob(areaplotlist, nrow = 4, ncol = 3)
  barplot.filename <- paste0("Bar_", sub("_XX", "", OUTPUT.BASENAME))
  ggsave(barplot.filename, all.bar.plots, height = 21,width = 21)
  areaplot.filename <- paste0("Area_", sub("_XX", "", OUTPUT.BASENAME))
  ggsave(barplot.filename, all.area,plots, height = 21,width = 21)
  print("Plots sucessfully output as single file")
}

dev.off()
print("Finish BarPlot_Clusters_ByRegion.R")

