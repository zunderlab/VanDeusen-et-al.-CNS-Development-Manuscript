#Adapted by Amy Van Deusen from script written by Corey Williams, University of Virginia
#2022
#Generate box plots by sample type (outputs both to see which more aesthetically pleasing)

print("Start ClusterReplicateAnalysis.R")

rm(list = ls(all = TRUE))
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(janitor)
library(gridExtra)

print("Libraries loaded")

## Input parameters ===============================================================================
INPUT.FOLDER <- getwd()
OUTPUT.FOLDER <- INPUT.FOLDER
METADATA.FILENAME <- "metadata_master.csv" # Metadata for combined dataset
ABUNDANCES.FILENAME <- "CellsPerCluster_File.csv" # Should be organized accodring to METADATA.TYPE
PLOT.PARAMETER.FILENAME <- "Scatter_Plot_Parameters.csv"
CELL.CUTOFF <- 5 # Clusters with less cells than this value will be excluded from variance calculations/heatmap, NULL = do not apply threshold              
THRESHOLD.METHOD <- "individual" # "cluster" = apply threshold to total cells in cluster for each replicate
                              # "indivdal" = apply threshold to number of cells in individual replicates
                              # NULL = do not apply threshold
OUTPUT.BASENAME <- "Plot_Cluster_Abundances_XX.pdf" #XX Must be used in name for outputs!
OUTPUT.TYPE <- "pdf" # choose png or pdf
CUSTOM.COLORS <- c("E11B"="#F242F5", "E12B"="#F242F5", "E13C"="#4E42F5", "E14C" = "#4E42F5",
                  "E15C"="#4E42F5", "E16C"="#4E42F5", "E17C"="#4E42F5", "E18C"="#4E42F5",
                  "P0C"="#4E42F5", "P1C"="#4E42F5", "P2C"="#4E42F5", "P3C"="#4E42F5",
                  "P4C"="#4E42F5", "E13D"="#18C92D", "E14D"="#18C92D", "E15D"="#18C92D",
                  "E16D"="#18C92D", "E17D"="#18C92D", "E18D"="#18C92D", "P0D"="#18C92D",
                  "P1D"="#18C92D","P2D"="#18C92D","P3D"="#18C92D","P4D"="#18C92D",
                  "E13M"="#FFC316","E14M"="#FFC316","E15M"="#FFC316","E16M"="#FFC316",
                  "E17M"="#FFC316","E18M"="#FFC316","P0M"="#FFC316","P1M"="#FFC316",
                  "P2M"="#FFC316","P3M"="#FFC316","P4M"="#FFC316","E13H"="#FF4E41",
                  "E14H"="#FF4E41","E15H"="#FF4E41","E16H"="#FF4E41","E17H"="#FF4E41",
                  "E18H"="#FF4E41","P0H"="#FF4E41","P1H"="#FF4E41","P2H"="#FF4E41",
                  "P3H"="#FF4E41","P4H"="#FF4E41")

## Read in necessary files =====================================================
metadata <- read.csv(paste0(INPUT.FOLDER,"/",METADATA.FILENAME))
abundances <- read.table(ABUNDANCES.FILENAME, header=TRUE, sep=",", check.names=FALSE)
plot.parameters <- fread(paste0(INPUT.FOLDER, "/", PLOT.PARAMETER.FILENAME), stringsAsFactors = F)

## Generate plotting dataframe =================================================
# Extract necessary metadata
metadata.in <- metadata[c("Total_Cells","Orig_File","Sample_Type")]

# Reformat abundances
abundances.in <- t(abundances)
abundances.in <- row_to_names(abundances.in, row_number = 1)
abundances.in.sorted <- abundances.in[order(as.numeric(row.names(abundances.in))), ]

# Combine metadata and cell abundances
plot.df.init <- cbind(metadata.in, abundances.in.sorted)

## Generate box plot of all clusters for each sample type =======================================================
# Initialize plots and set custom color palette
scatterplotlist = list() # Generate empty list in which to collect scatter plots
varlist = list() # Generate empty list in which to collect values to generate variance plots
varplotlist = list() # Generate empty list in which to collect variance plots

#Cycle through each cluster and build a line graph plot for each metadata type

n.clusters <- ncol(plot.df.init) - 3
n.cluster.cols <- ncol(plot.df.init)

for (i in 4:n.cluster.cols) {
  # Prepare plot dataframe
  cluster.id <- i - 3
  cluster.name <- plot.parameters[cluster.id,2]
  plot.df.i <- cbind(plot.df.init[1:3], plot.df.init[i])
  plot.df.i.pct <- mutate(plot.df.i, Abundance = (plot.df.i[,4]/plot.df.i[,1]*100)) # Calculate abundances (as percentage)
  plot.df <- as.data.frame(cbind(plot.df.i.pct[,3], plot.df.i.pct[,5])) # Take sample types and abundances only
  colnames(plot.df) <- c("Sample_Type", "Abundance")
  plot.df$Abundance <- as.numeric(plot.df$Abundance) # Convert abundance values to numeric
  
  # Set plot parameters
  plot.filename <- sub("XX", i, OUTPUT.BASENAME)
  
  # Generate scatter plots
  scatterplot <- ggplot(plot.df, aes(x=Sample_Type, y=Abundance, color = Sample_Type)) + 
    geom_point(shape=1, size=1, stroke=0.5) +
    scale_color_manual(values = CUSTOM.COLORS) + 
    scale_y_continuous(labels = function(x) format(x, nsmall = 2)) +
    theme_classic() + labs(title = cluster.name) +
    guides(shape = guide_legend(override.aes = list(size = 0.1))) +
    guides(fill = guide_legend(ncol=1)) +
    theme(legend.position="none", plot.title = element_text(size = 5.5),
          axis.text.x = element_text(size = 2.3, angle = 45, hjust=1, vjust=1.3), axis.title.x = element_blank(), 
          axis.text.y = element_text(size = 6, hjust=1.5), axis.title.y = element_blank()) 
    
  scatterplotlist[[i]] = scatterplot
  
  #Calculate variance in cluster abundance for each sample type
  sample.types <- unique(plot.df$Sample_Type)
  for (s in sample.types) {
    plot.df.var.init <- filter(plot.df.i.pct[3:5], Sample_Type == s)
    colnames(plot.df.var.init) <- c("Sample_Type","Cluster","Abundance")
    # Optional: Apply threshold for low-abudance clusters to remove from analysis
    if (!is.null(CELL.CUTOFF) & THRESHOLD.METHOD == "cluster") {
      # Apply threshold to whole cluster (i.e. total cells in cluster for all sample replicates < CELL.CUTOFF, whole cluster not analyzed for variance)
      if(sum(plot.df.var.init$Cluster) <= CELL.CUTOFF) {
        plot.df.var.init$Abundance <- NA # If 1 or more subthreshold clusters, throw out whole sample from variance analysis
      }
    }  
    if (!is.null(CELL.CUTOFF) & THRESHOLD.METHOD == "individual") {
      # Apply threshold to each cluster (i.e. if any cluster replicate < CELL.CUTOFF, whole cluster not analyzed for variance)
      plot.df.var.init$Cluster[plot.df.var.init$Cluster <= CELL.CUTOFF] <- NA # Identify subthreshold clusters and change to NA
      if (sum(is.na(plot.df.var.init$Cluster)) >= 1) {
        plot.df.var.init$Abundance <- NA # If 1 or more subthreshold clusters, throw out whole sample from variance analysis
      }
    }
    plot.df.var.s <- mutate(plot.df.var.init, Variance = var(plot.df.var.init[,2])) # Calculate variances (as percentage)
    varlist[[s]] <- c(s,cluster.id,plot.df.var.s[1,3])
  }
  varplot.df <- t(as.data.frame(varlist))
  colnames(varplot.df) <- c("Sample_Type","Cluster","Variance")
  rownames(varplot.df) <- NULL
  varplotlist[[i]] = varplot.df
}

# Output all scatter plots as single .pdf file
PLOT.ORDER <- plot.parameters$Dot_Plot_Name
scatterplotlist <- scatterplotlist[-c(1:3)] # Remove empty list components
names(scatterplotlist) <- sort(plot.parameters$Dot_Plot_Name)
scatterplotlist.ordered <- scatterplotlist[as.character(PLOT.ORDER)]
all.scatter.plots <- grid.arrange(grobs = scatterplotlist, ncol = 2, top="All Clusters")
scatterplot.filename <- paste0("ScatterPlot_", sub("_XX", "", OUTPUT.BASENAME))
#ggsave(scatterplot.filename, all.scatter.plots, height=24, width=12) #for All Plots
ggsave(scatterplot.filename, all.scatter.plots, height=48, width=5) #for Enlarged Plots
print("Scatter plots sucessfully output as single file")


## Output heatmap of cluster variance for all sample types ==============================================
# Generate heatmap dataframe
varplotlist <- varplotlist[-c(1:3)] # Remove empty list components
varplot.final <- as.data.frame(do.call(rbind, varplotlist)) # Unlist variances calculated for each cluster and sample type
varplot.final$Variance <- as.numeric(varplot.final$Variance)
varplot.final <- mutate(varplot.final, LogVar = log(Variance)) # Convert variance to numeric for continuous scale in ggplot
varplot.final$Cluster <- factor(varplot.final$Cluster, levels = rev(c(plot.parameters$Cluster))) # Apply custom ordering to clusters

#Generate and output heatmap with cluster variances per sample type
varplot <- ggplot(varplot.final, aes(x = Sample_Type, y = Cluster, fill = LogVar)) +
  geom_tile() + coord_fixed() +
  scale_fill_gradient2(low = "white", high = "red", midpoint = -1, na.value = "#F0F0F0") +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust=1), axis.title.x = element_blank())
varplot.filename <- paste0("VarianceHeatmap_", sub("_XX", "", OUTPUT.BASENAME))
ggsave(varplot.filename, varplot, height=21, width=21)
print("Heatmap sucessfully output as single file")

print("Finish ClusterReplicateAnalysis.R")

