#Corey Williams, University of Virginia
#15 July 2019
#Generate violin plot colored by clusters

print("Start Violin_Plot_R2_Clusters.R")

rm(list = ls())
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(ggfortify)
library(ggstance)
library(ggpubr)
library(forcats)
library(dplyr)
library(data.table)

print("Libraries loaded")

## Input parameters
INPUT_FOLDER <- getwd()
CLUSTER_ROUND <- 2
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
PANEL_FILENAME <- "panel.csv"
CLUSTERS_BASENAME <- "cluster_RX_assigns.csv"
MARKER_ORDER_FILENAME <- "marker_order.csv"
PLOT_PARAMETERS_FILENAME <- "Dot_Plot_Parameters.csv" #has to be individually crafted with labels and desired colors
OUTPUT_FILENAME <- "ViolinPlots_R2_Clusters.pdf"
VIOLIN_HEIGHT_FACTOR <- 5
print("Input parameters loaded")

## Read needed files
exprs_mat <- cbind(read.table(EXPRS_MAT_INFILE, header=TRUE, sep=",",
                              check.names=FALSE),
                   read.table(EXPRS_MAT_OTHER, header=TRUE, sep=",",
                              check.names=FALSE))
panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

clusters_filename <- sub("RX", paste0("R", CLUSTER_ROUND), CLUSTERS_BASENAME)
clusters_in <- read.table(clusters_filename, header=TRUE, sep=",",
                          check.names=FALSE)

plot_parameters <- fread(paste0(INPUT_FOLDER, "/", PLOT_PARAMETERS_FILENAME), stringsAsFactors = F)

marker_order <- data.frame(read.csv(paste0(MARKER_ORDER_FILENAME)))

print("Needed files successfully loaded")

## Read in custom marker order, cluster order, and color palette
#plot_vars <- panel[which(panel[,"Plot"]==1), "Fixed_Param"] # Old way, pull in form panel
plot_vars_ordered <- as.vector(marker_order$order) # Marker order
VIOLIN.ORDER <- as.vector(plot_parameters$Cluster) # Cluster order
CUSTOM.PALETTE <- as.vector(plot_parameters$Color) # Custom colors (in order)

print("Custom color palette loaded")

## Prep dataframe for plotting
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars_ordered],clusters_in))
colnames(plot_df)[ncol(plot_df)] <- "Cluster"

print("Dataframe prepared with standard order")

#Apply custom cluster order to dataframe
plot_df_ordered <- plot_df %>%
    mutate(Cluster=factor(Cluster,levels=VIOLIN.ORDER))

print("Custom order applied to dataframe, ready for plotting")

#Make list of violin plots by cluster
plist = sapply(plot_vars_ordered, function(marker_plot) {
  if (marker_plot == plot_vars_ordered[1]){
    ggplot(plot_df_ordered, aes(x = plot_df_ordered[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            legend.position = "none",axis.text.x = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
  else {
    ggplot(plot_df_ordered, aes(x = plot_df_ordered[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            legend.position = "none",axis.text = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
}, simplify=FALSE)

print("List of violin plots by cluster generated")

#save Violin plots
#ggsave(plist)

ggsave("ViolinPlots_Custom_Order_Color.pdf",
       annotate_figure(ggarrange(plotlist = plist,ncol=length(plist)),
                       left = text_grob("Cluster",rot=90)),
       height=max(clusters_in)/VIOLIN_HEIGHT_FACTOR,width=length(plot_vars_ordered))

print("End 12_ViolinPlots_Custom.R")


## To plot only select markers

ALL.MARKERS <- c("Nestin_144Nd","Sox2_160Gd","Pax6_154Sm","Prominin_148Nd",
                    "SSEA-1_157Gd","N-Cadherin_161Dy","CD24_159Tb","Sox1_145Nd",
                    "Tbr2_146Nd","NeuroD1_156Gd","Tbr1_171Yb","Cux1_170Er","Ctip2_168Er",
                    "Dcx_163Dy","TuJ1_89Y","MAP2_164Dy","PSA-NCAM_139La","NeuN_150Nd",
                    "GAD65_162Dy","BLBP_172Yb","GLAST_167Er","GFAP_169Tm","ALDH1A1_165Ho",
                    "A2B5_115In","Olig2_113In","PDGFRa_147Sm","Sox10_151Eu",
                    "OligoO4_153Eu","TrkB_174Yb","p75NTR_176Yb","Ly-6C_166Er",
                    "PECAM1_143Nd","VCAM1_142Nd","MCAM_155Gd","PDGFRb_141Pr","CD45_149Sm",
                    "CD11b_158Gd","F4_80_175Lu","Ki67_152Sm","cl-Casp3_173Yb","Intercalator_191Ir")
NEURON.MARKERS <- c("Pax6_154Sm","Sox2_160Gd","Sox1_145Nd","Nestin_144Nd",
                    "N-Cadherin_161Dy","CD24_159Tb","Prominin_148Nd","SSEA-1_157Gd",
                    "Tbr2_146Nd","NeuroD1_156Gd","Tbr1_171Yb","Cux1_170Er","Ctip2_168Er",
                    "Dcx_163Dy","TuJ1_89Y","MAP2_164Dy","NeuN_150Nd",
                    "GAD65_162Dy","BLBP_172Yb","GLAST_167Er","ALDH1A1_165Ho",
                    "A2B5_115In","Olig2_113In","TrkB_174Yb","p75NTR_176Yb",
                    "Ki67_152Sm","cl-Casp3_173Yb","Intercalator_191Ir", "Cluster")
plot_df_ordered_select <- plot_df_ordered[NEURON.MARKERS]
plot_vars_ordered_select <- NEURON.MARKERS
GLIA.MARKERS <- c("Nestin_144Nd","Sox2_160Gd","Pax6_154Sm","Prominin_148Nd",
                 "N-Cadherin_161Dy","CD24_159Tb","Sox1_145Nd","Cux1_170Er",
                 "Dcx_163Dy","TuJ1_89Y","MAP2_164Dy","GAD65_162Dy","BLBP_172Yb","GLAST_167Er",
                 "GFAP_169Tm","ALDH1A1_165Ho","A2B5_115In","Olig2_113In","PDGFRa_147Sm","Sox10_151Eu",
                 "OligoO4_153Eu","TrkB_174Yb","p75NTR_176Yb","VCAM1_142Nd","MCAM_155Gd",
                 "Ki67_152Sm","cl-Casp3_173Yb","Intercalator_191Ir","Cluster")
plot_df_ordered_select <- plot_df_ordered[GLIA.MARKERS]
plot_vars_ordered_select <- GLIA.MARKERS

HEMA.MARKERS <- c("CD45_149Sm","CD11b_158Gd","F4_80_175Lu","CD24_159Tb",
                  "Sox2_160Gd","Tbr2_146Nd","NeuroD1_156Gd","Tbr1_171Yb","Cux1_170Er","Ctip2_168Er",
                 "Dcx_163Dy","TuJ1_89Y","MAP2_164Dy","NeuN_150Nd","GAD65_162Dy",
                 "BLBP_172Yb","GLAST_167Er","Olig2_113In","OligoO4_153Eu","Ly-6C_166Er","Ki67_152Sm","Cluster")
plot_df_ordered_select <- plot_df_ordered[HEMA.MARKERS]
plot_vars_ordered_select <- HEMA.MARKERS

plist = sapply(plot_vars_ordered_select, function(marker_plot) {
  if (marker_plot == plot_vars_ordered_select[1]){
    ggplot(plot_df_ordered_select, aes(x = plot_df_ordered_select[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            legend.position = "none",axis.text.x = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
  else {
    ggplot(plot_df_ordered_select, aes(x = plot_df_ordered_select[,marker_plot], y = fct_rev(factor(Cluster)), fill = factor(Cluster))) +
      geom_violinh(trim=TRUE,scale = "width") + scale_fill_manual(values = CUSTOM.PALETTE) +
      xlab(marker_plot) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            legend.position = "none",axis.text = element_blank(),axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8))
  }
}, simplify=FALSE)

#save Violin plots
ggsave(OUTPUT_FILENAME,
       annotate_figure(ggarrange(plotlist = plist,ncol=length(plist)),
                       left = text_grob("Cluster",rot=90)),
       height=max(clusters_in)/VIOLIN_HEIGHT_FACTOR,width=length(plot_vars_ordered))

print("Finish Violin_Plot_R2_Clusters.R")
