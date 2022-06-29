#Corey Williams, University of Virginia
#15 July 2019
#Plot colored by expression of markers

print("Start UMAP_R1_Layour_ColoredBy_Marker.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

library(ggfortify)

## Input parameters
CLUSTER_ROUND <- 1
OUTPUT_FOLDER_NAME <- "_Marker_Expression_UMAP_XY"
EXPRS_MAT_INFILE <- "expression_matrix_analysis.csv"
EXPRS_MAT_OTHER <- "expression_matrix_other.csv"
LAYOUT_FILENAME <- "umap_xy_R1_shifted.csv"
PANEL_FILENAME <- "panel.csv"
POINT_SIZE <- 0.0001
PLOT_HEIGHT <- 14
PLOT_WIDTH <- 14

## Read needed files
exprs_mat <- cbind(read.table(EXPRS_MAT_INFILE, header=TRUE, sep=",",
                              check.names=FALSE),
                   read.table(EXPRS_MAT_OTHER, header=TRUE, sep=",",
                              check.names=FALSE))

layout_filename <- LAYOUT_FILENAME
layout_in <- read.table(layout_filename, header=TRUE, sep=",",
                        check.names=FALSE)
panel <- read.table(PANEL_FILENAME, header=TRUE, sep=",", check.names=FALSE)

## Prep dataframe for plotting
plot_vars <- panel[which(panel[,"Plot"]==1), "Fixed_Param"]
plot_df <- as.data.frame(cbind(exprs_mat[,plot_vars], layout_in))

colnames(plot_df)[c(ncol(plot_df)-1, ncol(plot_df))] <- c("umap_x", "umap_y")

## Randomize order for plotting, to avoid cells at the end (often the last file concatenated) going on top of everything else
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)),]
row.names(plot_df) <- NULL #remove old indices, so everything has a new order

## Save plots colored by each marker
#make output folder
time_now <- Sys.time()
output_dir <- paste0(substr(time_now,start=1,stop=10), "_",
                     substr(time_now,start=12,stop=13), ".",
                     substr(time_now,start=15,stop=16), ".",
                     substr(time_now,start=18,stop=19), OUTPUT_FOLDER_NAME)
dir.create(output_dir)
setwd(output_dir)
#loop through variables to plot
for (v in plot_vars){
  ggsave(paste0(v,".png"),
         plot=ggplot(plot_df, aes_string(x="umap_x", y="umap_y", color=paste0("`", v, "`"))) +
           geom_point(size = POINT_SIZE) +
           scale_color_gradientn(colours=c("#0B0688","#18A682","#FFFA00","#FFFB2A"),values=c(0,0.2,0.85,1)) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")),
         height = PLOT_HEIGHT, width = PLOT_WIDTH)
}
setwd("..")

print("Finish UMAP_R1_Layour_ColoredBy_Marker.R")
