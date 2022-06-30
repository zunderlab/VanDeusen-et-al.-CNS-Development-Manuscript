##li Zunder, University of Virginia
#2021

# Shift_XY_Shift.R
# 1. Generate plots to illustrate cell groups and how they shift: Shift_Plot_n.png
# 2. Generate plot to illustrate final XY after shifts applied: Shift_Plot_All.png
# 3. Generate CSV with new shifted XY coordinates: umap_xy_R1_shifted.csv
#
#
# Instructions:
# 1. Collect the following files in a single folder:
#      1. Shift_XY_Init.R
#      2. Shift_XY_Run.R
#      3. umap_xy_R1.csv (XY coordinates)
#      4. cluster_R1_assigns.csv (Cluster assignments - this is optional)
# 2. Run previous script, Shift_XY_Init.R
#      1. Outputs Initial_Plot.png and XY_shifts.csv
# 3. View Initial plot.png, and decide the rectangular bounds for cell
#    groups you want to shift, and how much you want to shift them by.
# 4. Enter rectangular bounds and x/y shifts into XY_shifts.csv. Add
#    a new numbered column for each additional group you want to shift.
# 5. After XY_shifts.csv is edited and saved, run this script: Shift_XY_Shift.R
#      1. Outputs Shift_Plot_1.png (through n), Shift_Plot_All.png (colored
#         by cluster if provided, and umap_xy_R1_shifted.csv)
#
#
#

print("Start UMAP_Shift_XY_Run.R")

rm(list = ls(all = TRUE))
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(data.table)

## Input parameters

LAYOUT_FILENAME <- "umap_xy_R1.csv"
CLUST_ASSIGNS_FILENAME <- "cluster_R2_assigns.csv"
SHIFTS_FILENAME <- "XY_shifts.csv"
SHIFT_PLOT_BASENAME <- "Shift_Plot"

PLOT_FMT <- "png" #png or pdf ---note that label background is broken in pdf
PLOT_WTH <- 1000 #1000 for png, 10 for pdf
PLOT_HGT <- 1000 #1000 for png, 10 for pdf
IND_COL <- "red"
#POINT_SIZE <- 20

xy <- read.table(LAYOUT_FILENAME, header=TRUE, sep=",", check.names=FALSE)

size_layout <- nrow(xy)

# Try and calculate the best point size --may want to try and find
# a better function that gives good plots for all datasets.
if (size_layout <= 70) {
  point_size <- 25
} else if (size_layout <= 100000000) {
  point_size <- 20/log10(0.1*size_layout) - 2.5  # meets at 23.66589
} else {
  point_size <- 2.5/log10(0.1*size_layout) # meets at 0.3571429
}

# Try and calculate the best point size --may want to try and find
# a better function that gives good plots for all datasets.
if (size_layout <= 70) {
  point_size <- 25
} else if (size_layout <= 100000000) {
  point_size <- 20/log10(0.1*size_layout) - 2.5  # meets at 23.66589
} else {
  point_size <- 2.5/log10(0.1*size_layout) # meets at 0.3571429
}

shifts <- read.table(SHIFTS_FILENAME, header=TRUE, row.names=1, sep=",", check.names=FALSE)

xy_shifted <- xy

for (i in 1:ncol(shifts)) {
  
  shift_bool <- xy[,1] > shifts["x_min",i] &
    xy[,1] < shifts["x_max",i] &
    xy[,2] > shifts["y_min",i] &
    xy[,2] < shifts["y_max",i]
  
  xy_to_shift <- xy[which(shift_bool),]
  xy_to_not_shift <- xy[which(!shift_bool),]
  
  xy_shifted[which(shift_bool),1] <- xy_shifted[which(shift_bool),1] + shifts["x_shift",i]
  xy_shifted[which(shift_bool),2] <- xy_shifted[which(shift_bool),2] + shifts["y_shift",i]
  
  fwrite(xy_shifted, sub(".csv","_shifted.csv",LAYOUT_FILENAME))
  
  shift_basename <- paste0(SHIFT_PLOT_BASENAME, "_", i)
  
  if (PLOT_FMT == "png") {
    png(paste0(shift_basename,".png"), width=PLOT_WTH*2, height=PLOT_HGT)
  } else if (PLOT_FMT == "pdf") {
    pdf(paste0(shift_basename,".pdf"), width=PLOT_WTH*2, height=PLOT_HGT)
  }
  
  par(mfrow=c(1,2))
  
  plot(x=xy[,1], y=xy[,2], pch=".", cex=point_size, col="grey",
       axes=FALSE, ylab="", xlab="", main="Pre-shift", cex.main=5)
  points(x=xy_to_shift[,1], y=xy_to_shift[,2], pch=".",
         cex=point_size*2, col=IND_COL)
  rect(xleft=shifts["x_min",i], ybottom=shifts["y_min",i],
       xright=shifts["x_max",i], ytop=shifts["y_max",i])
  box()
  
  plot(x=xy_shifted[,1], y=xy_shifted[,2], pch=".", cex=point_size,
       col="grey", axes=FALSE, ylab="", xlab="", main="Post-shift", cex.main=5)
  points(x=xy_shifted[which(shift_bool),1], y=xy_shifted[which(shift_bool),2],
         pch=".", cex=point_size*2, col=IND_COL)
  box()
  dev.off()
  
}

if (!is.null(CLUST_ASSIGNS_FILENAME)) {
  clust_color <- read.table(CLUST_ASSIGNS_FILENAME, header=TRUE, sep=",",check.names=FALSE)[[1]]
} else {
  clust_color <- "grey"
}

if (PLOT_FMT == "png") {
  png("Shift_Plot_All.png", width=PLOT_WTH, height=PLOT_HGT)
} else if (PLOT_FMT == "pdf") {
  pdf("Shift_Plot_All.png", width=PLOT_WTH, height=PLOT_HGT)
}
plot(x=xy_shifted[,1], y=xy_shifted[,2], pch=".", cex=point_size, col=clust_color, axes=FALSE, ylab="", xlab="")
box()
dev.off()

print("Finish UMAP_Shift_XY_Run.R")