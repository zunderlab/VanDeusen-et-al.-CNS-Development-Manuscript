#Eli Zunder, University of Virginia
#2021
# Shift_XY_Init.R
# 1. Generate the plot that will be used to determine rectangular bounds for
#    cells to shift, and how much to shift them by.
# 2. Output CSV file that can be used to enter rectangular bounds and shifts.
#
#
# Instructions:
# 1. Collect the following files in a single folder:
#      1. Shift_XY_Init.R (This script)
#      2. Shift_XY_Shift.R (To be run next)
#      3. umap_xy_R1.csv (XY coordinates)
#      4. cluster_R1_assigns.csv (Cluster assignments - this is optional)
# 2. Run this script, Shift_XY_Init.R
#      1. Outputs Initial_Plot.png and XY_shifts.csv
# 3. View Initial plot.png, and decide the rectangular bounds for cell
#    groups you want to shift, and how much you want to shift them by.
# 4. Enter rectangular bounds and x/y shifts into XY_shifts.csv. Add
#    a new numbered column for each additional group you want to shift.
# 5. After XY_shifts.csv is edited and saved, run Shift_XY_Shift.R
#      1. Outputs Shift_Plot_1.png (through n), Shift_Plot_All.png (colored
#         by cluster if provided, and umap_xy_R1_shifted.csv)
#
#
#


print("Start UMAP_Shift_XY_Initialize.R")

rm(list = ls(all = TRUE))
.libPaths(c(.libPaths(), "~/R/4.1.1"))

library(data.table)

## Input parameters
LAYOUT_FILENAME <- "umap_xy_R1.csv"
SHIFTS_FILENAME <- "XY_shifts.csv"
INIT_PLOT_BASENAME <- "Initial_Plot_"
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

# Write out the plot as PNG or PDF
if (PLOT_FMT == "png") {
  png(paste0(INIT_PLOT_BASENAME,".png"), width=PLOT_WTH, height=PLOT_HGT)
} else if (PLOT_FMT == "pdf") {
  pdf(paste0(INIT_PLOT_BASENAME,".pdf"), width=PLOT_WTH, height=PLOT_HGT)
}
x_ticks <- floor(min(xy[,1])):ceiling(max(xy[,1]))
y_ticks <- floor(min(xy[,2])):ceiling(max(xy[,2]))
plot(x=xy[,1], y=xy[,2], pch=".", cex=point_size, col="grey", xlab="", ylab="")
axis(side=1, at=x_ticks)
axis(side=2, at=y_ticks)
abline(h=y_ticks, v=x_ticks, col="lightgray", lty=3)
dev.off()

shift_csv_template <- matrix(data='', nrow=6, ncol=1, dimnames=list(c("x_min","x_max","y_min","y_max","x_shift","y_shift"), c("1")))

write.csv(shift_csv_template, file=SHIFTS_FILENAME, row.names=TRUE)

print("Finish UMAP_Shift_XY_Initialize.R")
