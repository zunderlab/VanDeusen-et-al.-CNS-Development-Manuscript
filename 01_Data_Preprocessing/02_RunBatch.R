#Zunder Lab, University of Virginia
#Script to input parameters for batch adjustment by 03_BatchAdjust.R
#Further details in Schuyler, et al. 2019 (https://10.3389/fimmu.2019.02367)

print("Start 02_RunBatch.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

INPUT.DIR <- "" #directory where FCS files are
print(INPUT.DIR)
source("03_BatchAdjust.R")

BatchAdjust(basedir=".", outdir=".", channelsFile = "ChannelsToAdjust.txt", batchKeyword="Set", anchorKeyword="Universal",
            method="65p",transformation=FALSE,addExt="_65p",plotDiagnostics=TRUE)

print("Finish 02_RunBatch.R")