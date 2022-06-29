#Zunder Lab, University of Virginia
#Add barcode-negative (bc_neg) and barcode-positive (bc_pos) parameters
#and output .fcs file containing these new parameters for enhanced barcode stringency
#Further details in Fread, et al. 2017 (https://doi.org/10.1142/9789813207813_0054)

print("Start 01_BC_Extra.R")

rm(list = ls(all = TRUE))
.libPaths(c( .libPaths(), "~/R/4.1.1"))

bc_keys <- c("102Pd", "104Pd", "105Pd", "106Pd", "108Pd", "110Pd")

cat("/n")

for (f in list.files(pattern=".fcs$")) {
  
  cat("Calculating BC-negative sums for ", f, "\n", sep="")
  
  input_frame <- read.FCS(f, transformation=FALSE)
  input_exprs <- exprs(input_frame)
  marker_names <- pData(parameters(input_frame))[,"desc"]
  marker_names[which(is.na(marker_names))] <- pData(parameters(input_frame))[which(is.na(marker_names)),"name"]
  colnames(input_exprs) <- marker_names
  
  bc_markers <- c()
  for (k in bc_keys) {
    bc_markers <- c(bc_markers, marker_names[grep(k,marker_names)])
  }
  
  bc_only <- input_exprs[,bc_markers]
  bc_only <- asinh(bc_only/5)
  #calculate the sum of 3-max and 3-min BC measurements for each cell, output as 2-column matrix
  new_cols <- t(apply(bc_only, 1, function(x) {
                                            y=sort(x)
                                            return(c(sum(y[1:3]), sum(y[4:6])))
                                            }))
  colnames(new_cols) <- c("bc_neg", "bc_pos")
  new_cols <- sinh(new_cols)*5
  
  #add 2 new columns to fcs file, and write out in the same direrctory with modified filename
  output_frame <- fr_append_cols(input_frame, new_cols)
  write.FCS(x=output_frame, filename=sub(".fcs","-bc_extra.fcs", f))
}

print("Finish 01_BC_Extra.R")
