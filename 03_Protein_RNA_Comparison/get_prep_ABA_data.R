# This code inspired by: https://community.brain-map.org/t/whole-mouse-brain-gene-expression-data/447/3

# For reference:
# expression density = sum of expressing pixels / sum of all pixels in division
# expression intensity = sum of expressing pixel intensity / sum of expressing pixels
# expression energy = expression intensity * expression density


library(dplyr)

PARENT_PATH <- getwd()

### FIRST!! run "/home/sg4dm/ZunderLab/project_Protein_RNA/subProject_ABAcomparison/ABA_2021_03_12.ipynb"
### to get section dataset ids needed for this pipeline

aba_summed_for_dotplot <- list()

get_summarized_sectionDataset <- function(section_ds_id, gene_acronym, age) {
  over_11_sa = c('M', 'H', 'D', 'RSP', 'CSP')
  age_tissue_df = data.frame("tissue" = I(list(c('F', 'M', 'H'), over_11_sa, over_11_sa, over_11_sa, over_11_sa)),
                             "age" = c("11", "13", "15", "18", "23"))
  #for (i in 1:nrow(my.genes.sds.ids)) {  
    #i = 1
    #section_ds_id = my.genes.sds.ids$id.1[i] 
    #gene_acronym = my.genes.sds.ids$acronym[i]
    #age = my.genes.sds.ids$age[i]
    
  print(age)
  print(gene_acronym)
  outputs <- RCurl::getURL(glue::glue(
    "http://api.brain-map.org/api/v2/data/SectionDataSet/query.xml?id={section_ds_id}&include=structure_unionizes(structure)"
  )) %>% #(XML::xmlParse) %>% (XML::xmlToDataFrame) %>% return()
    (XML::xmlParse) %>% 
    (XML::xmlToList)
  print("xml fetched")
  data = outputs$`section-data-sets`$`section-data-set`$`structure-unionizes`
  
  expression = data.frame(matrix(unlist(data), nrow=length(data), byrow=T), stringsAsFactors = FALSE)
  
  colnames(expression) <- names(unlist(data[[1]]))
  
  cols_keep <- c("sum-expressing-pixel-intensity.text",
                 "sum-expressing-pixels.text",
                 "sum-pixels.text",
                 "structure.acronym",
                 "structure.name")
  
  sa_ls = age_tissue_df$tissue[which(age_tissue_df$age == age)][[1]]
  expression <- expression[which(expression$structure.acronym %in% sa_ls), cols_keep]
  print("structures selected")
  cols_new_names <- c("sum_expressing_pixel_intensity",
                      "sum_expressing_pixels",
                      "sum_pixels",
                      "structure_acronym",
                      "structure_name")
  
  colnames(expression) <- cols_new_names
  
  expression$genename <- rep(gene_acronym, nrow(expression))
  
  expression$age <- rep(age, nrow(expression))
  
  expression$structure_acronym[which(expression$structure_acronym == '')]
  
  expression$sum_expressing_pixel_intensity = as.numeric(expression$sum_expressing_pixel_intensity)
  expression$sum_expressing_pixels = as.numeric(expression$sum_expressing_pixels)
  expression$sum_pixels = as.numeric(expression$sum_pixels)
  
  if (age == 11) {
    sepi_sum = sum(expression$sum_expressing_pixel_intensity)
    sep_sum = sum(expression$sum_expressing_pixels)
    sp_sum = sum(expression$sum_pixels)
    structure_acronym = "WB"
    structure_name = "whole brain"
    genename = expression$genename[1]
    age = expression$age[1]
    expression<-data.frame(sepi_sum,sep_sum,sp_sum,structure_acronym,structure_name,genename,age)
    names(expression)<-c("sum_expressing_pixel_intensity","sum_expressing_pixels","sum_pixels","structure_acronym","structure_name","genename","age")
  } else {
    expression <- aggregate(expression[,c("sum_expressing_pixel_intensity",
                                          "sum_expressing_pixels",
                                          "sum_pixels")],by=list(structure_acronym=expression$structure_acronym,
                                                                 structure_name=expression$structure_name,
                                                                 genename=expression$genename,
                                                                 age=expression$age),FUN="sum")
    
    sepi_sum = sum(expression$sum_expressing_pixel_intensity[which(expression$structure_acronym %in% c('D', 'RSP'))])
    sep_sum = sum(expression$sum_expressing_pixels[which(expression$structure_acronym %in% c('D', 'RSP'))])
    sp_sum = sum(expression$sum_pixels[which(expression$structure_acronym %in% c('D', 'RSP'))])
    structure_acronym = "Di"
    structure_name = "Diencephalon"
    genename = expression$genename[which(expression$structure_acronym %in% c('D', 'RSP'))][1]
    age = expression$age[which(expression$structure_acronym %in% c('D', 'RSP'))][1]
    de<-data.frame(sepi_sum,sep_sum,sp_sum,structure_acronym,structure_name,genename,age)
    names(de)<-c("sum_expressing_pixel_intensity","sum_expressing_pixels","sum_pixels","structure_acronym","structure_name","genename","age")
    expression <- rbind(expression,de)
    expression <- expression[which(expression$structure_acronym %in% c('Di', 'CSP', 'M', 'H')),]
  }
    
  expression$percent <- as.numeric(expression$sum_expressing_pixels) / as.numeric(expression$sum_pixels)
  
  expression$measurement <- as.numeric(expression$sum_expressing_pixel_intensity) / as.numeric(expression$sum_expressing_pixels)
  print("percent and measurement computed")
  
  expression <- expression[,c("genename", "age", "structure_acronym", "percent", "measurement")]
  
  print("aggregated")
  
  return(expression)
}


# all genes dataset info dowloaded from python script:
all.genes.from.python <- read.csv(paste0(PARENT_PATH,"/out_genes.csv"))

## read in my list of target gene names and ensembl ids
my.genes <- read.csv(paste0(PARENT_PATH,"/ensembleIDs_geneNames_to_entrez.txt"), sep = "\t")

## subset to just genes of interest
my.genes.mm.entrez <- unique(all.genes.from.python[which(all.genes.from.python$acronym %in% stringr::str_to_title(my.genes$Gene.name)),])

## subset to just time points of interest (via reference_space_id)
my.genes.mm.entrez <- my.genes.mm.entrez[which(my.genes.mm.entrez$reference_space_id %in% c(1,2,3,5,6)),]

## translate ref space id to age
ref.space.to.timepoint <- data.frame("ref_space" = c(1,2,3,5,6),
                                     "age" = c(11,13,15,18,23))
my.genes.mm.entrez$age <- ref.space.to.timepoint$age[match(my.genes.mm.entrez$reference_space_id, ref.space.to.timepoint$ref_space)]

## subset cols to only those needed for the function
my.genes.sds.ids <- my.genes.mm.entrez[,c("acronym","id.1","age")]


aba_summed_for_dotplot <- lapply(1:nrow(my.genes.sds.ids), function(x) { #
  get_summarized_sectionDataset(section_ds_id = my.genes.sds.ids$id.1[x], 
                                gene_acronym = my.genes.sds.ids$acronym[x],
                                age = my.genes.sds.ids$age[x])})


aba_summed_for_dotplot_df <- do.call(rbind.data.frame, aba_summed_for_dotplot)
aba_summed_for_dotplot_df <- dplyr::rename(aba_summed_for_dotplot_df, genename = gene_name)

normalit<-function(x){
  #x = x*100
  (x-min(x))/(max(x)-min(x))
}

aba_summed_for_dotplot_df_wb <- aba_summed_for_dotplot_df[which(aba_summed_for_dotplot_df$age == '11'),]
# Subset aba rna to telencephalon ############################################################
aba_summed_for_dotplot_df_F <- aba_summed_for_dotplot_df[which(aba_summed_for_dotplot_df$structure_acronym == 'CSP' & aba_summed_for_dotplot_df$age != '11'),]
aba_summed_for_dotplot_df_F <- rbind(aba_summed_for_dotplot_df_wb, aba_summed_for_dotplot_df_F)
# Normalize within-gene from 0 to 1
aba_summed_for_dotplot_df_F <- aba_summed_for_dotplot_df_F %>%
  dplyr::group_by(genename) %>%
  dplyr::summarise(genename=genename, age=age, percent=percent, measurement = normalit(measurement))
# Write to csv
write.csv(aba_summed_for_dotplot_df_F, paste0(PARENT_PATH,"/aba_rna_telencephalon_summarizedTissueAge.csv"), row.names = F)

# Subset aba rna to diencephalon ############################################################
aba_summed_for_dotplot_df_F <- aba_summed_for_dotplot_df[which(aba_summed_for_dotplot_df$structure_acronym == 'Di' & aba_summed_for_dotplot_df$age != '11'),]
aba_summed_for_dotplot_df_F <- rbind(aba_summed_for_dotplot_df_wb, aba_summed_for_dotplot_df_F)
# Normalize within-gene from 0 to 1
aba_summed_for_dotplot_df_F <- aba_summed_for_dotplot_df_F %>%
  dplyr::group_by(genename) %>%
  dplyr::summarise(genename=genename, age=age, percent=percent, measurement = normalit(measurement))
# Write to csv
write.csv(aba_summed_for_dotplot_df_F, paste0(PARENT_PATH,"aba_rna_diencephalon_summarizedTissueAge.csv"), row.names = F)

# Subset aba rna to mesencephalon ############################################################
aba_summed_for_dotplot_df_M <- aba_summed_for_dotplot_df[which(aba_summed_for_dotplot_df$structure_acronym == 'M' & aba_summed_for_dotplot_df$age != '11'),]
aba_summed_for_dotplot_df_M <- rbind(aba_summed_for_dotplot_df_wb, aba_summed_for_dotplot_df_M)
# Normalize within-gene from 0 to 1
aba_summed_for_dotplot_df_M <- aba_summed_for_dotplot_df_M %>%
  dplyr::group_by(genename) %>%
  dplyr::summarise(genename=genename, age=age, percent=percent, measurement = normalit(measurement))
# Write to csv
write.csv(aba_summed_for_dotplot_df_M, paste0(PARENT_PATH,"/aba_rna_mesencephalon_summarizedTissueAge.csv"), row.names = F)

# Subset aba rna to rhombencephalon ############################################################
aba_summed_for_dotplot_df_H <- aba_summed_for_dotplot_df[which(aba_summed_for_dotplot_df$structure_acronym == 'H' & aba_summed_for_dotplot_df$age != '11'),]
aba_summed_for_dotplot_df_H <- rbind(aba_summed_for_dotplot_df_wb, aba_summed_for_dotplot_df_H)
# Normalize within-gene from 0 to 1
aba_summed_for_dotplot_df_H <- aba_summed_for_dotplot_df_H %>%
  dplyr::group_by(genename) %>%
  dplyr::summarise(genename=genename, age=age, percent=percent, measurement = normalit(measurement))
# Write to csv
write.csv(aba_summed_for_dotplot_df_H, paste0(PARENT_PATH,"/aba_rna_rhombencephalon_summarizedTissueAge.csv"), row.names = F)

