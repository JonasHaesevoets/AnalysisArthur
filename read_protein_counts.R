

library(data.table)
library(R.utils)
library(readr)
unfiltered_csv_paths <- c( "../../Desktop/AnalysisArthur/rawData/Output_J8/BD-Analysis-BMachiels-J8_DBEC_MolsPerCell_Unfiltered.csv.gz",
                           "../../Desktop/AnalysisArthur/rawData/Output_J38/BD-Analysis-BMachiels-J38_DBEC_MolsPerCell_Unfiltered.csv")



exp_name <- c("day8", "day38")

for (i in seq_along(unfiltered_csv_paths)) {
  counts_fread = fread(unfiltered_csv_paths[i], select = c("Cell_Index",
                                                           "I-A_I-E|H2-Ab_Ad_Aq_Ed_Ek|AMM2019|pAbO",
                                                           "CD274|Cd274|AMM2038|pAbO"), 
                       showProgress = T)
  
  colnames(counts_fread) <- c("cell_index","H2_ia_ie_AbSeq","Cd274_AbSeq")
  write_csv(counts_fread,paste0("../../Desktop/AnalysisArthur/intermediate_data/",exp_name[i],"unfiltered_prot_counts_fread.csv") )
  
}