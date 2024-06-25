# This script loads BD Rhapsody pipeline multimodal count data (transcripts, AbSeq, cell hashtags)
# into separate Seurat files per library and also combines Seurat files per experiment.
# The Seurat objects are stored in the "intermediate_data" folder as .rds files.

# Set a seed for reproducibility
set.seed(2023)

# Install necessary packages
install.packages(c("Seurat", "dplyr", "janitor", "forcats", "tidyseurat", 
                   "Matrix", "vroom", "tidyfst", "readr", "stringr", "tidyverse"))
install.packages("easypackages")

# Load installed packages
easypackages::libraries("Seurat", "dplyr", "janitor", "forcats", "tidyseurat", 
                        "Matrix", "vroom", "tidyfst", "readr", "stringr", "tidyverse")

# Define the root directory for raw data
raw_data_root <- "raw_data/machiels_lab/viral"

# Define input and output file names
# Avoid using "_" for experiment names because Seurat can't handle it
file_names_tbl <- tribble(
  ~name_for_merged_seurat_file, ~raw_data_root, ~dbec_file_path, ~name_for_seurat_file, ~sample_tag_reads_per_cell,
  
  # Experiment 1
  "day8", # Name for merged Seurat file
  "../../Desktop/AnalysisArthur/rawData/", # Absolute path for raw_data_root
  "output_J8/Combined_BD-Analysis-BMachiels-J8_DBEC_MolsPerCell.csv", # Path to dbec file with relevant count data
  "seurat_obj_d8_raw_dbec", # Name for intermediate Seurat file output
  "output_J8/BD-Analysis-BMachiels-J8_Sample_Tag_ReadsPerCell.csv", # Path to sample tag reads per cell file
  
  # Experiment 2
  "day38", # Name for merged Seurat file
  "../../Desktop/AnalysisArthur/rawData/", # Absolute path for raw_data_root
  "output_J38/Combined_BD-Analysis-BMachiels-J38_DBEC_MolsPerCell.csv", # Path to dbec file with relevant count data
  "seurat_obj_d38_raw_dbec", # Name for intermediate Seurat file output
  "output_J38/BD-Analysis-BMachiels-J38_Sample_Tag_ReadsPerCell.csv" # Path to sample tag reads per cell file
)

# Define output paths
file_path <- vector("list")
file_path$output <- "../../Desktop/AnalysisArthur/output/"
file_path$intermediate_data <- "../../Desktop/AnalysisArthur/intermediate_data/"

# Function to read Rhapsody multi-assay CSV file into Seurat object
read_rhapsody_multi_assay_tibble_based <- function(dbec_counts_path, project_name, sample_tag_reads) {
  # Load the dbec corrected count data using vroom
  counts <- vroom(dbec_counts_path, skip = 7)
  
  # Inspect column names to ensure 'Cell_Index' exists
  print("Columns in counts data frame:")
  print(colnames(counts))
  
  if(!"Cell_Index" %in% colnames(counts)) {
    stop("The 'Cell_Index' column was not found in the counts data frame.")
  }
  
  barcodes <- counts %>% pull(Cell_Index)
  
  print("Count table loaded")
  
  # Extract and clean protein counts
  protein <- counts[grep(names(counts),pattern = "pAbO")] |> 
    clean_names() |>  #get rid of special characters that would interfere with downstream functions in feature names
    tidyfst::t_dt() #efficiently transposes dataframes 
  colnames(protein) <- barcodes # reapply cell index to put the column names back in the object
  
  # Extract transcriptome counts (everything except protein counts)
  transcriptome <- counts %>% select(-starts_with("pAbO"), -Cell_Index)
  print("Transposing transcriptome tibble")
  transcriptome <- t(as.matrix(transcriptome))
  colnames(transcriptome) <- barcodes
  rm(counts) # Free up memory
  gc() # Free up memory
  
  print("Creating Seurat object")
  seurat <- CreateSeuratObject(counts = transcriptome, project = project_name)
  
  print("Adding assay objects")
  seurat[['protein']] <- CreateAssayObject(counts = protein)
  
  print("Adding sample tags")
  sample_tag_reads <- read_csv(sample_tag_reads, skip = 7) %>%
    mutate(Cell_Index = as.character(Cell_Index)) %>%
    clean_names()
  
  # colnames(sample_tag_reads) <- str_remove_all(colnames(sample_tag_reads), "_mm_st_ab_o")
  # barcodes <- pull(sample_tag_reads, Cell_Index)
  # sample_tag_reads <- sample_tag_reads %>% select(-Cell_Index)
  # sample_tag_reads <- t(as.matrix(sample_tag_reads))
  # colnames(sample_tag_reads) <- barcodes
  # seurat[['sampletags']] <- CreateAssayObject(counts = sample_tag_reads)
  # 
  # 
  # sample_tag_reads <- read_csv(sample_tag_reads, skip = 7) |> ## skip 7 because first 7 lines are not the data themselves
  #   mutate(Cell_Index=as.character(Cell_Index)) |> clean_names()
  # 
  colnames(sample_tag_reads) <- str_remove_all(colnames(sample_tag_reads), "_mm_st_ab_o")
  barcodes <- pull(sample_tag_reads,cell_index)
  sample_tag_reads <- sample_tag_reads |> select(-cell_index) |> tidyfst::t_dt()
  colnames(sample_tag_reads) <- barcodes
  seurat[['sampletags']] <- CreateAssayObject(counts = sample_tag_reads)
  return(seurat)
}

# Main script

# Iterate through rows of "file_names_tbl" to get input files and store the output data files
for (line in 1:nrow(file_names_tbl)) {
  # Define absolute paths
  seurat_path_absolute <- paste0(file_names_tbl$raw_data_root[line], file_names_tbl$dbec_file_path[line])
  proj_name <- file_names_tbl$name_for_merged_seurat_file[line]
  sample_tag_reads <- paste0(file_names_tbl$raw_data_root[line], file_names_tbl$sample_tag_reads_per_cell[line])
  
  print(paste("Processing Seurat object for:", proj_name))
  
  seurat_obj <- read_rhapsody_multi_assay_tibble_based(dbec_counts_path = seurat_path_absolute, # Absolute path
                                                       project_name = proj_name,
                                                       sample_tag_reads = sample_tag_reads)
  
  # Save the individual Seurat object
  seurat_obj_path <- paste0(file_path$intermediate_data, file_names_tbl$name_for_seurat_file[line], ".rds")
  write_rds(seurat_obj, file = seurat_obj_path)
  
  gc() # Free up memory
}
