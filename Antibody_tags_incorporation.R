library(Seurat)
library(tidyverse)
library(tidyseurat)
library(Matrix)
library(data.table)
library(dsb)
set.seed(2023)

# read in all sequencing counts of all Rhapsody wells
# this means also those wells that the rhapsody pipeline does not call as containing a cell

## read in the empty
unfiltered_csv_paths <- c( "../../Desktop/AnalysisArthur/rawData/Output_J8/BD-Analysis-BMachiels-J8_DBEC_MolsPerCell_Unfiltered.csv.gz",
                           "../../Desktop/AnalysisArthur/rawData/Output_J38/BD-Analysis-BMachiels-J38_DBEC_MolsPerCell_Unfiltered.csv.gz")



exp_name <- c("day8", "day38")

for (i in seq_along(unfiltered_csv_paths)) {
  # we read only those proteins who were actually stained
  counts_fread = fread(unfiltered_csv_paths[i], select = c("Cell_Index",
                                                           "CD274|CD274|AHS0004|pAbO",
                                                           "I-A_I-E|H2-Ab_Ad_Aq_Ed_Ek|AMM2019|pAbO"), 
                       showProgress = T)
  #rename the proteins
  colnames(counts_fread) <- c("cell_index","H2_ia_ie_AbSeq","Cd274_AbSeq")
  write_csv(counts_fread,paste0("intermediate_data/",exp_name[i],"unfiltered_prot_counts_fread.csv") )
}

###################


# coming from the specific scripts for each of the data specific scripts which were ran in the QC dashboard
day_8 = read_rds("intermediate_data/seurat_obj_d8_afterQCdashboard.rds")
day_38 = read_rds("intermediate_data/seurat_obj_d8_afterQCdashboard.rds")

setup_chunks <- c(day_8, day_38)

# the counts also containing empty wells for the proteins
unfiltered_prot_counts <- c(
  "intermediate_data/day8unfiltered_prot_counts_fread.csv",
  "intermediate_data/day38unfiltered_prot_counts_fread.csv")

dataset_name <- c("d8", "d38")

# unfiltered_prot_counts <- unfiltered_prot_counts |> map_vec(\(x) paste0("intermediate_data/", x))


plot_list <- list()
# Loop through each element in 'unfiltered_prot_counts'
for (i in seq_along(unfiltered_prot_counts)) {
  # Print the current element
  print(unfiltered_prot_counts[i])
  
  # Read CSV file into a data frame 'all_BD_protein_counts'
  all_BD_protein_counts <- read_csv(unfiltered_prot_counts[i])
  
  # Convert data frame to tidy format and set column names
  protein_counts_cell_is_col <- all_BD_protein_counts |> tidyfst::t_dt() 
  colnames(protein_counts_cell_is_col) <- protein_counts_cell_is_col[1,] 
  protein_counts_cell_is_col <- protein_counts_cell_is_col[-1,]
  
  # Create a tibble of protein counts summed across cells
  protein_sum_tbl <- tibble(cell_sum = protein_counts_cell_is_col |> colSums(),
                            cell = colnames(protein_counts_cell_is_col))
  
  # Set cutoff value
  cutoff <- 6000
  
  # Read setup chunks file
  file = setup_chunks[i]
  
  # Update current Seurat object based on specified conditions
  
  setup_chunks[[i]] <- setup_chunks[[i]] |>
    mutate(
      kept_cell = case_when(
        sampletag_multiplets != "single_hashtag" ~ "no_singlet",
        TRUE ~ "keep"
      )
    )
  
  # Create 'meta_data_protein_sum_tbl' by joining Seurat object data with protein sum table
  setup_chunks[[i]] <- setup_chunks[[i]]
  meta_data_protein_sum_tbl <- setup_chunks[[i]] |>
    as_tibble() |>
    dplyr::rename("cell" = ".cell") |>
    full_join(protein_sum_tbl) |>
    mutate(kept_cell = if_else(
      is.na(kept_cell), "empty_by_BD", kept_cell)
    )
  
  # Clear 'seurat_obj' from memory
  seurat_obj <- NULL
  
  # Perform dsb normalization
  
  # Filter cells to keep
  cells <- meta_data_protein_sum_tbl |> filter(kept_cell == "keep") |> pull(cell)
  cells_tbl <- all_BD_protein_counts |> filter(cell_index %in% cells) 
  cells_matrix <- t(as.matrix(cells_tbl |> select(-cell_index)))
  colnames(cells_matrix) <- cells_tbl |> pull(cell_index)
  
  # Filter empty droplets
  empty_droplets <- meta_data_protein_sum_tbl |> filter(kept_cell == "empty_by_BD") |> pull(cell)
  empty_droplets_tbl <- all_BD_protein_counts |> filter(cell_index %in% empty_droplets) 
  empty_droplets_matrix <- t(as.matrix(empty_droplets_tbl |> select(-cell_index)))
  colnames(empty_droplets_matrix) <- empty_droplets_tbl |> pull(cell_index)
  
  # Perform DSB normalization
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix = cells_matrix, 
    empty_drop_matrix = empty_droplets_matrix, 
    denoise.counts = F, 
    use.isotype.control = F
  )
  
  # Write the result to an RDS file
  write_rds(cells.dsb.norm, paste0("intermediate_data/dsb_matrix_", dataset_name[i], ".rds"))
}


# write_rds(plot_list,paste0(".\\intermediate_data\\protein_reads_QC_plot_list.rds"))





##################
# Define the path to the Seurat object file
path_1 <- "intermediate_data/seurat_obj_integrated.rds"

# Read the Seurat object from the specified path
obj.v5 <- read_rds(path_1)

# Read and preprocess protein count matrices for each experiment
d8 <- "intermediate_data/dsb_matrix_d8.rds" |> read_rds() |> t() |> as_tibble(rownames = "cell") |> mutate(cell = paste0("exp_1_lung_", cell))
d38 <- "intermediate_data/dsb_matrix_d38.rds" |> read_rds() |> t() |> as_tibble(rownames = "cell") |> mutate(cell = paste0("exp_2_lung_", cell))


# Combine all protein count matrices into one
dsb_all <- bind_rows(d8, d38)

# Extract protein features
dsb_features <- colnames(dsb_all)[-1]

# Reduce protein count matrix to cells present in Seurat object
dsb_all <- dsb_all |> filter(cell %in% (obj.v5 |> colnames()))

# Append protein object with "empty" cells with missing cell names
missing_cellnames <- setdiff(colnames(obj.v5), pull(dsb_all, cell))
dsb_all <- dsb_all |> filter(!is.na(cell)) |> bind_rows(tibble(cell = missing_cellnames))

# Match protein count matrix columns with Seurat object cells
dsb_all <- dsb_all[match(obj.v5 |> colnames(), dsb_all$cell),]

# Convert protein count matrix to sparse matrix
dsb_all_t <- dsb_all |> dplyr::select(-cell) |> t() |> as.sparse()
colnames(dsb_all_t) <- pull(dsb_all, cell)

# Create an assay object from the sparse matrix
dsb_all_seur <- CreateAssay5Object(counts = dsb_all_t)

# Define protein markers
proteins <- c(
  "i-a-i-e-h2-ab-ad-aq-ed-ek-amm2019-p-ab-o" = "H2-ia-ie-AbSeq",
  "cd274-cd274-amm2038-p-ab-o" = "Cd274-AbSeq"
)

# Subset Seurat object protein counts based on defined protein markers
adt_counts <- obj.v5@assays$protein$counts[names(proteins),]
rownames(adt_counts) <- names(proteins)
# Create an assay object from the subsetted protein counts and sparse matrix
dsb_all_seur <- CreateAssay5Object(counts = adt_counts, data = dsb_all_t)

# Set Seurat object assay version to v5
options(Seurat.object.assay.version = "v5")

# Add assay object to Seurat object as "adt" assay
obj.v5[["adt"]] <- dsb_all_seur

# Change protein assay to v5 assay
DefaultAssay(obj.v5) <- "adt"

# Join layers in Seurat object
obj.v5 <- JoinLayers(obj.v5)

# Process and format protein data
dsb <- obj.v5@assays$adt$data |> t() |>
  as_tibble(rownames = ".cell") |>
  pivot_longer(cols = "H2-ia-ie-AbSeq":"Cd274-AbSeq", names_to = "marker") |>
  mutate(method = "prot_counts") |>
  left_join(obj.v5@meta.data |> as_tibble(rownames = ".cell") |> dplyr::select(.cell, orig.ident)) |>
  mutate(dsb_zero = abs(value)) |>
  group_by(marker, orig.ident) |>
  mutate(
    method_marker_0.9_quantile_dsb_zero = quantile(dsb_zero, probs = 0.9, na.rm = TRUE),
    method_marker_max_dsb_zero = quantile(dsb_zero, probs = 1, na.rm = TRUE),
    quantile_0.9_scaled_dsb_zero = dsb_zero / method_marker_0.9_quantile_dsb_zero,
    quantile_max_scaled_dsb_zero = dsb_zero / method_marker_max_dsb_zero
  ) |>
  ungroup()

# Pivot formatted protein data to wide format
scaled_data_tbl <- dsb |> dplyr::select(.cell, marker, quantile_max_scaled_dsb_zero) |>
  pivot_wider(names_from = marker, values_from = quantile_max_scaled_dsb_zero)

# Arrange protein data according to Seurat object cells
scaled_data_mtx <- scaled_data_tbl |> dplyr::select(-.cell) |> as.matrix()
rownames(scaled_data_mtx) <- pull(scaled_data_tbl, .cell)
scaled_data_mtx <- scaled_data_mtx[match(obj.v5 |> colnames(), scaled_data_tbl$.cell),]
scaled_data_mtx_t <- scaled_data_mtx |> t()

# Add scaled protein data to Seurat object as "protein" assay
obj.v5@assays$adt$scale.data <- scaled_data_mtx_t
obj.v5[["protein"]] <- NULL
new_rownames <- c("H2-ia-ie-AbSeq", "Cd274-AbSeq")
rownames(obj.v5@assays$adt$counts) <- new_rownames
print(rownames(obj.v5@assays$adt$counts))


# Write the modified Seurat object to a file
write_rds(x = obj.v5, file = "intermediate_data/seurat_obj_central.rds")

protein_names <- rownames(obj.v5@assays$protein)
protein_names
rownames(obj.v5@assays$adt$counts)