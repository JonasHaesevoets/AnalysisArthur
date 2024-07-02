obj = read_rds("../../Desktop/AnalysisArthur/intermediate_data/seurat_obj_central.rds")

filtered_metadata <- obj@meta.data[, !grepl("^RNA_snn", colnames(obj@meta.data))]
obj@meta.data <- filtered_metadata
print(colnames(obj@meta.data))
obj$harmony_clusters_0.22 = NULL


obj$day = obj$orig.ident
obj$orig.ident = NULL
obj$day_mock = obj@meta.data %>% mutate(day_mock = if_else(virus == "Mock", "day0", day))
obj@meta.data$sampletag_name = factor(obj@meta.data$sampletag_name,
                                             levels = c("Mock_r1", "Mock_r2","Mock_r3", "MuHV4_r1", "MuHV4_r2", "MuHV4_r3","del73_r1", "del73_r2", "del73_r3")
)
obj$day_mock_condition = obj@meta.data %>% mutate(day_mock_condition = paste(day_mock, virus, sep = "_"))
obj$virusFactor = factor(obj@meta.data$virusFactor, levels = c("Mock", "MuHV4", "del73"))


obj |> write_rds("../../Desktop/AnalysisArthur/intermediate_data/seurat_obj_central.rds")
