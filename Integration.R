
library(clustree)
library(tidyseurat)
library(readr)
library(Seurat)
set.seed(2023)

# load and merge data files
day_38 <- "intermediate_data/seurat_obj_d38_afterQCdashboard.rds" |> read_rds()
day_8 <- "intermediate_data/seurat_obj_d8_afterQCdashboard.rds" |> read_rds()


obj <- merge(day_8, day_38,add.cell.ids = c("day_8", "day_38"))
day_8 <- NULL
day_38 <- NULL
gc()
obj <- JoinLayers(obj, assay = "RNA")

### convert object to v5 object

meta_data <- obj[[]]

obj.v5 <- CreateSeuratObject(counts = obj[["RNA"]]$counts, meta.data = meta_data)
#obj.v5[[]] <- meta_data
#obj.v5$orig.ident <- meta_data$orig.ident
protein <- CreateAssay5Object(counts = obj[["protein"]]$counts)
sampletags <- CreateAssay5Object(counts = obj[["sampletags"]]$counts)
obj.v5[["protein"]] <- protein
obj.v5[["sampletags"]] <- sampletags

obj <- NULL
gc()

obj.v5[["RNA"]] <- split(obj.v5[["RNA"]], f = obj.v5$orig.ident)
gc()

obj.v5 <- NormalizeData(obj.v5) # individual size factors accoriding to Ahlmann-Eltze et al (2023) could be added here
obj.v5 <- FindVariableFeatures(obj.v5,
                               selection.method = "vst",
                               nfeatures = 2000, 
                               verbose = FALSE)
obj.v5 <- ScaleData(obj.v5)
obj.v5 <- RunPCA(obj.v5)



obj.v5 <- IntegrateLayers(
  object = obj.v5, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE
)
obj.v5 <- JoinLayers(obj.v5, assay= "RNA")
#obj <- JoinLayers(obj, assay = "protein")
#obj <- JoinLayers(obj, assay = "sampletags")

obj.v5 <- FindNeighbors(obj.v5, reduction = "integrated.harmony", dims = 1:40)


# Find clusters using a range of resolutions
resolution.range <- seq(from = 0, to = 1, by = 0.2)

obj.v5 <- Seurat::FindClusters(object = obj.v5, resolution = resolution.range)
resolution.range <- seq(from = 0, to = 0.2, by = 0.02)

obj.v5 <- Seurat::FindClusters(object = obj.v5, resolution = resolution.range)
clustree(obj.v5)
ggsave("clustree_integration.png", path = "../../Desktop/AnalysisArthur/CaseStudy_QC_integration/integration", width = 14, height = 14)

# based on this we decided to put the resolution at 0.18

# obj.v5 <- FindClusters(obj.v5, cluster.name = "harmony_clusters_0.18", resolution = 0.18)
# obj.v5 <- RunUMAP(obj.v5, reduction = "integrated.harmony", reduction.name = "umap_harmony", dims = 1:40)
# obj.v5 |> DimPlot(group.by = "harmony_clusters_0.18", label=T)
# ggsave("umap_0.18.png", path = "../../Documents/machiels_lab_viral/CaseStudy_QC_integration/integration", width = 12, height = 10)
# 
# obj.v5 |>  write_rds("../../Documents/machiels_lab_viral/intermediate_data/seurat_obj_experiment_1_2_integrated.rds")