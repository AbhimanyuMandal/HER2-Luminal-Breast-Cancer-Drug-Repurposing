library(Seurat)
library(dplyr)
library(ggplot2)

message("==============================================")
message("GSE113196 — NORMAL BREAST RPCA INTEGRATION")
message("==============================================")

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

input_file <- "data/processed/GSE113196_normal_breast_QC_v1.rds"

output_file <- "data/processed/GSE113196_normal_breast_RPCA_v1.rds"

output_dir <- "results/normal_breast"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---------------------------------------------------------
# 1. Load QC-filtered object
# ---------------------------------------------------------

message("")
message("Loading QC-filtered object...")

obj <- readRDS(input_file)

message("Cells: ", ncol(obj))
message("Genes: ", nrow(obj))

if (!"sample_id" %in% colnames(obj@meta.data)) {
  stop("sample_id metadata column not found.")
}

message("")
message("Individuals:")
print(table(obj$sample_id))

# ---------------------------------------------------------
# 2. Split object by individual
# ---------------------------------------------------------

message("")
message("Splitting object by individual...")

obj.list <- SplitObject(
  obj,
  split.by = "sample_id"
)

message("Number of individual objects: ",
        length(obj.list))

# ---------------------------------------------------------
# 3. Normalize each individual
# ---------------------------------------------------------

message("")
message("Normalizing individual datasets...")

obj.list <- lapply(
  obj.list,
  function(x) {
    
    x <- NormalizeData(
      x,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
    
    x <- FindVariableFeatures(
      x,
      selection.method = "vst",
      nfeatures = 3000,
      verbose = FALSE
    )
    
    return(x)
  }
)

message("Normalization and HVG selection complete.")

# ---------------------------------------------------------
# 4. Collect integration features
# ---------------------------------------------------------

message("")
message("Selecting integration features...")

integration_features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 3000
)

message(
  "Integration features: ",
  length(integration_features)
)

# ---------------------------------------------------------
# 5. Scale and PCA each individual
# ---------------------------------------------------------

message("")
message("Running PCA on each individual...")

obj.list <- lapply(
  obj.list,
  function(x) {
    
    x <- ScaleData(
      x,
      features = integration_features,
      verbose = FALSE
    )
    
    x <- RunPCA(
      x,
      features = integration_features,
      npcs = 30,
      verbose = FALSE
    )
    
    return(x)
  }
)

message("Individual PCA complete.")

# ---------------------------------------------------------
# 6. Find RPCA integration anchors
# ---------------------------------------------------------

message("")
message("Finding RPCA integration anchors...")

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = integration_features,
  reduction = "rpca",
  dims = 1:30,
  normalization.method = "LogNormalize"
)

message("RPCA anchors identified.")

# ---------------------------------------------------------
# 7. Integrate datasets
# ---------------------------------------------------------

message("")
message("Integrating individuals using RPCA...")

integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  normalization.method = "LogNormalize"
)

message("Integration complete.")

# ---------------------------------------------------------
# 8. Use integrated assay
# ---------------------------------------------------------

DefaultAssay(integrated) <- "integrated"

message("")
message("Default assay: ",
        DefaultAssay(integrated))

# ---------------------------------------------------------
# 9. Scale integrated data
# ---------------------------------------------------------

message("")
message("Scaling integrated data...")

integrated <- ScaleData(
  integrated,
  verbose = FALSE
)

# ---------------------------------------------------------
# 10. PCA on integrated data
# ---------------------------------------------------------

message("")
message("Running integrated PCA...")

integrated <- RunPCA(
  integrated,
  npcs = 30,
  verbose = FALSE
)

message("Integrated PCA complete.")

# ---------------------------------------------------------
# 11. Integrated UMAP
# ---------------------------------------------------------

message("")
message("Generating integrated UMAP...")

set.seed(1234)

integrated <- RunUMAP(
  integrated,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

message("Integrated UMAP complete.")

# ---------------------------------------------------------
# 12. UMAP by individual
# ---------------------------------------------------------

message("")
message("Generating UMAP by individual...")

p_sample <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "sample_id"
) +
  ggtitle(
    "GSE113196 — RPCA Integrated UMAP by Individual"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    output_dir,
    "UMAP_RPCA_by_individual.png"
  ),
  plot = p_sample,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------------------------------------------------------
# 13. Save integrated object
# ---------------------------------------------------------

message("")
message("Saving integrated object...")

saveRDS(
  integrated,
  output_file
)

# ---------------------------------------------------------
# 14. Final summary
# ---------------------------------------------------------

message("")
message("==============================================")
message("RPCA INTEGRATION COMPLETE")
message("==============================================")

message("")
message("Cells: ", ncol(integrated))
message("Genes: ", nrow(integrated))
message("Individuals: ", length(unique(integrated$sample_id)))

message("")
message("Integrated assay: ",
        DefaultAssay(integrated))

message("")
message("Files saved:")

message(output_file)

message(
  file.path(
    output_dir,
    "UMAP_RPCA_by_individual.png"
  )
)

message("")
message("==============================================")