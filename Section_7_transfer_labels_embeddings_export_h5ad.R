# Call from conda-environment MegaBundle
" Transfer the sub_sub_trajectory labels from reference to a query dataset (here from pooled WT to Sox9 Regulatory KO). Then save the new RDS file with filename appended with _subsubtrajectory

Usage: transfer_labels.R --file_sc_obj_ref=<file> --file_sc_obj_query=<file> --file_functions=<file>
Options:
  -h --help			Show this screen.
  --file_sc_obj_ref=<file>		The reference file to transfer labels from
  --file_sc_obj_query=<file>		The query file to transfer labels from
  --file_functions=<file>   The .R file containing all functions
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(tidyr))
suppressMessages(library(docopt))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj_ref <- arguments$file_sc_obj_ref
file_sc_obj_query <- arguments$file_sc_obj_query
file_functions <- arguments$file_functions

message("file_sc_obj_ref: ", file_sc_obj_ref)
message("file_sc_obj_query: ", file_sc_obj_query)
message("file_functions: ", file_functions)

# --- load all functions defined in the r-script
source(file_functions)

# --- Custom settings
subset_key <- "sub_trajectory"
subset_names <- c("Limb mesenchyme trajectory")
parameterset <- "mine"
assay="RNA"

# --- Read files
sc_obj_ref <- readRDS(file_sc_obj_ref)
sc_obj <- readRDS(file_sc_obj_query) # use simpler name

# Filter cells that fit the subsetting criteria defined at the start
cells_to_keep <- sc_obj@meta.data %>%
  filter(if_any(subset_key) %in% subset_names) %>%
  rownames()
sc_obj_subset <- sc_obj %>% subset(cells=cells_to_keep)

if(parameterset=="CX"){
  # specific parameter sets (my parameters)
    nfeatures <- 2000
  if(ncol(sc_obj_ref) > 50000){
    npcs = 50
    } else if (ncol(sc_obj_ref) > 1000) {
    npcs = 30
    } else {
    npcs = 10
    }
} else if(parameterset=="mine"){
  # specific parameter sets (my parameters)
  nfeatures <- 2000
  npcs <- 50
}

# --- Find anchors and transfer labels
anchors <- FindTransferAnchors(reference=sc_obj_ref, query=sc_obj_subset, reference.reduction="pca", dims=1:npcs)
predictions <- TransferData(anchorset=anchors, refdata=sc_obj_ref$sub_sub_trajectory)
sc_obj_subset <- AddMetaData(sc_obj_subset, metadata=predictions[, "predicted.id"], col.name="sub_sub_trajectory")

# Get percentile of the score=0.8
ecdf_prediction <- ecdf(predictions$prediction.score.max)
message("The percentile of 0.8 max prediction score is: ", ecdf_prediction(0.8))

# Transfer the labels of the WT cells from the reference dataset to the query dataset
common_wt_cells <- intersect(Cells(sc_obj_ref), Cells(sc_obj_subset))
sc_obj_subset@meta.data[common_wt_cells, "sub_sub_trajectory"] <- sc_obj_ref@meta.data[common_wt_cells, "sub_sub_trajectory"]

# --- transfer umap embeddings from refernce to the query
sc_obj_subset <- MapQuery(anchorset=anchors,
  query=sc_obj_subset,
  reference=sc_obj_ref,
  reference.reduction="pca",
  reduction.model="umap",
  integrateembeddings.args=list(dims.to.integrate=1:npcs),
  query.dims=1:npcs,
  reference.dims=1:npcs
  )

# Save the Sox9 Limb mesenchyme object as Seurat object (rds) and as anndata object (h5ad) for scVelo figures
file_rds <- "Sox9_Regulatroy_KO_WT_LimbMes_subsubtrajectory.rds"
saveRDS(sc_obj_subset, file=file_rds)
file_h5ad <- gsub(file_rds, pattern=".rds", replacement=".h5ad")
seurat_to_anndata(sc_obj=sc_obj_subset,
  assay=assay,
  cleanup_genenames=FALSE,
  file_h5ad=file_h5ad,
  reduction_umap="ref.umap",
  reduction_pca="ref.pca")

# Save the WT Limb mesenchyme object as anndataobject (h5ad) for scVelo figures
file_h5ad <- gsub(file_sc_obj_ref, pattern=".rds", replacement=".h5ad")
seurat_to_anndata(sc_obj=sc_obj_ref,
  assay=assay,
  cleanup_genenames=FALSE,
  file_h5ad=file_h5ad)

# --- Transfer these labels into the full query dataset (not Limb Mesenchyme alone)
sc_obj <- AddMetaData(sc_obj, metadata=sc_obj_subset$sub_sub_trajectory, col.name="sub_sub_trajectory")
# Transfer sub_trajectory names to sub_sub_trajectory, when sub_sub_trajectory is NA
sc_obj@meta.data <- sc_obj@meta.data %>% mutate(sub_sub_trajectory = case_when(
  is.na(sub_sub_trajectory) ~ sub_trajectory,
  TRUE ~ sub_sub_trajectory
))

saveRDS(sc_obj, file="Sox9_Regulatroy_KO_WT_subsubtrajectory.rds")
