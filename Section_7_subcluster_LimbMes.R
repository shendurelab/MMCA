# Call from conda-environment MegaBundle
"Use this code to find DE genes and annotate clusters within the Limb Mesenchymal trajectory on the allWT dataset

Usage: subcluster_LimbMes.R --file_sc_obj_WT=<file>  --res=<value> --file_refbiomart=<file>
Options:
  -h --help			Show this screen.
  --file_sc_obj_WT=<file>		The file containing RDS file of WT cells to subcluster
  --res=<value>   The resolution to subcluster
" -> doc

# --- Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the input parameters
file_sc_obj_WT <- arguments$file_sc_obj_WT
res <- as.numeric(arguments$res)
file_refbiomart <- arguments$file_refbiomart

message("file_sc_obj_WT: ", file_sc_obj_WT)
message("res: ", res)
message("file_refbiomart: ", file_refbiomart)

# --- Custom settings
subset_key <- "sub_trajectory"
subset_names <- c("Limb mesenchyme trajectory")
assay <- "RNA"
parameterset <- "mine"


# --- Read files
# get geneids vs name relations
refbiomart <- read.table(file_refbiomart, sep="\t", header=TRUE)
colnames(refbiomart) <- c("id","name")

# --- Read seurat object and subset to only selected cells
sc_obj_ref <- readRDS(file_sc_obj_WT)
cells_to_keep <- sc_obj_ref@meta.data %>%
  filter(all_of(if_any(subset_key)) %in% subset_names) %>%
  rownames()
sc_obj_ref <- sc_obj_ref %>% subset(cells=cells_to_keep)

if(parameterset=="CX"){
  # specific parameter sets (my parameters)
  k.param=20
  nfeatures <- 2000
  if(ncol(sc_obj_ref) > 50000){
    npcs = 50
    } else if (ncol(sc_obj_ref) > 1000) {
    npcs = 30
    } else {
    npcs = 10
    }
  npcs_clustering <- npcs
  algorithm <- 1 # 4=Leiden, 1=original Louvain default
  vars.to.regress <- c("MT_percent", "RIBO_percent", "log_umi")
} else if(parameterset=="mine"){
  # specific parameter sets (my parameters)
  k.param=20
  nfeatures <- 2000
  npcs <- 50
  npcs_clustering <- npcs
  algorithm <- 1 # 4=Leiden, 1=original Louvain default
  vars.to.regress <- NULL
}

# --- Do Seurat analysis and subclustering
DefaultAssay(sc_obj_ref) <- assay
set.seed(111)
sc_obj_ref <- sc_obj_ref %>%
  NormalizeData(normalization.method="LogNormalize") %>%
  FindVariableFeatures(selection.method="vst", nfeatures=nfeatures) %>%
  ScaleData(vars.to.regress=vars.to.regress) %>%
  RunPCA(npcs=npcs) %>%
  RunUMAP(dims=1:npcs, return.model=TRUE) %>%
  FindNeighbors(dims=1:npcs_clustering, k.param=k.param) %>%
  FindClusters(resolution=res, algorithm=algorithm)

# --- Calculate de genes for the resolutions of interest.
# Find DE Genes, save list and make heatmap
markers <- FindAllMarkers(sc_obj_ref, min.pct=0.2, only.pos=TRUE)
markers <- merge(x=markers, y=refbiomart[,c("name","id")], by.x="row.names", by.y="id", sort=FALSE) %>%
  tibble::column_to_rownames(var="Row.names")

filename <- paste0("markers_WT_res",res,".tsv")
write.table(markers, file=filename, row.name=FALSE, sep="\t")

# --- Annotate sub_sub_clusters
# With cell identities based on de Genes
# Undifferentiated Mesenchyme: Sox11, Sox4, Col1a1, Col1a2
# Condensating Mesenchyme Sox5, Sox6, Sox9, Col2a1
# Perichondrium Prrx1, Wnt5

# Verify that the markers are where they are supposed to be
markers %>% filter(name %in% c("Sox5", "Sox6", "Sox9"))
markers %>% filter(name %in% c("Prrx1", "Wnt5a"))
markers %>% filter(name %in% c("Sox11", "Sox4", "Col1a1", "Col1a2"))

sc_obj_ref <- RenameIdents(sc_obj_ref, "0"="Undifferentiated Mesenchyme", "1"="Perichondrium", "2"="Condensing Mesenchyme")
sc_obj_ref$sub_sub_trajectory <- as.character(Idents(sc_obj_ref))

# --- Make figures and save
plot <- ElbowPlot(sc_obj_ref, ndims=npcs)
ggsave(plot, filename="elbowplot_WT.pdf")

# Make Umap
plot <- DimPlot(sc_obj_ref, reduction="umap", label=TRUE)
ggsave(plot, filename=paste0("umap_after_clustering_WT_res",res,".pdf"), width=7, height=5)

# Do and Save heatmap
markers_to_plot <- markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)

plot <- DoHeatmap(sc_obj_ref, features=markers_to_plot$gene) +
  scale_y_discrete(breaks=markers_to_plot$gene, label=markers_to_plot$name)
filename <- paste0("heatmap_WT_res",res, ".png")
ggsave(plot=plot, filename=filename, width=10,height=7)

# Save the RDS file after sub-clustering and annotation
saveRDS(sc_obj_ref, file="allWT_LimbMes_subsubtrajectory.rds")
