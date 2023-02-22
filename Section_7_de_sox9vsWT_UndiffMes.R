# Call from conda-environment MegaBundle
" Find DE genes between the WT and the mutant for cells in the Undifferentiated Mesenchyme

Usage: de_sox9vsWT_UndiffMes.R --file_sc_obj=<file> --file_refbiomart=<file>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		The seurat object containing Limb Mesenchyme cells, with sub_sub_trajectory annotation
  --file_refbiomart=<file>    The refbiomart file containing geneid and genenames
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters
file_sc_obj <- arguments$file_sc_obj
file_refbiomart <- arguments$file_refbiomart

message("file_sc_obj: ", file_sc_obj)
message("file_refbiomart: ", file_refbiomart)

# custom parameters
assay <- "RNA"
parameterset <- "mine" # make umap and DE genes using this. Almost never used the "CX" parameterset
redOI <- "ref.umap"

# --- Read files
sc_obj <- readRDS(file_sc_obj)

# --- get geneids vs genenames from the refbiomart file
refbiomart <- read.table(file_refbiomart, sep="\t", header=TRUE)
colnames(refbiomart) <- c("id","name")

# -- subset only undifferentiated mesenchyme cells corresponding to the two mutants
cells_1 <- sc_obj@meta.data %>%
  filter(Mutant_id=="Sox9_Regulatroy_KO" & sub_sub_trajectory=="Undifferentiated Mesenchyme") %>%
  rownames()
cells_2 <- sc_obj@meta.data %>%
  filter(Mutant_id=="WT" & sub_sub_trajectory=="Undifferentiated Mesenchyme") %>%
  rownames()

# Create another Seurat object containing only Undifferentiated mesechyme from the two samples
sc_obj <- subset(sc_obj,cells=c(cells_1,cells_2))

# -- Setup parameters for Seurat analysis to be able to find markers
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

# --- Do Seurat analysis
DefaultAssay(sc_obj) <- assay
set.seed(111)
sc_obj <- sc_obj %>%
  NormalizeData(normalization.method="LogNormalize") %>%
  FindVariableFeatures(selection.method="vst", nfeatures=nfeatures) %>%
  ScaleData(vars.to.regress=vars.to.regress)

Idents(sc_obj) <- sc_obj$Mutant_id
markers <- FindAllMarkers(sc_obj, min.pct=0.2, only.pos=TRUE)
markers <- merge(x=markers, y=refbiomart[,c("name","id")], by.x="row.names", by.y="id", sort=FALSE) %>%
  tibble::column_to_rownames(var="Row.names")


write.table(x=markers, file="markers_de_sox9vsWT_undiff_mes.tsv", sep="\t")

goi <- markers %>% group_by(cluster) %>% slice_head(n=6)
plot <- FeaturePlot(sc_obj, features=goi$gene, split.by="Mutant_id", reduction=redOI,  keep.scale="all", pt.size=2) & scale_color_viridis()
ggsave(plot, file="featPlot_de_sox9vsWT_undiff_mes.png", width=8, height=48)

plot <- DimPlot(sc_obj, split.by="Mutant_id", reduction=redOI)
ggsave(plot, file="umap_de_sox9vsWT_undiff_mes.png", width=10, height=4)

# Save the RDS file (for eg. to do GSEA analysis)
saveRDS(sc_obj, file="Sox9_Regulatory_KOvsWT_UndiffMes.rds")
