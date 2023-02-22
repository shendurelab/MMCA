# Call from conda-environment MegaBundle
" Subcluster the two branches of Undifferentiated Mesenchyme in the Sox9 Mutant and find DE genes.

Usage: de_sox9_branches_UndiffMes.R --file_sc_obj=<file> --file_refbiomart=<file>

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
res <- 0.2
redOI <- "ref.umap"

# --- Read files
sc_obj <- readRDS(file_sc_obj)

# --- get geneids vs genenames from the refbiomart file
refbiomart <- read.table(file_refbiomart, sep="\t", header=TRUE)
colnames(refbiomart) <- c("id","name")

# -- subset only undifferentiated mesenchyme cells of the Sox9 Mutant - not the WT
cells_1 <- sc_obj@meta.data %>%
  filter(Mutant_id=="Sox9_Regulatroy_KO" & sub_sub_trajectory=="Undifferentiated Mesenchyme") %>%
  rownames()

# Create another Seurat object containing only Undifferentiated mesechyme from Sox9 Mutant
sc_obj <- subset(sc_obj,cells=c(cells_1))

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
  ScaleData(vars.to.regress=vars.to.regress) %>%
  RunPCA(npcs=npcs) %>%
  RunUMAP(dims=1:npcs, seed.use=111) %>%
  FindNeighbors(dims=1:npcs_clustering, k.param=k.param) %>%
  FindClusters(resolution=res, algorithm=algorithm)

sc_obj$branch <- paste0("Branch_",sc_obj@meta.data[, paste0("RNA_snn_res.",res)])
Idents(sc_obj) <- sc_obj$branch

# find marker genes
markers <- FindAllMarkers(sc_obj, min.pct=0.2, only.pos=TRUE)
markers <- merge(x=markers, y=refbiomart[,c("name","id")], by.x="row.names", by.y="id", sort=FALSE) %>%
  tibble::column_to_rownames(var="Row.names")

write.table(x=markers, file="markers_de_sox9_branches_undiff_mes.tsv", sep="\t")

# plot genes of interest
goi <- markers %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% slice_head(n=6)
plots <- FeaturePlot(sc_obj, features=goi$gene, reduction=redOI, combine=FALSE, keep.scale="all", pt.size=2)
for (i in seq(length(plots))){
    plots[[i]] <- plots[[i]] + ggtitle(goi$name[i]) + theme_void() +
    theme(text=element_text(family="sans", size=30),
      legend.position="none", legend.key.height=unit(0.1, units="npc"),
      legend.text=element_text(family="sans", size=30),
      plot.title=element_text(hjust=0.5, size=30, vjust=-2, family="sans"),
      aspect.ratio=1) +
    scale_color_viridis()
}
ggsave(plot=plot_grid(plotlist=plots, ncol=6), file="featPlot_de_sox9_branches_undiff_mes.png", width=5*6, height=5*2, limitsize = FALSE)

goi <- markers %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% slice_head(n=6)
plot <- DotPlot(sc_obj, features=goi$gene) +
  scale_x_discrete(breaks=goi$gene, labels=goi$name) +
  ylab("Undifferentiated mesenchyme branches") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # +
  # scale_color_viridis()
ggsave(plot, filename="DotPlot_de_sox9_branches_undiff_mes.pdf", width=5.2, height=2)

plot <- DimPlot(sc_obj, group.by="branch", reduction=redOI) + theme(aspect.ratio=1)
ggsave(plot, file="umap_de_sox9_branches_undiff_mes.png", width=5+2, height=5)

# Save the RDS file (for eg. to do GSEA analysis)
saveRDS(sc_obj, file="Sox9_Regulatory_KO_UndiffMes_branches.rds")
