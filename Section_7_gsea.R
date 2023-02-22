# Call from conda-environment MegaBundle
" Perform a GSEA analysis on the submitted Seurat object, and output the Seurat Object with enrichment scores as well as a  dittoheatmap of enrichment scores (clustered and unclustered). Using ESCAPE package. Also do some statistics

Usage: runGSEA.R --file_sc_obj=<file> --file_refbiomart=<file> --geneSet=<value> --groupby=<value> --assay=<value>

Options:
  -h --help			Show this screen.
  --file_sc_obj=<file>		Desired Seurat Object
  --file_refbiomart=<file>  The refbiomart file containing geneid and genenames
  --geneSet=<value>  Library,NamePattern Check www.gsea-msigdb.org/gsea/msigdb/collections.jsp
  --groupby=<value>   Which metadata column of Seurat object to group for statistical testing and violin plotting
  --assay=<value>   Assay in Seurat object to use <RNA>, <spliced>, or <unspliced>
" -> doc

# --- load libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressPackageStartupMessages(library(escape))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(ggrepel))

# suppressMessages(library(SeuratWrappers)) # convert seurat to monocle3
# suppressMessages(library(monocle3))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the output of the sc_multi_sample pipline
file_sc_obj <- arguments$file_sc_obj
file_refbiomart <- arguments$file_refbiomart
geneSet <- strsplit(arguments$geneSet, split=",") %>% unlist
geneSetLibrary <- geneSet[1]
geneSetNamePattern <- geneSet[2]
if(is.na(geneSetNamePattern) || geneSetNamePattern=='NULL' || geneSetNamePattern =='null'){
  geneSetNamePattern <- NULL
}
groupby <- arguments$groupby
assay <- arguments$assay

message("file_sc_obj: ", file_sc_obj)
message("file_refbiomart: ", file_refbiomart)
message("geneSetLibrary: ", geneSetLibrary)
message("geneSetNamePattern: ", geneSetNamePattern)
message("groupby: ", groupby)
message("assay: ", assay)

# parameter set
res <- 0.2
parameterset <- "mine" # make umap and DE genes using this. Almost never used the "CX" parameterset
redOI <- "ref.umap"
ncores <- parallel::detectCores()
message("ncores: ",ncores)

# --- Read files
sc_obj <- readRDS(file_sc_obj)

# --- get geneids vs genenames from the refbiomart file
refbiomart <- read.table(file_refbiomart, sep="\t", header=TRUE)
colnames(refbiomart) <- c("id","name")
rownames(refbiomart) <- refbiomart$id

# -- In this step, replace the ensembl ids with gene symbols. Do this by creating a new Seurat object
counts <- GetAssayData(sc_obj, assay=assay, slot="counts")
common_genes <- intersect(rownames(counts), refbiomart$id)

counts <- counts[common_genes,]
gene_names <- refbiomart[rownames(counts), "name"]
rownames(counts) <- gene_names
meta.data <- sc_obj@meta.data

# Normalize the data
sc_obj_genenames <- CreateSeuratObject(counts=counts, meta.data=meta.data) %>%
  NormalizeData(normalization.method = "LogNormalize")

# -- start the escape analysis
all_gene_sets = msigdbr(species = "Mus musculus") #Get all gene sets in the MSIG database

if(is.null(geneSetNamePattern)){
  selected_gene_sets <- NULL
} else {
  selected_gene_sets <- all_gene_sets %>% filter(gs_cat==geneSetLibrary) %>% filter(grepl(gs_name,pattern=geneSetNamePattern,ignore.case=TRUE)) %>% distinct(gs_name) %>% .$gs_name
}

gene.sets <- getGeneSets(library = geneSetLibrary, gene.sets=selected_gene_sets, species="Mus musculus") # Get the user defined gene sets

# calculate enrichment scores for the cells
ES <- enrichIt(obj <- sc_obj_genenames,
  gene.sets = gene.sets,
  groups = 10000, cores = ncores,
  min.size = NULL)
colnames(ES) <- paste0("ESCAPE_",colnames(ES)) #This is done so that these columns can be selected in the plotting R-script

# Add the enrichment scores as metadata to the original sc_obj
sc_obj <- AddMetaData(sc_obj, ES)

#save the Seuratobject
filename=gsub(file_sc_obj, pattern=".rds", replacement=paste0("_",assay,"_",geneSetLibrary,"_",geneSetNamePattern,".rds"))
saveRDS(sc_obj, file=filename)

# Make pheatmap
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF",
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

filename=gsub(file_sc_obj, pattern=".rds", replacement=paste0("_",assay,"_",geneSetLibrary,"_",geneSetNamePattern,".pdf"))
pdf(filename, width=10,height=10)
dittoHeatmap(sc_obj, genes = NULL, metas = names(ES),
             heatmap.colors = rev(colorblind_vector(50)),
             annot.by = c(groupby, "sub_sub_trajectory"),
             cluster_cols = TRUE,
             fontsize = 7)

dittoHeatmap(sc_obj, genes = NULL, metas = names(ES),
  heatmap.colors = rev(colorblind_vector(50)),
  annot.by = c(groupby),
  cluster_cols = FALSE,
  fontsize = 7)
dev.off()

# --- Do statistics on the enrichment values. First create a new dataframe containing the enrichments and the user-defined grouping variable
enrichments <- cbind(ES, select(sc_obj@meta.data, matches(groupby)))

output <- getSignificance(enrichments, group=groupby, fit="T.test")
output1 <- output %>% arrange(FDR)
# Save the filtered and ordered t-test output to a tsv file
filename=paste0("significance_",gsub(file_sc_obj, pattern=".rds", replacement=paste0("_",assay,"_",geneSetLibrary,"_",geneSetNamePattern,".tsv")))
write.table(output1, file=filename, sep="\t")

# Make violin plots for the top 6 significant gene sets
selected <- output %>% filter(FDR < 0.05) %>% arrange(FDR) %>% rownames() %>% head(n=12)
plots <- VlnPlot(sc_obj, features=selected, group.by=groupby, combine=FALSE)
for (i in seq(length(plots))){
    plots[[i]] <- plots[[i]] + ggtitle(gsub(selected[i], pattern="_", replacement="\n")) +
    theme(text=element_text(family="sans", size=5), legend.position="none",
      plot.title=element_text(hjust=0.5, size=10, family="sans"))
}
filename=paste0("Violin_",gsub(file_sc_obj, pattern=".rds", replacement=paste0("_",assay,"_",geneSetLibrary,"_",geneSetNamePattern,".png")))
ggsave(plot=plot_grid(plotlist=plots, ncol=ceiling(length(selected)/2)), filename=filename,
  width=5*length(selected)/2, height=5*2, limitsize = FALSE)
