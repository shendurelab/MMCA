
library(Seurat)
library(plotly)
library(htmlwidgets)
library(plyr)
library(FNN)
library(argparse)
library(dplyr)

parser = argparse::ArgumentParser(description='Script that takes finds MMCA kNNs in MOCA.')
parser$add_argument('run_id', help='mutant type to run analysis on (1-10).')
args = parser$parse_args()
j = as.integer(args$run_id)

pd = readRDS('../data/final/pd.rds')
mtx = readRDS('../data/final/count.rds')
mt_meta = readRDS('../data/lochness_meta.rds')
mt_meta$mt1 = gsub(" ", "_", mt_meta$mt1)
mt_meta$mt2 = gsub(" ", "_", mt_meta$mt2)
pd$Background_genotype = gsub("/", "", pd$Background_genotype)

moca_exp = readRDS('https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_count_cleaned.RDS')
moca_gene = read.csv('https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/gene_annotate.csv',header=T)
moca_meta = read.csv('https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/cell_annotate.csv',header=T)
moca_meta_new = read.csv('https://shendure-web.gs.washington.edu/content/members/cao1025/public/mouse_embryo_atlas/cell_annotate_20200119.csv',header=T)
moca_meta = moca_meta[moca_meta$sample %in% colnames(moca_exp),]

mt_list = c("SScn11A GoF","ZRS limb enhancer KO","Sox9 Regulatroy KO","Tbx3 TAD boundary KO","Gorab KO",
            "Atp6v0a2 KO","Gli2 KO","Atp6v0a2RQ KI","TAD boundary KI","Tbx5 TAD boundary KO")
moca_ct_list = list(c("Endothelial trajectory"),c("Epithelial trajectory"),c("Haematopoiesis trajectory"),c("Mesenchymal trajectory"),c("Neural crest 1"),c("Neural crest 2"),c("Neural crest 3"),c("Neural tube and notochord trajectory"))
mmca_ct_list = list(c("Endothelial trajectory"),c("Epithelial trajectory"),c("Haematopoiesis trajectory"),c("Mesenchymal trajectory","Myoblast trajectory","Myotube trajectory","Cardiomyocyte trajectory"),c("Neural crest (PNS neuron) trajectory","Olfactory sensory neuron trajectory"),c("Neural crest (PNS glia) trajectory"),c("Melanocyte trajectory"),c("Neural tube and notochord trajectory"))
ct_names = c("Endothelial trajectory","Epithelial trajectory","Haematopoiesis trajectory","Mesenchymal trajectory","Neural crest (PNS neuron) trajectory","Neural crest (PNS glia) trajectory","Melanocyte trajectory","Neural tube trajectory")
tps = c("9.5","10.5","11.5","12.5","13.5")

for (i in 1:length(mt_list)) {
  mt = mt_list[i]
  moca_ct = moca_ct_list[[j]]
  mmca_ct = mmca_ct_list[[j]]
  ct_name = ct_names[j]
  
  mmca_pd = pd[pd$main_trajectory %in% mmca_ct & pd$Mutant %in% c(mt,'WT'),]
  mmca_meta = mmca_pd[,c("Background_genotype","sub_trajectory")]
  colnames(mmca_meta) = c("sample","anno")
  mmca_mtx = mtx[,rownames(mmca_meta)]
  rownames(mmca_mtx) = gsub("\\..*","",rownames(mmca_mtx))
  mmca_seurat = CreateSeuratObject(mmca_mtx,meta.data=mmca_meta)
  
  moca_meta_sub = moca_meta[moca_meta$Main_trajectory %in% moca_ct,]
  rownames(moca_meta_sub) = moca_meta_sub$sample
  moca_meta_sub = moca_meta_sub[,c("development_stage","Main_trajectory")]
  colnames(moca_meta_sub) = c("sample","anno")
  moca_mtx = moca_exp[,rownames(moca_meta_sub)]
  rownames(moca_mtx) = gsub("\\..*","",rownames(moca_mtx))
  moca_seurat = CreateSeuratObject(moca_mtx,meta.data=moca_meta_sub)
  
  gene.intersect = intersect(rownames(moca_seurat),rownames(mmca_seurat))
  merged_seurat = merge(moca_seurat[gene.intersect,],mmca_seurat[gene.intersect,],add.cell.ids = c("MOCA","MMCA"))
  
  merged_seurat <- NormalizeData(object = merged_seurat)
  merged_seurat <- FindVariableFeatures(object = merged_seurat)
  merged_seurat <- ScaleData(object = merged_seurat)
  merged_seurat <- RunPCA(object = merged_seurat)
  merged_seurat <- FindNeighbors(object = merged_seurat)
  merged_seurat <- FindClusters(object = merged_seurat)
  merged_seurat <- RunUMAP(object = merged_seurat, dims=1:30,n.components = 3L)
  
  df = cbind(Embeddings(object = merged_seurat, reduction = "umap"),merged_seurat@meta.data)
  
  saveRDS(merged_seurat,paste0(sd,'/integrate_origMOCA_',ct_name,'_',mt,'_WT_seurat.rds'))
  
  ######knn
  moca_pca <- Embeddings(object = merged_seurat[,merged_seurat$sample %in% tps], reduction = "pca")
  mmca_pca <- Embeddings(object = merged_seurat[,!merged_seurat$sample %in% tps], reduction = "pca")
  moca_meta_knn = merged_seurat[,merged_seurat$sample %in% tps]@meta.data
  moca_meta_knn$time = as.numeric(gsub("E","",moca_meta_knn$sample))
  
  k_neigh = 10
  neighbors <- get.knnx(moca_pca, mmca_pca, k = k_neigh)$nn.index
  tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
  for(kk in 1:k_neigh){
    tmp1[,kk] <- moca_meta_knn$time[neighbors[,kk]]
  }
  rownames(tmp1) = rownames(mmca_pca)
  neighbors_df = data.frame(tmp1)
  names(neighbors_df) = paste0("x", 1:k_neigh)
  
  mmca_meta = merged_seurat[,!merged_seurat$sample %in% tps]@meta.data
  mmca_meta$time_scores = rowMeans(neighbors_df)
  
  mu <- ddply(mmca_meta, "sample", summarise, grp.mean=mean(time_scores))
  p<-ggplot(mmca_meta, aes(x=time_scores, fill=sample, color=sample)) +
    geom_histogram(position="identity", alpha=0.5) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=sample),linetype="dashed") + theme_classic()
  p+ggsave(paste0(sd,'/time_hist_origMOCA_',ct_name,'_',mt,'_WT.png'))
  
  saveRDS(mmca_meta,paste0(sd,'/integrate_origMOCA_',ct_name,'_',mt,'_WT_meta.rds'))
}
