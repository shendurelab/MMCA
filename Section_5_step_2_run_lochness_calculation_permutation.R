###permutation
library(Seurat)
library(Matrix)
library(irlba)
library(ggplot2)
library(dplyr)
library(plyr)
library(monocle3)
library(argparse)
library(FNN)
set.seed(2020)

safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}

wd = './'
sd = './'

# Load data
print('Loading data...')
pca_object_list = readRDS('../data/sub_trajectory_aligned_pca_coor.rds')
umap_object = readRDS('../data/sub_trajectory_umap_coor.rds')
pd = readRDS('../data/pd.rds')
mtx = readRDS('../data/count.rds')
pd$Background[pd$Background=="C57BL/6"]<-"C57BL6"
pd$Background[pd$Background=="BALB/C"]<-"BALBC"
pd$Background_Mutant = paste0(pd$Background,':',pd$Mutant)

full_seurat_obj = CreateSeuratObject(mtx,meta.data=pd) 
full_seurat_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=as.matrix(umap_object), key='UMAP_', assay='RNA')

parser = argparse::ArgumentParser(description='Script that calculates lochness under permutation 10 times.')
parser$add_argument('id', help='the id of the background mutant pair (1-26).')
args = parser$parse_args()

ii = as.integer(args$id)
nn = 10

i = unique(full_seurat_obj$Background_Mutant)[ii]
mt = sub(".*\\:", "", i)
bg = sub("\\:.*", "", i)
if (mt != 'WT') {
  for (tj in names(pca_object_list)) {
    print(paste0(bg,'_',mt,'_',tj))
    seurat_obj = full_seurat_obj[,full_seurat_obj$main_trajectory==tj]
    pca_object = pca_object_list[[tj]]
    pca_dims = dim(pca_object)[2]
    seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(pca_object), key='PC_', assay='RNA')
    seurat_object = seurat_obj[,seurat_obj$Background_Mutant==i | seurat_obj$Mutant == "WT"]
    kadj = round(0.5 * sqrt(ncol(seurat_object)))
    ####PCA space
    seurat_object = seurat_object %>%Seurat::L2Dim(reduction='pca')
    for (shuffle_ind in c(1:nn)) {
      set.seed(shuffle_ind+2022)
      ##shuffle
      seurat_object$sMutant = as.vector(seurat_object$Mutant[sample(length(seurat_object$Mutant))])
      df = data.frame()
      for (j in unique(seurat_object$RT_group)) {
        mt_cells = colnames(seurat_object)[seurat_object$RT_group==j]
        wt_cells = colnames(seurat_object)[seurat_object$RT_group!=j]
        mt.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[mt_cells, ]))
        wt.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[wt_cells, ]))
        # Run KNN and assign labels based on majority of NNs
        knn = get.knnx(t(wt.mtx), t(mt.mtx), k=kadj, algorithm="kd_tree")
        knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(wt.mtx)),to=as.character(seurat_object[,wt_cells]$sMutant),warn_missing=FALSE)
        # MT scores
        mt_counts = as.data.frame(table(seurat_object[,wt_cells]$sMutant))
        mutant_mask = ifelse(knn$label=="WT", 0, 1)
        mutant_neighbors = rowSums(mutant_mask)
        pca_mutant_score = as.data.frame((mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / sum(mt_counts$Freq)) - 1)
        rownames(pca_mutant_score) = mt_cells
        colnames(pca_mutant_score) = c("raw_mt_score")
        df = rbind(df,pca_mutant_score)
      }
      #saving
      df_umap = as.data.frame(seurat_object[['umap']]@cell.embeddings)
      df_umap$raw_mt_score = df[rownames(df_umap),]
      df_umap$main_trajectory = seurat_object$main_trajectory
      saveRDS(df_umap,paste0(sd,bg,'_',mt,'_',tj,'_lochness_excludeself_pool_shuffle',shuffle_ind,'.rds'))
    }
  }
}

for (ii in c(1:10)) {
  df_list = list()
  for (i in unique(full_seurat_obj$Background_Mutant)) {
    mt = sub(".*\\:", "", i)
    bg = sub("\\:.*", "", i)
    if (mt != 'WT') {
      for (tj in names(pca_object_list)) {
        df_umap=readRDS(paste0(sd,bg,'_',mt,'_',tj,'_lochness_excludeself_pool_shuffle',ii,'.rds'))
        df_umap$mt = mt
        df_umap$bg = bg
        df_umap$cell_id = rownames(df_umap)
        df_list[[paste0(tj,i)]] = df_umap
      }}}
  dff = do.call(rbind,df_list)
  print(head(dff$raw_mt_score))
  saveRDS(dff,paste0(sd,'lochness_excludeself_pool_shuffle',ii,'.rds'))
}
