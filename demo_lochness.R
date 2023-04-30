library(ggplot2)
library(dplyr)
library(tidyverse)
library(Seurat)
library(FNN)
library(plyr)
library(scales)

safe_column_scale = function(bmat, scale_vector) {
  bmat@x <- bmat@x / rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}
safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}
wd = 'demo_data/'
seurat_object = readRDS(paste0(wd, "demo_Seurat.rds"))

# The processed seurat object should have precomtuted PCA and UMAP coordinates
kadj = round(0.5 * sqrt(ncol(seurat_object)))
seurat_object = seurat_object %>% Seurat::L2Dim(reduction='pca')

df = data.frame()
for (j in unique(seurat_object$RT_group)) {
  #split mutant and wildtype cells
  mt_cells = colnames(seurat_object)[seurat_object$RT_group==j]
  wt_cells = colnames(seurat_object)[seurat_object$RT_group!=j]
  mt.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[mt_cells, ]))
  wt.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[wt_cells, ]))
  
  # Run KNN and assign labels based on majority of NNs
  knn = get.knnx(t(wt.mtx), t(mt.mtx), k=kadj, algorithm="kd_tree")
  knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(wt.mtx)),to=as.character(seurat_object[,wt_cells]$Mutant),warn_missing=FALSE)
  
  # calculate lochNESS
  mt_counts = as.data.frame(table(seurat_object[,wt_cells]$Mutant))
  mutant_mask = ifelse(knn$label=="WT", 0, 1)
  mutant_neighbors = rowSums(mutant_mask)
  pca_mutant_score = as.data.frame((mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / sum(mt_counts$Freq)) - 1)
  rownames(pca_mutant_score) = mt_cells
  colnames(pca_mutant_score) = c("lochNESS")
  df = rbind(df,pca_mutant_score)
}

#saving and plotting
df_umap = as.data.frame(seurat_object[['umap']]@cell.embeddings)
df_umap$lochNESS = df[rownames(df_umap),]
df_umap$sub_trajectory = seurat_object$sub_trajectory
#saveRDS(df_umap,paste0(wd,'demo_lochness_scores.rds'))

plot <- ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2,color=lochNESS)) +
  geom_point(size=0.05) +
  scale_color_gradientn(
    colours = c(muted("red"), "white", muted("blue")),
    values = rescale(c(min(df_umap$lochNESS),0,max(df_umap$lochNESS))),
    guide = "colorbar", limits=c(min(df_umap$lochNESS),max(df_umap$lochNESS))
  ) + theme_classic() 
#ggsave(plot,filename = paste0(wd,'demo_lochness_umap2.pdf'))

plot2 <- ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2,color=sub_trajectory)) + theme_classic() 
