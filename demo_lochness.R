library(ggplot2)
library(dplyr)
library(tidyverse)
library(Seurat)
library(FNN)
library(plyr)
library(scales)
library(reshape2)

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

###########################
#Manuscriot implementation#
###########################
# The current implementation of lochNESS is a simple point estimate, 
# considering all wildtype cells as one entity. 

dff = data.frame()
for (j in unique(seurat_object$RT_group)) {
  #split cells from one embryos and other embryos
  mt_cells = colnames(seurat_object)[seurat_object$RT_group==j] 
  ot_cells = colnames(seurat_object)[seurat_object$RT_group!=j]
  mt.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[mt_cells, ]))
  ot.mtx = as.matrix(t(as.data.frame(seurat_object@reductions$pca.l2@cell.embeddings)[ot_cells, ]))
  
  # Run KNN and assign labels based on majority of NNs
  knn = get.knnx(t(ot.mtx), t(mt.mtx), k=kadj, algorithm="kd_tree")
  knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(ot.mtx)),to=as.character(seurat_object[,ot_cells]$Mutant),warn_missing=FALSE)
  
  # calculate lochNESS
  mt_counts = as.data.frame(table(seurat_object[,ot_cells]$Mutant))
  mutant_mask = ifelse(knn$label=="WT", 0, 1)
  mutant_neighbors = rowSums(mutant_mask)
  pca_mutant_score = as.data.frame((mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / sum(mt_counts$Freq)) - 1)
  rownames(pca_mutant_score) = mt_cells
  colnames(pca_mutant_score) = c("lochNESS")
  dff = rbind(dff,pca_mutant_score)
}

#saving..
df_umap = as.data.frame(seurat_object[['umap']]@cell.embeddings)
df_umap$lochNESS = dff[rownames(df_umap),]
df_umap$sub_trajectory = seurat_object$sub_trajectory
#saveRDS(df_umap,paste0(wd,'demo_lochness_scores.rds'))

############################
#Alternative implementation#
############################
# Alternatively, if the dataset contains sufficient cells, one can
# calculate lochNESS_i (one score per wildtype sample), and assess
# the variability of lochNESS when using different wildtype samples.
# Wildtype RT IDs are 42-44, and mutant RT IDs are 94-97.

wildtype_ids = unique(seurat_object[,seurat_object$Mutant=="WT"]$RT_group)
mutant_ids = unique(seurat_object[,seurat_object$Mutant!="WT"]$RT_group)

df_list = list()
for (ot_rt in wildtype_ids) { # calculate a score for one wildtype sample
  df = data.frame()
  seurat_obj = seurat_object[,seurat_object$RT_group %in% c(mutant_ids,ot_rt)] # specifiy which samples go into the analysis
  for (j in mutant_ids) {
    #split mutant and wildtype cells
    mt_cells = colnames(seurat_obj)[seurat_obj$RT_group==j]
    ot_cells = colnames(seurat_obj)[seurat_obj$RT_group!=j]
    mt.mtx = as.matrix(t(as.data.frame(seurat_obj@reductions$pca.l2@cell.embeddings)[mt_cells, ]))
    ot.mtx = as.matrix(t(as.data.frame(seurat_obj@reductions$pca.l2@cell.embeddings)[ot_cells, ]))
    
    # Run KNN and assign labels based on majority of NNs
    knn = get.knnx(t(ot.mtx), t(mt.mtx), k=kadj, algorithm="kd_tree")
    knn$label <- mapvalues(knn$nn.index,from=seq(1, ncol(ot.mtx)),to=as.character(seurat_obj[,ot_cells]$Mutant),warn_missing=FALSE)
    
    # calculate lochNESS
    mt_counts = as.data.frame(table(seurat_obj[,ot_cells]$Mutant))
    mutant_mask = ifelse(knn$label=="WT", 0, 1)
    mutant_neighbors = rowSums(mutant_mask)
    pca_mutant_score = as.data.frame((mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / sum(mt_counts$Freq)) - 1)
    rownames(pca_mutant_score) = mt_cells
    colnames(pca_mutant_score) = c(paste0("lochNESS_",ot_rt))
    df = rbind(df,pca_mutant_score)
  }
  df_list[[ot_rt]] = df
}
cells_to_plot = colnames(seurat_object)[seurat_object$RT_group %in% mutant_ids]
df_umap = as.data.frame(seurat_object[['umap']]@cell.embeddings[cells_to_plot,])

#plotting and comparing the two versions 
df_umap$lochNESS = dff[rownames(df_umap),]
df_umap$sub_trajectory = seurat_object$sub_trajectory[cells_to_plot]
plot1 <- ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2,color=lochNESS)) +
  geom_point(size=0.05) +
  scale_color_gradientn(
    colours = c(muted("red"), "white", muted("blue")),
    values = rescale(c(min(df_umap$lochNESS),median(df_umap$lochNESS),max(df_umap$lochNESS))),
    guide = "colorbar", limits=c(min(df_umap$lochNESS),max(df_umap$lochNESS))
  ) + theme_classic() 
#ggsave(plot1,filename = paste0(wd,'demo_lochness_umap2d.pdf'),width=6,height=6)

plot2 <- ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2,color=sub_trajectory)) + 
  geom_point(size=0.2) + theme_classic() 
#ggsave(plot2,filename = paste0(wd,'demo_subtrajectory_umap2d.pdf'),width=8,height=6)

for (ot_rt in wildtype_ids) { 
  temp = df_list[[ot_rt]][cells_to_plot,paste0("lochNESS_",ot_rt)]
  df_umap[[paste0("lochNESS_",ot_rt)]] = temp
  plot <- ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2,color=temp)) +
    geom_point(size=0.05) +
    scale_color_gradientn(
      colours = c(muted("red"), "white", muted("blue")),
      values = rescale(c(min(temp),median(temp),max(temp))),
      guide = "colorbar", limits=c(min(temp),max(temp))
    ) + theme_classic() 
  #ggsave(plot,filename = paste0(wd,'demo_lochness_',ot_rt,'_alt_umap2d.pdf'),width=6,height=6)
}

#plotting correlation matrix
cormat <- round(cor(df_umap[,grep("lochNESS", names(df_umap), value = TRUE)]),2)
melted_cormat <- melt(cormat)
plot <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme_classic()
#ggsave(plot,filename = paste0(wd,'demo_lochness_corr.pdf'),width=7,height=6)



