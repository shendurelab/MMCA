library(Seurat)
library(Matrix)
library(irlba)
library(ggplot2)
library(dplyr)
library(plyr)
library(monocle3)
library(argparse)
set.seed(2020)

safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}

wd = './'
sd = './'

main_tjs = c("Mesenchymal trajectory","Neural tube and notochord trajectory","Epithelial trajectory","Neural crest (PNS neuron) trajectory", 
             "Hepatocyte trajectory","Haematopoiesis trajectory","Endothelial trajectory","Melanocyte trajectory",                
             "Neural crest (PNS glia) trajectory","Myoblast trajectory","Myotube trajectory","Olfactory sensory neuron trajectory",
             "Cardiomyocyte trajectory")

# Load data
print('Loading data...')
pca_object_list = readRDS('../data/sub_trajectory_aligned_pca_coor.rds')
umap_object = readRDS('../data/sub_trajectory_umap_coor.rds')
pd = readRDS('../data/pd.rds')
mtx = readRDS('../data/count.rds')
pd$Background[pd$Background=="C57BL/6"]<-"C57BL6"
pd$Background[pd$Background=="BALB/C"]<-"BALBC"
pd$Background_Mutant = paste0(pd$Background,':',pd$Mutant)

seurat_obj_full = CreateSeuratObject(mtx,meta.data=pd) 
seurat_obj_full[['umap']] = Seurat::CreateDimReducObject(embeddings=as.matrix(umap_object), key='UMAP_', assay='RNA')

#####
#####mutant similarities
#####
for (i in main_tjs) {
  for (j in c("G4","C57BL6","FVB")) {
    seurat_object = seurat_obj[,seurat_obj$main_trajectory==i & seurat_obj$Background==j]
    seurat_object
    
    kadj = round(0.5 * sqrt(ncol(seurat_object)))
    
    ####PCA space
    seurat_object = seurat_object %>%
      Seurat::L2Dim(reduction='pca') %>%
      Seurat::FindNeighbors(reduction='pca.l2', dims=1:50, k=kadj) # nn.eps = 0.25   
    
    target_counts = as.data.frame(table(seurat_object$Mutant))
    targets = unique(seurat_object$Mutant)
    target_score <- list(length=c(length(seurat_object$Mutant)))
    for (j in 1:length(targets)) {
      target_mask = ifelse(seurat_object$Mutant==targets[j], 1, 0)
      target_neighbors = Matrix::rowSums(safe_column_multiply(seurat_object@graphs$RNA_nn, target_mask))
      target_score[[j]] = target_neighbors / kadj / target_counts$Freq[target_counts$Var1==targets[j]] * length(target_neighbors)
    }
    target_score_mtx = do.call(cbind,target_score)
    colnames(target_score_mtx) = targets
    df = as.data.frame(target_score_mtx)
    df$Mutant = seurat_object$Mutant
    
    df2 = df %>% group_by(Mutant) %>% summarise_all(mean)
    rnames = df2$Mutant
    df2 = as.matrix(df2[2:ncol(df2)])
    rownames(df2) = rnames
    heatmap(df2)
    saveRDS(df2,paste0(sd,i,'_',j,'_similarity_mtx.rds'))
  }
}
