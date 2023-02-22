#on each main trajectory
#7. run PCA on WT cells only, and save the PC space;
#8. for each of the mutant type, projecting their cells into the PC space
#9. run align_cds on the embedding to get the aligned PC space.
#10. do umap on the WT of each trajectory and then projecting MU cells on to that space

#########################################################################
### 7 - running PCA on WT cells only and save that space

library(monocle3)
library(Matrix)
library(dplyr)

main_trajectory_list = c("Neural tube and notochord trajectory",
                         "Endothelial trajectory",
                         "Haematopoiesis trajectory",
                         "Myotube trajectory",
                         "Mesenchymal trajectory",
                         "Neural crest (PNS neuron) trajectory",
                         "Myoblast trajectory",
                         "Epithelial trajectory",
                         "Hepatocyte trajectory",
                         "Melanocyte trajectory",
                         "Neural crest (PNS glia) trajectory",
                         "Olfactory sensory neuron trajectory",
                         "Cardiomyocyte trajectory")

count = readRDS(paste0(work_path, "/data/backup/count.rds"))
pd = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

args = commandArgs(trailingOnly=TRUE)
kk = main_trajectory_list[as.numeric(args[1])]
name = gsub('[(|)]', '', kk)
name = gsub(' ', '_', name)
print(kk)

pd_sub = pd[pd$Mutant == "WT" & pd$main_trajectory == kk,]
count_sub = count[,rownames(pd_sub)]
gene_sub = gene[rownames(count_sub),]

cds <- new_cell_data_set(count_sub,
                         cell_metadata = pd_sub,
                         gene_metadata = gene_sub)

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
fm_rowsums = Matrix::rowSums(FM)
FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

if(ncol(cds) > 50000){
    num_dim = 50
} else if (ncol(cds) > 1000) {
    num_dim = 30
} else {
    num_dim = 10
}

scaling = TRUE
set.seed(2016)
irlba_res <- my_sparse_prcomp_irlba(Matrix::t(FM), 
                                    n = min(num_dim, min(dim(FM)) - 1), 
                                    center = scaling, 
                                    scale. = scaling)
preproc_res <- irlba_res$x
row.names(preproc_res) <- colnames(cds)

gene_use = gene[rownames(FM),]
saveRDS(gene_use, paste0(work_path, "/data/sub_trajectory/", name, "/gene_use_pca.rds"))
saveRDS(preproc_res, paste0(work_path, "/data/sub_trajectory/", name, "/WT_pca_coor.rds"))
saveRDS(irlba_res, paste0(work_path, "/data/sub_trajectory/", name, "/pca_base.rds"))

################################
### 8. for each mutant type, projecting their cells onto the pca_base

mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

for(jj in 1:length(mutant_list)){
    mutant_i = mutant_list[jj]
    print(mutant_i)
    
    pd_sub = pd[pd$Mutant_id == mutant_i & pd$main_trajectory == kk,]
    count_sub = count[,rownames(pd_sub)]
    gene_sub = gene[rownames(count_sub),]
 
    cds <- new_cell_data_set(count_sub,
                             cell_metadata = pd_sub,
                             gene_metadata = gene_sub)
    
    set.seed(2016)
    FM = monocle3:::normalize_expr_data(cds, 
                                        norm_method = "log", 
                                        pseudo_count = 1)
    FM = FM[rownames(gene_use),]
    
    block_size = 50000
    if(ncol(FM) > block_size){
        FM_1 = FM[,1:block_size]
        FM_2 = FM[,(block_size+1):ncol(FM)]
        
        preproc_res_query_1 <- scale(t(FM_1), irlba_res$center, irlba_res$scale) %*% irlba_res$rotation
        preproc_res_query_2 <- scale(t(FM_2), irlba_res$center, irlba_res$scale) %*% irlba_res$rotation
        
        preproc_res_query = rbind(preproc_res_query_1, preproc_res_query_2)
    } else {
        preproc_res_query <- scale(t(FM), irlba_res$center, irlba_res$scale) %*% irlba_res$rotation
    }
    
    row.names(preproc_res_query) <- colnames(cds)
    
    saveRDS(preproc_res_query, paste0(work_path, "/data/sub_trajectory/", name,'/', mutant_i, "_pca_coor.rds"))    
    
}



#########################################################################
### 9  running align_cds on the pca_coor to regress out Mito% and Ribo% and log_umi
### for mesenchymal trajectory, add g2m and s as well

library(monocle3)
library(Matrix)
library(dplyr)

main_trajectory_list = c("Neural tube and notochord trajectory",
                         "Endothelial trajectory",
                         "Haematopoiesis trajectory",
                         "Myotube trajectory",
                         "Mesenchymal trajectory",
                         "Neural crest (PNS neuron) trajectory",
                         "Myoblast trajectory",
                         "Epithelial trajectory",
                         "Hepatocyte trajectory",
                         "Melanocyte trajectory",
                         "Neural crest (PNS glia) trajectory",
                         "Olfactory sensory neuron trajectory",
                         "Cardiomyocyte trajectory")

args = commandArgs(trailingOnly=TRUE)
kk = main_trajectory_list[as.numeric(args[1])]
name = gsub('[(|)]', '', kk)
name = gsub(' ', '_', name)
print(kk)

pd = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
pd = pd[pd$main_trajectory == kk,]
mutant_list = as.vector(unique(pd$Mutant_id))

pca_coor = NULL
for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(mutant_i)
    
    tmp = readRDS(paste0(work_path, "/data/sub_trajectory/", name,'/', mutant_i, "_pca_coor.rds"))
    pca_coor = rbind(pca_coor, tmp)
}
if (sum(!rownames(pca_coor) %in% rownames(pd)) != 0){
    FUN()
}
pca_coor = pca_coor[rownames(pd),]
saveRDS(pca_coor, paste0(work_path, "/data/sub_trajectory/", name, "/combined_pca_coor.rds"))

set.seed(2016)
if(kk == "Mesenchymal trajectory"){
    print(kk)
    residual_model_formula_str = "~RIBO_percent + MT_percent + log_umi + g2m_score + s_score"
} else {
    residual_model_formula_str = "~RIBO_percent + MT_percent + log_umi"
}

X.model_mat <- Matrix::sparse.model.matrix(stats::as.formula(residual_model_formula_str), 
                                           data = pd, drop.unused.levels = TRUE)
fit <- limma::lmFit(Matrix::t(pca_coor), X.model_mat)

beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0
aligned_coor <- Matrix::t(as.matrix(Matrix::t(pca_coor)) - 
                              beta %*% Matrix::t(X.model_mat[, -1]))

saveRDS(aligned_coor, paste0(work_path, "/data/sub_trajectory/", name, "/combined_aligned_coor.rds"))



#####################################
### step 10, doing UMAP on WT and save the space using uwot
### then projecting mutant cells into that space

library(monocle3)
library(Matrix)
library(dplyr)

main_trajectory_list = c("Neural tube and notochord trajectory",
                         "Endothelial trajectory",
                         "Haematopoiesis trajectory",
                         "Myotube trajectory",
                         "Mesenchymal trajectory",
                         "Neural crest (PNS neuron) trajectory",
                         "Myoblast trajectory",
                         "Epithelial trajectory",
                         "Hepatocyte trajectory",
                         "Melanocyte trajectory",
                         "Neural crest (PNS glia) trajectory",
                         "Olfactory sensory neuron trajectory",
                         "Cardiomyocyte trajectory")

args = commandArgs(trailingOnly=TRUE)
kk = main_trajectory_list[as.numeric(args[1])]
name = gsub('[(|)]', '', kk)
name = gsub(' ', '_', name)
print(kk)

pd = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
pd = pd[pd$main_trajectory == kk,]
aligned_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_aligned_coor.rds"))

pd_sub = pd[pd$Mutant == "WT",]
aligned_coor_sub = aligned_coor[rownames(pd_sub),]

set.seed(2016)
emb_train_umap = uwot::umap(as.matrix(aligned_coor_sub), 
                            n_components = 3,
                            n_neighbors = 15,
                            min_dist = 0.1,
                            metric = "cosine",
                            fast_sgd = FALSE,
                            nn_method = "annoy",
                            ret_model = TRUE,
                            n_threads = 1,
                            verbose = TRUE)

umap_coor = emb_train_umap$embedding
rownames(umap_coor) = rownames(pd_sub)
colnames(umap_coor) = paste0("UMAP_", 1:3)

saveRDS(umap_coor, paste0(work_path, "/data/sub_trajectory/", name, "/WT_umap_coor.rds"))
uwot::save_uwot(emb_train_umap, paste0(work_path, "/data/sub_trajectory/", name, "/umap_base"))

###################
### projecting cells from each mutant type and get their umap coordinates

mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(paste0(i, "/", mutant_i))
    
    pd_sub = pd[pd$Mutant_id == mutant_i,]
    aligned_coor_sub = aligned_coor[rownames(pd_sub),]
    
    set.seed(2016)
    emb_test_coor = uwot::umap_transform(as.matrix(aligned_coor_sub),
                                         emb_train_umap)
    
    rownames(emb_test_coor) = rownames(pd_sub)
    colnames(emb_test_coor) = paste0("UMAP_", 1:3)
    
    saveRDS(emb_test_coor, paste0(work_path, "/data/sub_trajectory/", name, "/", mutant_i, "_umap_coor.rds"))
    
}


### combining all the mutants as well as WT
mutant_list = as.vector(unique(pd$Mutant_id))
umap_coor = NULL
for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(mutant_i)
    
    tmp = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/", mutant_i, "_umap_coor.rds"))
    umap_coor = rbind(umap_coor, tmp)
}
print(sum(!rownames(umap_coor) %in% rownames(pd)))
umap_coor = umap_coor[rownames(pd),]
saveRDS(umap_coor, paste0(work_path, "/data/sub_trajectory/", name, "/combined_umap_coor.rds"))


