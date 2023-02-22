# on each main trajectory
#11. run louvain clustering on WT cells, 
#12. and annotating the cell types;


#####################################
### step 11, doing louain clustering on WT cells based on 3d UMAP embedding
### and then assigning the main trajectories

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

count = readRDS(paste0(work_path, "/data/backup/count.rds"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
pd = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
pd_sub = pd[pd$main_trajectory == kk & pd$Mutant_id == "WT",]
count_sub = count[,rownames(pd_sub)]
gene_sub = gene[rownames(count_sub),]

pca_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/WT_pca_coor.rds"))
pca_coor = pca_coor[rownames(pd_sub),]
aligned_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_aligned_coor.rds"))
aligned_coor = aligned_coor[rownames(pd_sub),]
umap_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/WT_umap_coor.rds"))
umap_coor = umap_coor[rownames(pd_sub),]

res = my_cluster_cells(umap_coor, pd_sub)

cds <- new_cell_data_set(count_sub,
                         cell_metadata = pd_sub,
                         gene_metadata = gene_sub)

if(ncol(cds) > 50000){
    num_dim = 50
} else if (ncol(cds) > 1000) {
    num_dim = 30
} else {
    num_dim = 10
}

cds <- preprocess_cds(cds, num_dim = num_dim)
reducedDims(cds)$PCA = pca_coor

if (kk == "Mesenchymal trajectory"){
    cds = align_cds(cds = cds,
                    alignment_group = NULL,
                    preprocess_method = "PCA",
                    residual_model_formula_str = "~RIBO_percent + MT_percent + log_umi + g2m_score + s_score")
} else {
    cds = align_cds(cds = cds,
                    alignment_group = NULL,
                    preprocess_method = "PCA",
                    residual_model_formula_str = "~RIBO_percent + MT_percent + log_umi")
}
reducedDims(cds)$Aligned = aligned_coor

cds <- reduce_dimension(cds, 
                        umap.n_neighbors = 15, 
                        umap.min_dist = 0.1,
                        reduction_method = "UMAP", 
                        preprocess_method = "Aligned", 
                        max_components = 3)
reducedDims(cds)$UMAP = umap_coor

cds = cluster_cells(cds, cluster_method = "louvain", reduction_method = "UMAP", random_seed = 123)
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")

if(sum(pData(cds)$my_cluster != res$clusters) == 0){
    saveRDS(cds, paste0(work_path, "/data/sub_trajectory_WT/", name, ".rds"))
} else {
    print("ERROR!")
}


