### on the global level

### 4, doing UMAP on WT and save the space using uwot; projecting cells from each mutant type to this space and get their umap coordinates
### 5, doing louain clustering on WT cells based on 3d UMAP embedding
### 6, performing knn on each mutant type to annotate their main trajectories

#####################################
### step 4, doing UMAP on WT and save the space using uwot

library(monocle3)
library(Matrix)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
aligned_coor = readRDS(paste0(work_path, "/data/main_trajectory/combined_aligned_coor.rds"))
print(sum(rownames(pd) != rownames(aligned_coor)))

pd_sub = pd[pd$Mutant == "WT",]
aligned_coor_sub = aligned_coor[rownames(pd_sub),]


set.seed(2016)
emb_train_umap = uwot::umap(as.matrix(aligned_coor_sub), 
                            n_components = 3,
                            n_neighbors = 50,
                            min_dist = 0.01,
                            metric = "cosine",
                            fast_sgd = FALSE,
                            nn_method = "annoy",
                            ret_model = TRUE,
                            n_threads = 1,
                            verbose = TRUE)

set.seed(2016)
emb_test_coor = uwot::umap_transform(as.matrix(aligned_coor_sub),
                                     emb_train_umap)

umap_coor = emb_train_umap$embedding
rownames(umap_coor) = rownames(pd_sub)
colnames(umap_coor) = paste0("UMAP_", 1:3)

saveRDS(umap_coor, paste0(work_path, "/data/main_trajectory/WT_umap_coor.rds"))
uwot::save_uwot(emb_train_umap, paste0(work_path, "/data/main_trajectory/umap_base"))

###################
### projecting cells from each mutant type and get their umap coordinates

library(monocle3)
library(Matrix)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
aligned_coor = readRDS(paste0(work_path, "/data/main_trajectory/combined_aligned_coor.rds"))
print(sum(rownames(pd) != rownames(aligned_coor)))

#emb_train_umap = uwot::load_uwot(paste0(work_path, "/data/main_trajectory/umap_base"))

pd$Mutant_id = gsub(" ", "_", pd$Mutant)
mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

for(i in 1:length(mutant_list)){
    kk = mutant_list[i]
    print(paste0(i, "/", kk))
    
    pd_sub = pd[pd$Mutant_id == kk,]
    aligned_coor_sub = aligned_coor[rownames(pd_sub),]
    
    set.seed(2016)
    emb_test_coor = uwot::umap_transform(as.matrix(aligned_coor_sub),
                                         emb_train_umap)
    
    rownames(emb_test_coor) = rownames(pd_sub)
    colnames(emb_test_coor) = paste0("UMAP_", 1:3)
    
    saveRDS(emb_test_coor, paste0(work_path, paste0("/data/main_trajectory/", kk, "_umap_coor.rds")))

}


### combining all the mutants as well as WT
mutant_list = as.vector(unique(pd$Mutant_id))
umap_coor = NULL
for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(mutant_i)
    
    tmp = readRDS(paste0(work_path, "/data/main_trajectory/", mutant_i, "_umap_coor.rds"))
    umap_coor = rbind(umap_coor, tmp)
}
print(sum(!rownames(umap_coor) %in% rownames(pd)))
umap_coor = umap_coor[rownames(pd),]
saveRDS(umap_coor, paste0(work_path, "/data/main_trajectory/combined_umap_coor.rds"))



######################
### making the 3d plots
library(htmlwidgets)
library(plotly)

pd = readRDS(paste0(work_path, "/data/main_trajectory/WT_umap_coor.rds"))
pd = data.frame(pd)
pd_update = readRDS(paste0(work_path, "/data/backup/pd_update.rds"))
pd_update = pd_update[rownames(pd),]
pd$main_trajectory = as.vector(pd_update$main_trajectory)
pd$sub_trajectory = as.vector(pd_update$sub_trajectory)

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(25), color = ~main_trajectory)

fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))

saveWidget(fig, paste0(work_path, "/data/main_trajectory/plot/", "global_WT_main_trajectory.html"))

fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(25), color = ~sub_trajectory)

fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))

saveWidget(fig, paste0(work_path, "/data/main_trajectory/plot/", "global_WT_sub_trajectory.html"))



#####################################
### step 5, doing louain clustering on WT cells based on 3d UMAP embedding
### and then assigning the main trajectories

rm(list = ls())
library(monocle3)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
umap_coor = readRDS(paste0(work_path, "/data/main_trajectory/WT_umap_coor.rds"))
pd_sub = pd[rownames(umap_coor),]

res = my_cluster_cells(umap_coor, pd_sub)

pd_sub$my_cluster = as.vector(res$clusters)
pd_sub$my_partition = as.vector(res$partitions)
pd_sub$UMAP_1 = umap_coor[,1]
pd_sub$UMAP_2 = umap_coor[,2]
pd_sub$UMAP_3 = umap_coor[,3]

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(25), color = ~my_cluster)

fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))

saveWidget(fig, paste0(work_path, "/data/main_trajectory/plot/", "global_WT_my_cluster.html"))

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(25), color = ~my_partition)

fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))

saveWidget(fig, paste0(work_path, "/data/main_trajectory/plot/", "global_WT_my_partition.html"))

main_trajectory = rep(NA, nrow(pd_sub))
main_trajectory[pd_sub$my_partition == 1] = "Neural tube and notochord trajectory"
main_trajectory[pd_sub$my_partition == 2] = "Endothelial trajectory"
main_trajectory[pd_sub$my_partition == 3] = "Haematopoiesis trajectory"
main_trajectory[pd_sub$my_partition == 4] = "Myotube trajectory"
main_trajectory[pd_sub$my_partition == 5] = "Mesenchymal trajectory"
main_trajectory[pd_sub$my_partition == 6] = "Neural crest (PNS neuron) trajectory"
main_trajectory[pd_sub$my_partition == 7] = "Myoblast trajectory"
main_trajectory[pd_sub$my_partition == 8] = "Epithelial trajectory"
main_trajectory[pd_sub$my_partition == 9] = "Hepatocyte trajectory"
main_trajectory[pd_sub$my_partition == 10] = "Melanocyte trajectory"
main_trajectory[pd_sub$my_partition == 11] = "Neural crest (PNS glia) trajectory"
main_trajectory[pd_sub$my_partition == 12] = "Olfactory sensory neuron trajectory"
main_trajectory[pd_sub$my_partition == 13] = "Cardiomyocyte trajectory"

pd_sub$main_trajectory = as.vector(main_trajectory)

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(25), color = ~main_trajectory)

fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))

saveWidget(fig, paste0(work_path, "/data/main_trajectory/plot/", "global_WT_main_trajectory.html"))

saveRDS(pd_sub, paste0(work_path, "/data/main_trajectory/WT_pd.rds"))




##################################################################################
### step 6, performing knn on each mutant type to annotate their main trajectories

library(monocle3)
library(Matrix)
library(dplyr)
library(FNN)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
pd$Mutant_id = gsub(" ", "_", pd$Mutant)
mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

WT_pd = readRDS(paste0(work_path, "/data/main_trajectory/WT_pd.rds"))
WT_coor = readRDS(paste0(work_path, "/data/main_trajectory/WT_umap_coor.rds"))
print(sum(rownames(WT_pd) != rownames(WT_coor)))
y_train = as.vector(WT_pd$main_trajectory)

k_num = 15

for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(paste0(i, "/", mutant_i))
    
    Mutant_pd = pd[pd$Mutant_id == mutant_i,]
    Mutant_coor = readRDS(paste0(work_path, "/data/main_trajectory/", mutant_i, "_umap_coor.rds"))
    err_num = sum(!rownames(Mutant_pd) %in% rownames(Mutant_coor))
    if(err_num != 0){
        break
    }
    Mutant_pd = Mutant_pd[rownames(Mutant_coor),]
    
    neighbors <- get.knnx(WT_coor, Mutant_coor, k = k_num)$nn.index
    
    y_test = rep(NA, nrow(Mutant_pd))
    y_score = rep(NA, nrow(Mutant_pd))
    for(j in 1:nrow(neighbors)){
        tmp = table(y_train[neighbors[j,]])
        y_test[j] = names(tmp)[which.max(tmp)]
        y_score[j] = max(tmp)/k_num
    }
    
    Mutant_pd$main_trajectory = as.vector(y_test)
    Mutant_pd$main_trajectory_score = as.vector(y_score)
    
    saveRDS(Mutant_pd, paste0(work_path, "/data/main_trajectory/", mutant_i, "_pd.rds"))
    
}


#########
### combined each pd to create the combined pd
### and make the UMAP on downsampling cells

library(monocle3)
library(Matrix)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
pd$Mutant_id = gsub(" ", "_", pd$Mutant)
mutant_list = as.vector(unique(pd$Mutant_id))

pd_all = NULL
for(i in 1:length(mutant_list)){
    print(paste0(i, "/", mutant_list[i]))
    pd_tmp = readRDS(paste0(work_path, "/data/main_trajectory/", mutant_list[i], "_pd.rds"))
    if (mutant_list[i] == "WT"){
        pd_tmp$main_trajectory_score = 1
    }
    pd_all = rbind(pd_all, pd_tmp[,c("Mutant","main_trajectory","main_trajectory_score")])
}
pd_all = pd_all[rownames(pd),]
pd$main_trajectory = as.vector(pd_all$main_trajectory)
pd$main_trajectory_score = as.vector(pd_all$main_trajectory_score)

saveRDS(pd, paste0(work_path, "/data/main_trajectory/combined_pd.rds"))


#############################
### the top marker gene of each main_trajectory calculated by WT cells

library(monocle3)
library(Matrix)
library(dplyr)

count = readRDS(paste0(work_path, "/data/backup/count.rds"))
pd = readRDS(paste0(work_path, "/data/main_trajectory/WT_pd.rds"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

count_sub = count[,rownames(pd)]
gene_sub = gene[rownames(count_sub),]

cds <- new_cell_data_set(count_sub,
                         cell_metadata = pd,
                         gene_metadata = gene_sub)


marker_test_res <- top_markers(cds, group_cells_by="main_trajectory", 
                               reference_cells=1000, cores=8)

markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(10, pseudo_R2)

saveRDS(marker_test_res, paste0(work_path, "/data/main_trajectory/WT_top_marker.rds"))


library(Seurat)
obj = CreateSeuratObject(count_sub, meta.data = pd)
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) = as.vector(obj$main_trajectory)

res = FindMarkers(obj, ident.1 = "Myocyte trajectory 1", ident.2 = "Myocyte trajectory 2", only.pos = T)
res1 = res %>% 
    mutate(ID = unlist(lapply(rownames(res), function(x) strsplit(x,"[.]")[[1]][1]))) %>%
    left_join(mouse_gene, by="ID")

res = FindMarkers(obj, ident.1 = "Myocytes trajectory 2", ident.2 = "Myocytes trajectory 1", only.pos = T)
res2 = res %>% 
    mutate(ID = unlist(lapply(rownames(res), function(x) strsplit(x,"[.]")[[1]][1]))) %>%
    left_join(mouse_gene, by="ID")
