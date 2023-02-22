### processing the whole data set to obtain the main-trajectory information ###

#############################################################
### 1: process the data set to obatin the main-trajectory ###
#############################################################

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

### 1,685,704 cells, which have filtered out cells labeled as doublets (by Scrublet) or from doublet-derived subclusters

gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
count = readRDS(paste0(work_path, "/orig_data/gene_count.RDS"))

mouse_gene = read.table("mouse.geneID.loc.txt", sep="\t", as.is=T)
mouse_gene_sex = as.vector(mouse_gene$V5[mouse_gene$V1 %in% c("chrX", "chrY")])
gene = gene[gene$gene_type %in% c("protein_coding","lincRNA","pseudogene"),]

count_sub = count[!rownames(count) %in% mouse_gene_sex,]
count_sub = count_sub[rownames(count_sub) %in% rownames(gene),]

pd_seurat = readRDS(paste0(work_path, "/seurat/pd_all_exclude_cluster.rds"))
pd_seurat = pd_seurat[pd_seurat$exclude_cluster | pd_seurat$Doublet == "doublet",]
count_sub = count_sub[,!colnames(count_sub) %in% rownames(pd_seurat)]

obj = CreateSeuratObject(count_sub, assay = "RNA", min.cells=10, min.features=100)
print(dim(obj))

saveRDS(GetAssayData(object = obj, slot = "counts"), paste0(work_path, "/data/count.rds"))

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
print("Done Normalization")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
print("Done Selecting highly variable genes")

saveRDS(VariableFeatures(obj), paste0(work_path, "/data/top_5000_genes.rds"))

### create cds file
pd = readRDS(paste0(work_path, "/orig_data/df_cell.RDS"))
rownames(pd) = as.vector(pd$sample)
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
count = readRDS(paste0(work_path, "/data/count.rds"))
include_gene = readRDS(paste0(work_path, "/data/top_5000_genes.rds"))

pd = pd[colnames(count),]
gene = gene[rownames(count),]

cds <- new_cell_data_set(count,
                         cell_metadata = pd,
                         gene_metadata = gene)

cds <- preprocess_cds(cds, num_dim = 50, use_genes = include_gene)
cds <- reduce_dimension(cds, umap.n_neighbors=50)
cds = cluster_cells(cds, resolution=1e-6, cluster_method = "leiden")
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")

saveRDS(cds, file=paste0(work_path, "/data/cds_all.RDS"))
saveRDS(data.frame(pData(cds)), file=paste0(work_path, "/data/pd_all_1e-06.RDS"))

jpeg(paste0(work_path, "/data/plot/", "cluster_umap", ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_cluster", group_cells_by="cluster", group_label_size = 4))
dev.off()

pData(cds)$my_partition = partitions(cds, reduction_method = "UMAP")
jpeg(paste0(work_path, "/data/plot/", "partition_umap", ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_partition", group_cells_by="partition", group_label_size = 4))
dev.off()

pd_old = readRDS(file=paste0(work_path, "/backup/pd_all.rds"))
trajectory_list = names(table(pd_old$update_trajectory_name))
for(i in 1:length(trajectory_list)){
  print(trajectory_list[i])
  pd_tmp = pd_old[pd_old$update_trajectory_name == trajectory_list[i],]
  pData(cds)$tmp = rownames(pd) %in% rownames(pd_tmp)
  jpeg(paste0(work_path, "/data/plot/", trajectory_list[i], ".jpeg"), width=12, height=12, units = 'in', res=300)
  print(plot_cells(cds, color_cells_by="tmp", group_cells_by="cluster", group_label_size = 4))
  dev.off()
}

main_cluster = as.vector(pData(cds)$my_partition)
main_cluster[pData(cds)$my_cluster == 9] = 11
main_cluster[pData(cds)$my_cluster == 19] = 12
pData(cds)$main_cluster = main_cluster
saveRDS(data.frame(pData(cds)), file=paste0(work_path, "/data/pd_all_1e-06.RDS"))

############ find marker genes of each main cluster ###############

library(monocle3)
library(dplyr)

cds = readRDS(paste0(work_path, "/data/cds_all.RDS"))
pd = readRDS(paste0(work_path, "/data/pd_all_1e-06.RDS"))
pd = pd[rownames(pData(cds)),]
pData(cds)$main_cluster = pd$main_cluster

marker_test_res <- top_markers(cds, group_cells_by="main_cluster", 
                               reference_cells=1000, cores=8)

markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

saveRDS(markers, file=paste0(work_path, "/data/marker_top10.rds"))


############################################################
### 2: dig deeper each main_cluster to identify doublets ###
############################################################

# To further remove such doublet cells, we took the cell clusters identified by Monocle 3 and 
# first computed differentially expressed genes across cell clusters with the top_markers() 
# function of Monocle 3. We then selected a gene set combining the top ten gene markers for 
# each cell cluster (ordered by pseudo_R2).

# Subclusters showing low expression of target cell cluster-specific markers and enriched 
# expression of non-target cell cluster-specific markers were annotated as doublets derived 
# subclusters and filtered out in visualization and downstream analysis (1,671,270 cells were left). 


args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
library(monocle3)
library(dplyr)
library(ggplot2)

pd_all = readRDS(file=paste0(work_path, "/data/pd_all_1e-06.RDS"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
count = readRDS(paste0(work_path, "/data/count.rds"))
markers = readRDS(file=paste0(work_path, "/data/marker_top10.rds"))

pd = pd_all[pd_all$main_cluster == kk,]
count = count[,colnames(count) %in% rownames(pd)]
pd = pd[colnames(count),]
gene = gene[rownames(count),]

cds <- new_cell_data_set(count,
                         cell_metadata = pd,
                         gene_metadata = gene)

saveRDS(cds, file=paste0(work_path, "/data/", kk, ".RDS"))

cds <- preprocess_cds(cds, num_dim = 10, use_genes = as.vector(markers$gene_id))
cds <- reduce_dimension(cds, umap.n_neighbors=50)
cds = cluster_cells(cds, resolution=1e-4, cluster_method = "leiden")
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")

saveRDS(cds, file=paste0(work_path, "/data/", kk, ".RDS"))

jpeg(paste0(work_path, "/data/plot/", kk, ".jpeg"), width=10, height=10, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_cluster", group_cells_by="cluster", group_label_size = 4))
dev.off()

for(j in 1:12){
  marker_sub = markers %>% filter(cell_group==j)
  jpeg(paste0(work_path, "/data/plot/", kk, "_", j, ".jpeg"), width=15, height=15, units = 'in', res=300)
  print(plot_cells(cds, genes=marker_sub$gene_short_name))
  dev.off()
}

df = data.frame(pData(cds))
jpeg(paste0(work_path, "/data/plot/", kk, "_boxplot.jpeg"), width=10, height=10, units = 'in', res=300)
print(ggplot(df, aes(factor(my_cluster), doublet_score)) + geom_boxplot())
dev.off()


### individual main cluster ###
library(monocle3)
library(dplyr)
library(ggplot2)

kk = 2
cds = readRDS(file=paste0(work_path, "/data/", kk, ".RDS"))
table(pData(cds)$my_cluster)

jpeg(paste0(work_path, "/data/plot/", kk, "_doublet_score.jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="doublet_score", group_cells_by="cluster", group_label_size = 4))
dev.off()

cds = cluster_cells(cds, resolution=1e-3, random_seed = 1)
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")
table(pData(cds)$my_cluster)

jpeg(paste0(work_path, "/data/plot/", kk, ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_cluster", group_cells_by="cluster", group_label_size = 4))
dev.off()

df = data.frame(pData(cds))
jpeg(paste0(work_path, "/data/plot/", kk, "_boxplot.jpeg"), width=12, height=12, units = 'in', res=300)
print(ggplot(df, aes(factor(my_cluster), doublet_score)) + geom_boxplot())
dev.off()

exclude = NULL

pData(cds)$exclude_subcluster = pData(cds)$my_cluster %in% exclude
res = data.frame(pData(cds))[,c("main_cluster","my_cluster","exclude_subcluster")]

saveRDS(res, paste0(work_path, "/data/pd_", kk, ".RDS"))

pd = NULL
for(kk in 1:12){
  x = readRDS(paste0(work_path, "/data/pd_", kk, ".RDS"))  
  pd = rbind(pd, x)
}

saveRDS(pd, paste0(work_path, "/data/pd_removeDoublets.RDS"))



#######################################
### 3: finally find main trajectory ###
#######################################

# Finally, genes expressed in less than 10 cells and cells expressing less than 100 genes 
# were further filtered out (1,671,269 cells vs. 22,050 genes were left). 
# The downstream dimension reduction and clustering analysis were done by Monocle 3. 
# The dimensionality of the data was reduced by PCA (50 components) first and 
# then with UMAP (max_components = 2, n_neighbors = 50, min_dist = 0.1, metric = 'cosine'). 
# Cell clusters were identified using the Leiden algorithm implemented in Monocle 3 (resolution = 1e-06). 
# 15 isolated partitions were identified.

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)

mouse_gene = read.table("mouse.geneID.loc.txt", sep="\t", as.is=T)
mouse_gene_sex = as.vector(mouse_gene$V5[mouse_gene$V1 %in% c("chrX", "chrY")])

gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
gene = gene[gene$gene_type %in% c("protein_coding","lincRNA","pseudogene"),]
gene = gene[!rownames(gene) %in% mouse_gene_sex,]
pd_removeDoublets = readRDS(paste0(work_path, "/remove_doublets/step2_removeByMonocle.rds"))

count = readRDS(paste0(work_path, "/orig_data/gene_count.RDS"))
count_sub = count[rownames(count) %in% rownames(gene), colnames(count) %in% rownames(pd_removeDoublets)[!pd_removeDoublets$exclude_subcluster]]

obj = CreateSeuratObject(count_sub, assay = "RNA", min.cells=10, min.features=100)
print(dim(obj))

saveRDS(GetAssayData(object = obj, slot = "counts"), paste0(work_path, "/data/main_trajectory/count.rds"))
print("Done write count matrix")

### create CDS object

count = readRDS(paste0(work_path, "/data/main_trajectory/count.rds"))

gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
pd = readRDS(paste0(work_path, "/data/backup/pd_all.RDS"))
rownames(pd) = as.vector(pd$sample)

pd = pd[colnames(count),]
gene = gene[rownames(count),]

cds <- new_cell_data_set(count,
                         cell_metadata = pd,
                         gene_metadata = gene)

cds <- preprocess_cds(cds, num_dim = 50)

colData(cds)$log_n.umi = 
  log(Matrix::colSums(count))

fData(cds)$gene_id = unlist(lapply(as.vector(fData(cds)$gene_id), function(x) strsplit(x,"[.]")[[1]][1]))

count = t(t(count) / Matrix::colSums(count)) * 100000
count@x = log(count@x + 1)

MT_gene = as.vector(mouse_gene[grep("^mt",mouse_gene$gene_short_name),]$ID)
Rpl_gene = as.vector(mouse_gene[grep("^Rpl",mouse_gene$gene_short_name),]$ID)
Mrpl_gene = as.vector(mouse_gene[grep("^Mrpl",mouse_gene$gene_short_name),]$ID)
Rps_gene = as.vector(mouse_gene[grep("^Rps",mouse_gene$gene_short_name),]$ID)
Mrps_gene = as.vector(mouse_gene[grep("^Mrps",mouse_gene$gene_short_name),]$ID)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

colData(cds)$MT_percent = Matrix::colSums(count[fData(cds)$gene_id %in% MT_gene, ])/Matrix::colSums(count)
colData(cds)$RIBO_percent = Matrix::colSums(count[fData(cds)$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)

cds = align_cds(cds = cds,
            alignment_group = NULL,
            preprocess_method = "PCA",
            residual_model_formula_str = "~log_n.umi + MT_percent + RIBO_percent")

cds <- reduce_dimension(cds, umap.n_neighbors=50)
cds = cluster_cells(cds, resolution=1e-6, cluster_method = "leiden")
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")
pData(cds)$my_partition = partitions(cds, reduction_method = "UMAP")

saveRDS(cds, file=paste0(work_path, "/data/main_trajectory/cds_all.RDS"))
saveRDS(data.frame(pData(cds)), file=paste0(work_path, "/data/main_trajectory/pd_all.RDS"))

tmp = as.vector(pData(cds)$my_cluster)
names(tmp) = colnames(exprs(cds))
cds@clusters[["UMAP"]][['clusters']] = factor(tmp)

jpeg(paste0(work_path, "/data/main_trajectory/plot/", "cluster_umap", ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_cluster", group_cells_by="cluster", group_label_size = 4))
dev.off()

tmp = as.vector(pData(cds)$my_partition)
names(tmp) = colnames(exprs(cds))
cds@clusters[["UMAP"]][['clusters']] = factor(tmp)

jpeg(paste0(work_path, "/data/main_trajectory/plot/", "partition_umap", ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="my_partition", group_cells_by="cluster", group_label_size = 4))
dev.off()

jpeg(paste0(work_path, "/data/main_trajectory/plot/", "doublet_score.jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="doublet_score", group_cells_by="cluster", group_label_size = 4))
dev.off()

pd_old = readRDS(file=paste0(work_path, "/data/backup/pd_all.RDS"))
trajectory_list = names(table(pd_old$main_trajectory))
for(i in 1:length(trajectory_list)){
  print(trajectory_list[i])
  pd_tmp = pd_old[pd_old$main_trajectory == trajectory_list[i],]
  pData(cds)$tmp = rownames(pd) %in% rownames(pd_tmp)
  jpeg(paste0(work_path, "/data/main_trajectory/plot/", trajectory_list[i], ".jpeg"), width=12, height=12, units = 'in', res=300)
  print(plot_cells(cds, color_cells_by="tmp", group_cells_by="cluster", group_label_size = 4))
  dev.off()
}

### show markers of each moca main trajectory
library(monocle3)
library(dplyr)
library(ggplot2)

cds = readRDS(file=paste0(work_path, "/data/main_trajectory/cds_all.RDS"))
moca_marker = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/data/sci3_trajectory/main_trajectory_summary_monocle3/cds_moca_marker_top10.rds")

trajectory_list = as.vector(unique(moca_marker$cell_group))
for(i in 1:length(trajectory_list)){
  print(trajectory_list[i])
  gene = as.vector(moca_marker$gene_short_name[moca_marker$cell_group==trajectory_list[i]])
  for(j in 1:length(gene)){
    print(paste0(j, "/", length(gene)))
    if(gene[j] %in% fData(cds)$gene_short_name){
      jpeg(paste0(work_path, "/data/main_trajectory/plot/", trajectory_list[i], "_", gene[j], ".jpeg"), width=12, height=12, units = 'in', res=300)
      print(plot_cells(cds, genes = gene[j], group_label_size = 4))
      dev.off()
    }
  }
}


### set main trajectory ###
cds = readRDS(file=paste0(work_path, "/data/main_trajectory/cds_all.RDS"))
pd = readRDS(paste0(work_path, "/data/main_trajectory/pd_all.RDS"))

main_trajectory = rep(NA, nrow(pd))
main_trajectory[pd$my_partition %in% c(1)] = "1_neural_tube"
main_trajectory[pd$my_partition %in% c(2)] = "2_mesenchyme"
main_trajectory[pd$my_partition %in% c(3)] = "3_epithelial"
main_trajectory[pd$my_partition %in% c(4)] = "4_myocytes"
main_trajectory[pd$my_partition %in% c(5)] = "5_blood"
main_trajectory[pd$my_partition %in% c(6)] = "6_endothelial"
main_trajectory[pd$my_partition %in% c(7)] = "7_neural_crest_neuron"
main_trajectory[pd$my_partition %in% c(8)] = "8_liver"
main_trajectory[pd$my_partition %in% c(9)] = "9_neural_crest_glia"
main_trajectory[pd$my_partition %in% c(10)] = "10_white_blood_1"
main_trajectory[pd$my_partition %in% c(11)] = "11_white_blood_2"
main_trajectory[pd$my_partition %in% c(12)] = "12_neural_crest_melanocyte"
main_trajectory[pd$my_partition %in% c(13)] = "13_megakaryocyte"
main_trajectory[pd$my_partition %in% c(14)] = "14_cardiac_muscle"
main_trajectory[pd$my_partition %in% c(15)] = "15_len"

pData(cds)$main_trajectory = main_trajectory
jpeg(paste0(work_path, "/data/main_trajectory/plot/", "main_trajectory", ".jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_cells(cds, color_cells_by="main_trajectory", group_cells_by="partition", group_label_size = 4))
dev.off()

pd$main_trajectory = main_trajectory
saveRDS(pd, file=paste0(work_path, "/data/main_trajectory/pd_all.RDS"))


##########################################################
### 4: identify key marker genes for main-trajectories ###
##########################################################

library(monocle3)
library(dplyr)
library(ggplot2)

cds = readRDS(file=paste0(work_path, "/data/main_trajectory/cds_all.RDS"))
pData(cds)$my_cluster = clusters(cds, reduction_method = "UMAP")
pData(cds)$my_partition = partitions(cds, reduction_method = "UMAP")

marker_test_res <- top_markers(cds, group_cells_by="my_partition", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(2, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

jpeg(paste0(work_path, "/data/main_trajectory/plot/cds_main_trajectory_marker_top2.jpeg"), width=12, height=12, units = 'in', res=300)
print(plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="my_partition",
                    ordering_type="cluster_row_col",
                    max.size=3))
dev.off()

markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

saveRDS(markers, file=paste0(work_path, "/data/main_trajectory/cds_main_trajectory_marker_top10.rds"))








