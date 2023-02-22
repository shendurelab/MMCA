### Section - 3: After removing cells with poor quality and potential doublets, we plan to perform such analysis for
### cell clustering and further annotation

#on the global level
#1. run PCA on WT cells only, and save the PC space;
#2. for each of the mutant type, projecting their cells into the PC space
#3. run align_cds on the embedding to get the aligned PC space.

#########################################################################
### 1 - first step, calculating some scores which will be used for align_cds as covariables

library(monocle3)
library(Matrix)

count = readRDS(paste0(work_path, "/orig_data/gene_count.RDS"))
pd_orig = readRDS(paste0(work_path, "/data/backup/pd_orig.rds"))
print(sum(!rownames(pd_orig) %in% colnames(count)))
count = count[,rownames(pd_orig)]
print(sum(rownames(pd_orig) != colnames(count)))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))
print(sum(rownames(gene) != rownames(count)))

### log_umi
pd_orig$log_umi = log(Matrix::colSums(count))

### RT_percent and Ribo_percent
gene$gene_id = unlist(lapply(as.vector(gene$gene_id), function(x) strsplit(x,"[.]")[[1]][1]))

MT_gene = as.vector(gene[grep("^mt-",gene$gene_short_name),]$gene_id)
Rpl_gene = as.vector(gene[grep("^Rpl",gene$gene_short_name),]$gene_id)
Mrpl_gene = as.vector(gene[grep("^Mrpl",gene$gene_short_name),]$gene_id)
Rps_gene = as.vector(gene[grep("^Rps",gene$gene_short_name),]$gene_id)
Mrps_gene = as.vector(gene[grep("^Mrps",gene$gene_short_name),]$gene_id)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd_orig$MT_percent = 100 * Matrix::colSums(count[gene$gene_id %in% MT_gene, ])/Matrix::colSums(count)
pd_orig$RIBO_percent = 100 * Matrix::colSums(count[gene$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)

### correct s/g2m if necessary
count = t(t(count) / Matrix::colSums(count)) * 100000
count@x = log(count@x + 1)

cellcycle_genes <- readRDS('cellcycle.gene.rds')
s.genes <- cellcycle_genes$s.genes
g2m.genes <- cellcycle_genes$g2m.genes

pd_orig$g2m_score = Matrix::colSums(count[gene$gene_id %in% g2m.genes,])
pd_orig$s_score = Matrix::colSums(count[gene$gene_id %in% s.genes,])

saveRDS(pd_orig, paste0(work_path, "/data/backup/pd_orig.rds"))


### removing cells with high mito%, ribo%, from sample_104, 41
pd = pd_orig[pd_orig$MT_percent <= 10 & pd_orig$RIBO_percent <=5 & !pd_orig$RT_group %in% c("104", "41"),]
count = count[,rownames(pd)]

saveRDS(pd, paste0(work_path, "/data/backup/pd.rds"))
saveRDS(count, paste0(work_path, "/data/backup/count.rds"))

#########################################################################
### 2 - second step, running PCA on WT cells only and save that space

library(monocle3)
library(Matrix)
library(dplyr)

count = readRDS(paste0(work_path, "/data/backup/count.rds"))
pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

count_sub = count[,pd$Mutant == "WT"]
pd_sub = pd[pd$Mutant == "WT",]
gene_sub = gene[rownames(count_sub),]
print(sum(rownames(pd_sub) != colnames(count_sub)))
print(nrow(pd_sub))
print(nrow(pd))

cds <- new_cell_data_set(count_sub,
                         cell_metadata = pd_sub,
                         gene_metadata = gene_sub)

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
fm_rowsums = Matrix::rowSums(FM)
FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]

num_dim = 50
scaling = TRUE
set.seed(2016)
irlba_res <- my_sparse_prcomp_irlba(Matrix::t(FM), 
                                 n = min(num_dim, min(dim(FM)) - 1), 
                                 center = scaling, 
                                 scale. = scaling)
preproc_res <- irlba_res$x
row.names(preproc_res) <- colnames(cds)

saveRDS(gene[rownames(FM),], paste0(work_path, "/data/main_trajectory/gene_use_pca.rds"))
saveRDS(preproc_res, paste0(work_path, "/data/main_trajectory/WT_pca_coor.rds"))
saveRDS(irlba_res, paste0(work_path, "/data/main_trajectory/pca_base.rds"))


### for each mutant type, projecting their cells onto the pca_base ###
library(monocle3)
library(Matrix)
library(dplyr)

count = readRDS(paste0(work_path, "/data/backup/count.rds"))
pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

irlba_res = readRDS(paste0(work_path, "/data/main_trajectory/pca_base.rds"))
gene_use = readRDS(paste0(work_path, "/data/main_trajectory/gene_use_pca.rds"))

pd$Mutant_id = gsub(" ", "_", pd$Mutant)
mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

args = commandArgs(trailingOnly=TRUE)
kk = mutant_list[as.numeric(args[1])]

count_sub = count[,pd$Mutant_id == kk]
pd_sub = pd[pd$Mutant_id == kk,]
gene_sub = gene[rownames(count_sub),]
print(sum(rownames(pd_sub) != colnames(count_sub)))

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

saveRDS(preproc_res_query, paste0(work_path, "/data/main_trajectory/", kk, "_pca_coor.rds"))


#########################################################################
### 3 - third step, running align_cds on the pca_coor to regress out Mito% and Ribo%

library(monocle3)
library(Matrix)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/backup/pd.rds"))

pd$Mutant_id = gsub(" ", "_", pd$Mutant)
mutant_list = as.vector(unique(pd$Mutant_id))

pca_coor = NULL
for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(mutant_i)
    
    tmp = readRDS(paste0(work_path, "/data/main_trajectory/", mutant_i, "_pca_coor.rds"))
    pca_coor = rbind(pca_coor, tmp)
}
print(sum(!rownames(pca_coor) %in% rownames(pd)))
pca_coor = pca_coor[rownames(pd),]
saveRDS(pca_coor, paste0(work_path, "/data/main_trajectory/combined_pca_coor.rds"))

set.seed(2016)
residual_model_formula_str = "~RIBO_percent + MT_percent + log_umi"
X.model_mat <- Matrix::sparse.model.matrix(stats::as.formula(residual_model_formula_str), 
                                           data = pd, drop.unused.levels = TRUE)
fit <- limma::lmFit(Matrix::t(pca_coor), X.model_mat)

beta <- fit$coefficients[, -1, drop = FALSE]
beta[is.na(beta)] <- 0
aligned_coor <- Matrix::t(as.matrix(Matrix::t(pca_coor)) - 
                             beta %*% Matrix::t(X.model_mat[, -1]))

saveRDS(aligned_coor, paste0(work_path, "/data/main_trajectory/combined_aligned_coor.rds"))



