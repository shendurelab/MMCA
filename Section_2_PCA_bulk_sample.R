
##########################################
### make PCA plot on individual sample ###
##########################################

library(monocle3)
library(ggplot2)

work_path = "~/work/mutant"
count = readRDS(paste0(work_path, "/bulk_sample/embryo.rds"))
pd = readRDS(paste0(work_path, "/bulk_sample/pd_mutant_embryo.rds"))
pd$Mutant = gsub(" ", "_", pd$Mutant)
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

color_plate_read = read.table(paste0(work_path, "/code/color_plate.txt"), as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

pd = pd[colnames(count),]
gene = gene[rownames(count),]

cds = new_cell_data_set(count,
                        cell_metadata = pd,
                        gene_metadata = gene)

num_dim = 15
cds = preprocess_cds(cds, num_dim = num_dim)
print(cds@preprocess_aux$prop_var_expl)
PCA_coor = data.frame(reducedDims(cds)$PCA[,1:3]); names(PCA_coor) = c('PC_1','PC_2','PC_3')
pd = cbind(pd, PCA_coor)


library(htmlwidgets)
library(plotly)
library(dplyr)

fig = plot_ly(pd, x=~PC_1, y=~PC_2, z=~PC_3, color = ~Background, colors=color_plate) %>%
    layout(scene =list(xaxis = list(title = "PC_1 (39.0%)"),
           yaxis = list(title = "PC_2 (16.7%)"),
           zaxis = list(title = "PC_3 (8.4%)")))
saveWidget(fig, paste0(work_path, "/bulk_sample/PCA_mutant.html"))


df = data.frame(x = rep(1, length(color_plate)),
                y = rep(1, length(color_plate)),
                group = factor(names(color_plate), levels = names(color_plate)))
p = ggplot(df, aes(x,y,color=group)) + 
    geom_point(size = 5) + 
    scale_colour_manual(values = color_plate) 




##############################
### projecting PCA on MOCA ###
##############################


library(monocle3)
library(ggplot2)
library(plotly)
library(htmlwidgets)
source("~/work/scripts/tome/utils.R")

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/mutant"
count = readRDS(paste0(work_path, "/bulk_sample/embryo.rds"))
pd = readRDS(paste0(work_path, "/bulk_sample/pd_mutant_embryo.rds"))
pd$Mutant = gsub(" ", "_", pd$Mutant)
gene = readRDS(paste0(work_path, "/orig_data/df_gene.RDS"))

pd = pd[colnames(count),]
gene = gene[rownames(count),]

pd$development_stage = "Mutant"
pd$project = "Mutant"

rownames(count) = unlist(lapply(rownames(count), function(x) strsplit(x,"[.]")[[1]][1]))
rownames(gene) = unlist(lapply(rownames(gene), function(x) strsplit(x,"[.]")[[1]][1]))

cds_mutant = new_cell_data_set(count,
                               cell_metadata = pd,
                               gene_metadata = gene)

count = readRDS(paste0(work_path, "/bulk_sample/moca/embryo.rds"))
pd = readRDS(paste0(work_path, "/bulk_sample/moca/pd_moca_embryo.rds"))
gene = readRDS(paste0(work_path, "/bulk_sample/moca/gene_moca.rds"))
pd = pd[colnames(count),]
gene = gene[rownames(count),]

rownames(count) = unlist(lapply(rownames(count), function(x) strsplit(x,"[.]")[[1]][1]))
rownames(gene) = unlist(lapply(rownames(gene), function(x) strsplit(x,"[.]")[[1]][1]))

moca = new_cell_data_set(count,
                         cell_metadata = pd,
                         gene_metadata = gene)

pData(moca)$project = "MOCA"

set.seed(2016)
FM = monocle3:::normalize_expr_data(moca, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
fm_rowsums = Matrix::rowSums(FM)
FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
FM = FM[rownames(FM) %in% rownames(cds_mutant),]
gene_use = rownames(FM)

num_dim = 15
scaling = TRUE
set.seed(2016)
irlba_res <- my_sparse_prcomp_irlba(Matrix::t(FM), 
                                    n = min(num_dim, min(dim(FM)) - 1), 
                                    center = scaling, 
                                    scale. = scaling)
preproc_res <- irlba_res$x
row.names(preproc_res) <- colnames(moca)

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds_mutant, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
FM = FM[gene_use,]
preproc_res_query <- scale(t(FM), irlba_res$center, irlba_res$scale) %*% irlba_res$rotation
row.names(preproc_res_query) <- colnames(cds_mutant)

my_cols = c("9.5"="#F8766D",
            "10.5"="#A3A500",
            "11.5"="#00BF7D",
            "12.5"="#00B0F6",
            "13.5"="#E76BF3",
            "Mutant"="darkgrey")

library(dplyr)
library(ggplot2)
df = data.frame(ID = c(rownames(preproc_res), rownames(preproc_res_query)),
                PC_1 = c(preproc_res[,1], preproc_res_query[,1]),
                PC_2 = c(preproc_res[,2], preproc_res_query[,2]),
                PC_3 = c(preproc_res[,3], preproc_res_query[,3]),
                development_stage = c(pData(moca)$development_stage, pData(cds_mutant)$development_stage))
df$development_stage = factor(df$development_stage, levels = c("9.5","10.5","11.5","12.5","13.5","Mutant"))

p <- ggplot(df, aes(PC_1,PC_2,color=development_stage,label=ID)) + geom_point(size=2)
p = p + theme_classic(base_size = 12) + scale_color_manual(values=my_cols) + labs(x = "PC_1 (54.3%)", y = "PC_2 (12.3%)")
p = p + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 



