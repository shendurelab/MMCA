### heatmap of nnls with moca

rm(list=ls())
work_path = "~/work/mutant"
source("~/work/tome/code/utils.R")
library(dplyr)
library(tidyr)
library("gplots")
library("reshape2")
library(viridis)
source("~/work/tome/code/backup/Jun_correlation.R")

###########################
### main trajectory

moca_main_trajectory = readRDS(paste0(work_path, "/bulk_sample/moca/main_trajectory.rds"))
moca_main_trajectory_names = colnames(moca_main_trajectory)
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 1"] = "Neural crest (PNS neuron) trajectory 1"
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 2"] = "Neural crest (PNS glia) trajectory 2"
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 3"] = "Neural crest (Melanocyte) trajectory 3"
colnames(moca_main_trajectory) = moca_main_trajectory_names

mutant_main_trajectory = readRDS(paste0(work_path, "/bulk_sample/main_trajectory.rds"))
colnames(mutant_main_trajectory)[colnames(mutant_main_trajectory) == "Neural Crest (Ganglia) trajectory"] = "Olfactory sensory neuron trajectory"

gene_overlap = intersect(rownames(moca_main_trajectory), rownames(mutant_main_trajectory))
moca_main_trajectory_sub = moca_main_trajectory[gene_overlap,]
mutant_main_trajectory_sub = mutant_main_trajectory[gene_overlap,]

conn <- correlation_analysis_bidirection(as.matrix(mutant_main_trajectory_sub), as.matrix(moca_main_trajectory_sub), fold.change = 1.5, top_gene_num = 3000, spec_gene_num = 3000)
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

pdf("~/work/mutant/bulk_sample/nnls_moca/nnls_main_trajectory.pdf", 8, 8)
heatmap.2(as.matrix(dat), 
          col = viridis,
          scale="none", 
          Rowv = TRUE, 
          Colv = TRUE, 
          keysize = 1.5,
          lhei=c(1,3.5), 
          lwid=c(1,3.5),
          density.info = "none", 
          trace = "none",
          cexRow = 1, 
          cexCol = 1,
          margins = c(16, 16))
dev.off()


###########################
### main trajectory subset

moca_main_trajectory = readRDS(paste0(work_path, "/bulk_sample/moca/main_trajectory.rds"))
moca_main_trajectory_names = colnames(moca_main_trajectory)
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 1"] = "Neural crest (PNS neuron) trajectory 1"
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 2"] = "Neural crest (PNS glia) trajectory 2"
moca_main_trajectory_names[moca_main_trajectory_names == "Neural crest 3"] = "Neural crest (Melanocyte) trajectory 3"
colnames(moca_main_trajectory) = moca_main_trajectory_names

mutant_main_trajectory = readRDS(paste0(work_path, "/bulk_sample/main_trajectory_subset.rds"))
colnames(mutant_main_trajectory)[colnames(mutant_main_trajectory) == "Neural Crest (Ganglia) trajectory"] = "Olfactory sensory neuron trajectory"

gene_overlap = intersect(rownames(moca_main_trajectory), rownames(mutant_main_trajectory))
moca_main_trajectory_sub = moca_main_trajectory[gene_overlap,]
mutant_main_trajectory_sub = mutant_main_trajectory[gene_overlap,]

conn <- correlation_analysis_bidirection(as.matrix(mutant_main_trajectory_sub), as.matrix(moca_main_trajectory_sub), fold.change = 1.5, top_gene_num = 3000, spec_gene_num = 3000)
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

dat = dat[rev(c("Endothelial trajectory",
            "Epithelial trajectory",
            "Haematopoiesis trajectory",
            "Liver hepatocyte trajectory",
            "Len epithelial trajectory",
            "Mesenchymal/Myotube/Myoblast/Cardiomyocyte trajectory",
            "Melanocyte trajectory",
            "Neural crest (PNS glia) trajectory",
            "Neural crest (PNS neuron)/Olfactory sensory neuron trajectory",
            "Neural tube and notochord trajectory")),]

pdf("~/work/mutant/bulk_sample/nnls_moca/nnls_main_trajectory_subset.pdf", 8, 8)
heatmap.2(as.matrix(dat), 
          col = viridis,
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          keysize = 1.5,
          lhei=c(1,3.5), 
          lwid=c(1,3.5),
          density.info = "none", 
          trace = "none",
          cexRow = 1, 
          cexCol = 1,
          margins = c(16, 16))
dev.off()





######################
### sub trajectory
moca_sub_trajectory = readRDS(paste0(work_path, "/bulk_sample/moca/sub_trajectory.rds"))
moca_sub_trajectory_names = colnames(moca_sub_trajectory)
colnames(moca_sub_trajectory) = moca_sub_trajectory_names

mutant_sub_trajectory = readRDS(paste0(work_path, "/bulk_sample/sub_trajectory.rds"))

gene_overlap = intersect(rownames(moca_sub_trajectory), rownames(mutant_sub_trajectory))
moca_sub_trajectory_sub = moca_sub_trajectory[gene_overlap,]
mutant_sub_trajectory_sub = mutant_sub_trajectory[gene_overlap,]

conn <- correlation_analysis_bidirection(as.matrix(mutant_sub_trajectory_sub), as.matrix(moca_sub_trajectory_sub), fold.change = 1.5, top_gene_num = 3000, spec_gene_num = 3000)
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

heatmap.2(as.matrix(dat), 
          col = viridis,
          scale="none", 
          Rowv = TRUE, 
          Colv = TRUE, 
          keysize = 1.5,
          lhei=c(1,3.5), 
          lwid=c(1,3.5),
          density.info = "none", 
          trace = "none",
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(16, 16))



######################
### separate different groups of sub-trajectories
library(Matrix)
moca_sub_trajectory = colnames(readRDS("~/work/mutant/bulk_sample/moca/sub_trajectory.rds"))

main_trajectory_list = c("blood",
                        "endothelial",
                        "epithelial",
                        "len",
                        "liver",
                        "mesenchyme",
                        "neural_crest_glia",
                        "neural_crest_melanocyte",
                        "neural_crest_neuron",
                        "neural_tube")

moca_list = NULL
for(i in 1:length(main_trajectory_list)){
    main_trajectory_i = main_trajectory_list[i]
    print(main_trajectory_i)
    tmp = readRDS(paste0("/net/shendure/vol10/projects/cxqiu/nobackup/data/sci3_trajectory/main_trajectory_summary_monocle3/", main_trajectory_i, "_cds.pd.rds"))
    tmp_list = as.vector(unique(tmp$subtrajectory_name))
    tmp_list = tmp_list[tmp_list %in% moca_sub_trajectory]
    moca_list = rbind(moca_list, data.frame(main_trajectory = rep(main_trajectory_i, length(tmp_list)),
                                            sub_trajectory = tmp_list))
}
moca_sub_trajectory[!moca_sub_trajectory %in% as.vector(moca_list$sub_trajectory)]
moca_list = rbind(moca_list, data.frame(main_trajectory = c("neural_tube", "neural_tube", "epithelial", "epithelial"),
                                        sub_trajectory = c("Granule neuron trajectory","Wnt8b positive cell trajectory","Squamous epithelium trajectory", "Foregut epithelial trajectory")))                
moca_list$main_trajectory = as.vector(moca_list$main_trajectory)
moca_list$sub_trajectory = as.vector(moca_list$sub_trajectory)
moca_list$main_trajectory[22:24] = "other"
moca_list$main_trajectory[37:38] = "other"

mutant_sub_trajectory = colnames(readRDS("~/work/mutant/bulk_sample/sub_trajectory.rds"))
pd = readRDS("~/work/mutant/data/final/pd.rds")
mutant = unique(pd[,c("main_trajectory", "sub_trajectory")])

mutant_list = NULL
i = 1
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Haematopoiesis trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 2
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Endothelial trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 3
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Epithelial trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 5
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Hepatocyte trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 6
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory %in% c("Mesenchymal trajectory","Myoblast trajectory","Myotube trajectory","Cardiomyocyte trajectory")])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 7
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Neural crest (PNS glia) trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 8
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Melanocyte trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 9
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory %in% c("Neural crest (PNS neuron) trajectory", "Olfactory sensory neuron trajectory")])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))

i = 10
main_trajectory_i = main_trajectory_list[i]
tmp_list = as.vector(mutant$sub_trajectory[mutant$main_trajectory == "Neural tube and notochord trajectory"])
mutant_list = rbind(mutant_list, data.frame(
    main_trajectory = rep(main_trajectory_i, length(tmp_list)),
    sub_trajectory = tmp_list
))
sum(!mutant_list$sub_trajectory %in% mutant_sub_trajectory)
mutant_list$main_trajectory = as.vector(mutant_list$main_trajectory)
mutant_list$sub_trajectory = as.vector(mutant_list$sub_trajectory)
mutant_list$main_trajectory[27:28] = "other"
mutant_list$main_trajectory[41] = "other"

saveRDS(list(moca_list, mutant_list), "~/work/mutant/bulk_sample/sub_trajectory_list.rds")


################
### locally

dat = readRDS("~/work/mutant/bulk_sample/sub_trajectory_list.rds")
moca_list = dat[[1]]
mutant_list = dat[[2]]

main_trajectory_list = c("blood",
                         "endothelial",
                         "epithelial",
                         "other",
                         "mesenchyme",
                         "neural_crest_glia",
                         "neural_crest_neuron",
                         "neural_tube")

moca_sub_trajectory = readRDS(paste0(work_path, "/bulk_sample/moca/sub_trajectory.rds"))
mutant_sub_trajectory = readRDS(paste0(work_path, "/bulk_sample/sub_trajectory.rds"))
gene_overlap = intersect(rownames(moca_sub_trajectory), rownames(mutant_sub_trajectory))
moca_sub_trajectory = moca_sub_trajectory[gene_overlap,]
mutant_sub_trajectory = mutant_sub_trajectory[gene_overlap,]


for(i in c(5,8)){
    
print(i)
    
mutant_sub_trajectory_sub = mutant_sub_trajectory[,colnames(mutant_sub_trajectory) %in% as.vector(mutant_list$sub_trajectory[mutant_list$main_trajectory == main_trajectory_list[i]])]
moca_sub_trajectory_sub = moca_sub_trajectory[,colnames(moca_sub_trajectory) %in% as.vector(moca_list$sub_trajectory[moca_list$main_trajectory == main_trajectory_list[i]])]

if(main_trajectory_list[i] %in% c("neural_tube", "mesenchyme")){
    kk = 1500
} else {
    kk = 1500
}

conn <- correlation_analysis_bidirection(as.matrix(mutant_sub_trajectory_sub), as.matrix(moca_sub_trajectory_sub), fold.change = 1.5, top_gene_num = kk, spec_gene_num = kk)
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

pdf(paste0("~/work/mutant/bulk_sample/nnls_moca/", main_trajectory_list[i],".pdf"), 8, 8)
print(heatmap.2(as.matrix(dat), 
          col = viridis,
          scale="row", 
          Rowv = TRUE, 
          Colv = TRUE, 
          keysize = 1.5,
          lhei=c(1,3.5), 
          lwid=c(1,3.5),
          density.info = "none", 
          trace = "none",
          cexRow = 1, 
          cexCol = 1,
          margins = c(16, 16)))
dev.off()

}



