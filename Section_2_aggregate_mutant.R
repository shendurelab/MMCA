### aggregating embryos of mutant

library(Matrix)
library(dplyr)

pd = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/mutant/orig_data/df_cell_annotation.RDS")
xx = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/mutant/backup/remove_doublets/step2_removeByMonocle.rds")
pd_RT_group = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/work/mutant/orig_data/df_cell_annotated_20200929_RT_added.RDS")

pd = pd[pd$sample %in% rownames(xx)[!xx$exclude_subcluster],]
rownames(pd_RT_group) = as.vector(pd_RT_group$sample)
pd_RT_group = pd_RT_group[as.vector(pd$sample),]
pd$RT_group = as.vector(pd_RT_group$RT_group)
pd$Mutant = as.vector(pd_RT_group$Mutant)
pd$Background = as.vector(pd_RT_group$Background)
pd = pd[pd$RT_group != "NA",]
pd$embryo = paste0("mutant_", pd$RT_group)
rownames(pd) = as.vector(pd$sample)

saveRDS(pd, "/net/gs/vol2/home/cxqiu/work/mutant/data/final/pd_after_remove_doublets.rds")

pd_res = unique(pd[,c("embryo", "Background", "Mutant")])
rownames(pd_res) = as.vector(pd_res$embryo)
saveRDS(pd_res, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/pd_mutant_embryo.rds")

exp = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/orig_data/gene_count.RDS")
exp = exp[,colnames(exp) %in% rownames(pd)]
pd = pd[colnames(exp),]

### embryo
res = NULL
embryo_list = as.vector(unique(pd$embryo))
for(i in 1:length(embryo_list)){
    embryo_i = embryo_list[i]
    print(paste0(i, "/", embryo_i))
    
    exp_sub = Matrix::rowSums(exp[,pd$embryo == embryo_i])
    res = cbind(res, exp_sub)
    
}
colnames(res) = embryo_list
saveRDS(res, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/embryo.rds")


##############################
### embryo_main_trajectory ###
##############################

obs = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/data/final/pd.rds")
exp = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/orig_data/gene_count.RDS")
exp = exp[,colnames(exp) %in% rownames(obs)]
obs = obs[colnames(exp),]
obs$embryo = paste0("mutant_", obs$RT_group)
main_trajectory_list = as.vector(unique(obs$main_trajectory))
main_trajectory_list = main_trajectory_list[!main_trajectory_list %in% c("Cardiomyocyte trajectory",
                                                                         "Myoblast trajectory",
                                                                         "Myotube trajectory",
                                                                         "Olfactory sensory neuron trajectory",
                                                                         "Hepatocyte trajectory")]
main_trajectory_list = c(main_trajectory_list, "Len epithelial trajectory")
main_trajectory_list = c(main_trajectory_list, "Liver hepatocyte trajectory")

result = list()
embryo_list = as.vector(unique(obs$embryo))
for(j in 1:length(main_trajectory_list)){
    main_trajectory_j = main_trajectory_list[j]
    res = NULL
    for(i in 1:length(embryo_list)){
        embryo_i = embryo_list[i]
        print(paste0(j, "/", i, "/", embryo_i))
        
        if(main_trajectory_j == "Mesenchymal trajectory"){
            keep = obs$embryo == embryo_i & obs$main_trajectory %in% c("Cardiomyocyte trajectory",
                                                                       "Myoblast trajectory",
                                                                       "Myotube trajectory",
                                                                       "Mesenchymal trajectory")
        } else if (main_trajectory_j == "Neural crest (PNS neuron) trajectory"){
            keep = obs$embryo == embryo_i & obs$main_trajectory %in% c("Neural crest (PNS neuron) trajectory",
                                                                       "Olfactory sensory neuron trajectory")
        } else if (main_trajectory_j == "Len epithelial trajectory"){
            keep = obs$embryo == embryo_i & obs$sub_trajectory == "Len epithelial trajectory"
        } else if (main_trajectory_j == "Liver hepatocyte trajectory"){
            keep = obs$embryo == embryo_i & obs$sub_trajectory == "Liver hepatocyte trajectory"
        } else {
            keep = obs$embryo == embryo_i & obs$main_trajectory == main_trajectory_j
        }
        
        exp_sub = tryCatch( 
            {Matrix::rowSums(exp[,keep])}, 
            error = function(e) { 
                return(rep(0, nrow(exp))) 
            } 
        ) 
        
        res = cbind(res, exp_sub)
        
    }
    colnames(res) = embryo_list
    
    if(main_trajectory_j == "Melanocyte trajectory"){
        name = "Neural crest 3"
    } else if (main_trajectory_j == "Neural crest (PNS glia) trajectory"){
        name = "Neural crest 2"
    } else if (main_trajectory_j == "Neural crest (PNS neuron) trajectory"){
        name = "Neural crest 1"
    } else if (main_trajectory_j == "Len epithelial trajectory"){
        name = "Lens trajectory"
    } else if (main_trajectory_j == "Liver hepatocyte trajectory"){
        name = "Hepatocyte trajectory"
    } else {
        name = main_trajectory_j
    }
    
    result[[name]] = res
    
}
saveRDS(result, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/embryo_main_trajectory.rds")



### main_trajectory
library(Matrix)
pd = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/data/final/pd.rds")
exp = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/orig_data/gene_count.RDS")
exp = exp[,colnames(exp) %in% rownames(pd)]
pd = pd[colnames(exp),]

count = exp
count = t(t(count) / Matrix::colSums(count)) * 100000
count@x = log(count@x + 1)

res = NULL
main_trajectory_list = as.vector(unique(pd$main_trajectory))
for(i in 1:length(main_trajectory_list)){
    main_trajectory_i = main_trajectory_list[i]
    print(paste0(i, "/", main_trajectory_i))
    
    exp_sub = Matrix::rowSums(count[,pd$main_trajectory == main_trajectory_i])
    print(sum(pd$main_trajectory == main_trajectory_i))
    res = cbind(res, exp_sub)
    
}
colnames(res) = main_trajectory_list
saveRDS(res, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/main_trajectory.rds")


### main_trajectory
library(Matrix)
pd = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/data/final/pd.rds")
exp = readRDS("/net/gs/vol2/home/cxqiu/work/mutant/orig_data/gene_count.RDS")
exp = exp[,colnames(exp) %in% rownames(pd)]
pd = pd[colnames(exp),]

count = exp
count = t(t(count) / Matrix::colSums(count)) * 100000
count@x = log(count@x + 1)

res = NULL
main_trajectory_list = c("Neural tube and notochord trajectory",
                         "Endothelial trajectory",
                         "Haematopoiesis trajectory",
                         "Mesenchymal/Myotube/Myoblast/Cardiomyocyte trajectory",
                         "Neural crest (PNS neuron)/Olfactory sensory neuron trajectory",
                         "Epithelial trajectory",
                         "Liver hepatocyte trajectory",
                         "Melanocyte trajectory",
                         "Neural crest (PNS glia) trajectory",
                         "Len epithelial trajectory")

for(i in 1:length(main_trajectory_list)){
    main_trajectory_i = main_trajectory_list[i]
    print(paste0(i, "/", main_trajectory_i))
    
    if(main_trajectory_i == "Mesenchymal/Myotube/Myoblast/Cardiomyocyte trajectory"){
        exp_sub = Matrix::rowSums(count[,pd$main_trajectory %in% c("Mesenchymal trajectory", "Myotube trajectory", "Myoblast trajectory", "Cardiomyocyte trajectory")])
    } else if(main_trajectory_i == "Neural crest (PNS neuron)/Olfactory sensory neuron trajectory"){
        exp_sub = Matrix::rowSums(count[,pd$main_trajectory %in% c("Neural crest (PNS neuron) trajectory", "Olfactory sensory neuron trajectory")])
    } else if(main_trajectory_i == "Liver hepatocyte trajectory"){
        exp_sub = Matrix::rowSums(count[,pd$sub_trajectory == "Liver hepatocyte trajectory"])
    } else if(main_trajectory_i == "Len epithelial trajectory"){
        exp_sub = Matrix::rowSums(count[,pd$sub_trajectory == "Len epithelial trajectory"])
    } else {
        exp_sub = Matrix::rowSums(count[,pd$main_trajectory == main_trajectory_i])
    }

    res = cbind(res, exp_sub)
    
}
colnames(res) = main_trajectory_list
saveRDS(res, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/main_trajectory_subset.rds")



### sub_trajectory
res = NULL
sub_trajectory_list = as.vector(unique(pd$sub_trajectory))
for(i in 1:length(sub_trajectory_list)){
    sub_trajectory_i = sub_trajectory_list[i]
    print(paste0(i, "/", sub_trajectory_i))
    
    exp_sub = Matrix::rowSums(count[,pd$sub_trajectory == sub_trajectory_i])
    res = cbind(res, exp_sub)
    
}
colnames(res) = sub_trajectory_list
saveRDS(res, "/net/gs/vol2/home/cxqiu/work/mutant/bulk_sample/sub_trajectory.rds")



