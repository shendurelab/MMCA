
########################################################################################################
### In this demo, we sought to repeat the cell-composition analysis which we present in the Fig. 2a. ###
########################################################################################################

### Last update: May-1, 2023

### It includes two part of analyses:
### 1) beta-binomial regression
### 2) log2-fold change of cell proportions

### The data used are provided from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/mmca/pd.rds

library(dplyr)
library(VGAM)
library(monocle3)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(scales)
library(corrplot)

########################################################
### Section - 1: performing beta-binomial regression ###
########################################################

pd = readRDS("pd.rds")
pd$genotype = pd$Mutant
pd$sample = paste0("sample_", pd$RT_group)
dat = pd %>% group_by(Background, sample, genotype, sub_trajectory) %>% summarise(celltype_n = n())

celltype_list = as.vector(unique(dat$sub_trajectory))
mutant_list = as.vector(unique(dat$genotype))

### size factor normalization using Monocle3
df = dcast(dat %>% select(sample, sub_trajectory, celltype_n), sub_trajectory ~ sample)
rownames(df) = as.vector(df[[1]])
df = df[,-1]
df[is.na(df)] = 0

cds = new_cell_data_set(as.matrix(df))
cds = monocle3::estimate_size_factors(cds)
cds_pd = pData(cds) %>%
    as.data.frame() %>%
    mutate(sample = rownames(pData(cds))) %>% 
    select(sample, Size_Factor)

dat = dat %>%
    left_join(cds_pd, by = "sample") %>%
    mutate(celltype_n_norm = round(celltype_n/Size_Factor))

dat_sub = dat %>%
    group_by(sample) %>%
    summarise(total_n_norm = sum(celltype_n_norm))

dat = dat %>%
    left_join(dat_sub, by = "sample") %>%
    mutate(cell_frac = celltype_n_norm/total_n_norm)

dat$genotype = factor(dat$genotype, levels = mutant_list)

### should we exclude cell type in which the mean cell num < 10?
celltype_list = dat %>%
    group_by(sub_trajectory) %>%
    summarise(celltype_n_norm_mean = mean(celltype_n_norm)) %>%
    filter(celltype_n_norm_mean >= 10) %>%
    pull(sub_trajectory)

dat = dat[dat$sub_trajectory %in% celltype_list,]
mutant_list = as.vector(unique(dat$genotype))

### the beta-binomial analysis was performed on samples within individual background
background_list = c("C57BL/6", "FVB", "G4")

for(xx in 1:3){
    
    res = NULL
    res_error = NULL
    
    background_i = background_list[xx]
    
    for(i in 1:length(celltype_list)){
        
        celltype_i = celltype_list[i]
        print(paste0(i, " / ", celltype_i))
        
        dat_sub = dat %>% filter(sub_trajectory == celltype_i, Background == background_i)
        count_df = cbind(dat_sub$celltype_n_norm, dat_sub$total_n_norm - dat_sub$celltype_n_norm)
        
        fit = tryCatch(
            {vglm(count_df ~ genotype, data = dat_sub, family = betabinomial, trace = FALSE)},
            error = function(e) {
                return(NA)
            }
        )
        
        if(!is.na(fit)){
            tmp = data.frame(coef(summary(fit)))
            tmp = tmp[,c(1,4)]
            names(tmp) = c("estimate", "pval")
            rownames(tmp) = gsub('genotype', '', rownames(tmp))
            tmp = tmp[rownames(tmp) %in% mutant_list,]
            tmp$mutant = rownames(tmp)
            tmp$celltype = celltype_i
            rownames(tmp) = NULL
            res = rbind(res, tmp)
        } else {
            res_error = c(res_error, celltype_i)
        }
        
    }
    
    if(background_i == "C57BL/6"){
        name = "C57BL"
    } else {
        name = background_i
    }
    
    saveRDS(list(res, res_error), paste0(name, "_pd_vgam_res.rds"))
    
    res_sub = res %>%
        mutate(fdr = p.adjust(pval, method = 'fdr')) %>%
        filter(fdr < 0.05)
    
    write.csv(res_sub, paste0(name, "_pd_vgam_res_fdr_0.05.csv"))
    
}


#############################################################################
### Section - 2: calculating log2-scaled fold changes of cell proportions ###
#############################################################################

pd = readRDS("pd.rds")
df = pd %>%
    group_by(Mutant_id, sub_trajectory, Background) %>%
    tally() %>%
    dplyr::rename(cell_num = n)

mutant_list = unique(df[,c("Mutant_id","Background")])
mutant_list = mutant_list[mutant_list$Mutant_id != "WT",]

df_mutant = NULL
for(i in 1:nrow(mutant_list)){
    mutant_i = mutant_list$Mutant_id[i]
    background_i = mutant_list$Background[i]
    print(paste0(mutant_i, ":", background_i))
    
    df_1 = df %>% filter(Background == background_i, Mutant_id == mutant_i)
    
    df_1$cell_frac = df_1$cell_num/sum(df_1$cell_num)
    
    df_2 = df %>% filter(Background == background_i, Mutant_id != mutant_i) %>%
        group_by(sub_trajectory) %>% summarise(cell_num_all = sum(cell_num))
    
    df_2$wt_cell_frac = df_2$cell_num_all/sum(df_2$cell_num_all)
    
    df_1 = df_1 %>% left_join(df_2, by = "sub_trajectory") %>%
        mutate(log_2_ratio = log2(cell_frac/wt_cell_frac))
    
    df_mutant = rbind(df_mutant, df_1)
    
}

saveRDS(df_mutant, "dat_log2fc.rds")


########################################################
### Section - 3: making the heatmap shown in Fig. 2a ###
########################################################

pd = readRDS("pd.rds")

mutant_order = unique(pd %>% filter(Mutant_id != "WT") %>% select(Mutant_id, Background)) %>%
    as.data.frame() %>%
    arrange(Background, Mutant_id)

dat = readRDS("dat_log2fc.rds")
FC = dat[,c("Mutant_id", "sub_trajectory", "log_2_ratio")]

dat_1 = readRDS("C57BL_pd_vgam_res.rds")[[1]]
dat_2 = readRDS("FVB_pd_vgam_res.rds")[[1]]
dat_3 = readRDS("G4_pd_vgam_res.rds")[[1]]
SIG = rbind(dat_1, dat_2, dat_3)

res = SIG
res$Mutant_id = gsub(" ", "_", res$mutant)
res$sub_trajectory = res$celltype

dat = FC %>% left_join(res %>% select(-estimate), by = c("Mutant_id", "sub_trajectory")) %>%
    filter(!is.na(pval))

dat$log_2_ratio[dat$log_2_ratio > 2] = 2
dat$log_2_ratio[dat$log_2_ratio < (-2)] = -2

dat$z_score = rescale(dat$log_2_ratio, to=c(-1, 1))
M = dcast(dat[,c("Mutant_id", "sub_trajectory", "z_score")], Mutant_id~sub_trajectory)
rownames(M) = as.vector(M$Mutant_id)
M = M[,-1]
M[is.na(M)] = 0
P = dcast(dat[,c("Mutant_id", "sub_trajectory", "pval")], Mutant_id~sub_trajectory)
rownames(P) = as.vector(P$Mutant_id)
P = P[,-1]
P[is.na(P)] = 10

celltype_order = pd %>%
    group_by(sub_trajectory) %>%
    tally() %>%
    arrange(desc(n))
celltype_order = celltype_order[celltype_order$sub_trajectory %in% colnames(M),]
celltype_order$sub_trajectory = factor(celltype_order$sub_trajectory, levels = as.vector(celltype_order$sub_trajectory))
p = ggplot(data=celltype_order, aes(x=sub_trajectory, y=log2(n))) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

M = M[as.vector(mutant_order$Mutant_id),as.vector(celltype_order$sub_trajectory)]
P = P[as.vector(mutant_order$Mutant_id),as.vector(celltype_order$sub_trajectory)]

pdf(paste0("Heatmap.pdf"), 10, 10)
corrplot(as.matrix(M), method="pie",
         p.mat = as.matrix(P), sig.level = 0.05, insig = "blank") 
dev.off()


















