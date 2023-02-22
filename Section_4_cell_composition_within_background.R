
##############################################################################
### within each Background
### try VGAM and size factor normalization, which were suggested by Lauren ###
##############################################################################


library(dplyr)
library(VGAM)
library(monocle3)
library(reshape2)

pd = readRDS(paste0(work_path, "/data/final/pd.rds"))
pd$genotype = pd$Mutant
pd$sample = paste0("sample_", pd$RT_group)
dat = pd %>% group_by(Background, sample, genotype, sub_trajectory) %>% summarise(celltype_n = n())

celltype_list = as.vector(unique(dat$sub_trajectory))
mutant_list = as.vector(unique(dat$genotype))

### size factor normalization ###
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
### 10 cell types were excluded
exclude = dat %>%
    group_by(sub_trajectory) %>%
    summarise(celltype_n_norm_mean = mean(celltype_n_norm)) %>%
    filter(celltype_n_norm_mean < 10)

dat = dat[!dat$sub_trajectory %in% as.vector(exclude$sub_trajectory),]
celltype_list = as.vector(unique(dat$sub_trajectory))



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
    
    saveRDS(list(res, res_error), paste0(work_path, "/cell_composition/within_background/", background_i, "_pd_vgam_res.rds"))
    
    res_sub = res %>%
        mutate(fdr = p.adjust(pval, method = 'fdr')) %>%
        filter(fdr < 0.05)
    
    write.csv(res_sub, paste0(work_path, "/cell_composition/within_background/", background_i, "_pd_vgam_res_fdr_0.05.csv"))
    
    
    
}



