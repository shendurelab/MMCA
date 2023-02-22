### check cell composition changes for each mutant type ###


#######################################################################
### 1 - first plot
### for each major/sub-trajectories, we make the log(mutant% / wt%) for each mutant type 
### (wt with the same background are used)
### generally, seperated by major/sub-trajectories for each PDF
#######################################################################

library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/final/pd.rds"))

### main trajectory
df = pd %>%
    group_by(Mutant_id, main_trajectory, Background) %>%
    tally() %>%
    dplyr::rename(cell_num = n)

df_sub1 = pd %>%
    group_by(Mutant_id, Background) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)

df = df %>%
    left_join(df_sub1, by = c("Mutant_id", "Background")) %>%
    mutate(cell_frac = cell_num/cell_num_total)

df_wt = df %>%
    filter(Mutant_id == "WT") %>%
    select(Mutant_id, main_trajectory, Background, wt_cell_frac = cell_frac)

df_mutant = df %>%
    filter(Mutant_id != "WT") %>%
    left_join(df_wt[,c("main_trajectory", "wt_cell_frac", "Background")], by = c("main_trajectory", "Background")) %>%
    mutate(log_2_ratio = log2(cell_frac/wt_cell_frac)) 

### creating color plate ###
color_plate_read = read.table(paste0(work_path, "/code/color_plate.txt"), as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

main_trajectory_list = as.vector(unique(df_mutant$main_trajectory))
for(i in 1:length(main_trajectory_list)){
    
    xx = 'C57BL/6'
    df_sub = df_mutant %>%
        filter(main_trajectory == main_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x=main_trajectory_list[i], y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p1 = p + scale_fill_manual(values = color_plate)
    
    
    xx = 'G4'
    df_sub = df_mutant %>%
        filter(main_trajectory == main_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x='', y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p2 = p + scale_fill_manual(values = color_plate)
    
    
    xx = 'FVB'
    df_sub = df_mutant %>%
        filter(main_trajectory == main_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x='', y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p3 = p + scale_fill_manual(values = color_plate)
    
    name = main_trajectory_list[i]
    name = gsub('[/]','_',name)
    name = gsub('[(|)]','',name)
    name = gsub('[-]','_',name)
    name = gsub(' ','_',name)
    pdf(paste0(work_path, '/cell_composition/main_trajectory/', name, '.pdf'), 10, 5)
    grid.arrange(p1, p2, p3, nrow=1, ncol=3)
    dev.off()
    
}



### sub trajectory
df = pd %>%
    group_by(Mutant_id, sub_trajectory, Background) %>%
    tally() %>%
    dplyr::rename(cell_num = n)

df_sub1 = pd %>%
    group_by(Mutant_id, Background) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)

df = df %>%
    left_join(df_sub1, by = c("Mutant_id", "Background")) %>%
    mutate(cell_frac = cell_num/cell_num_total)

df_wt = df %>%
    filter(Mutant_id == "WT") %>%
    select(Mutant_id, sub_trajectory, Background, wt_cell_frac = cell_frac)

df_mutant = df %>%
    filter(Mutant_id != "WT") %>%
    left_join(df_wt[,c("sub_trajectory", "wt_cell_frac", "Background")], by = c("sub_trajectory", "Background")) %>%
    mutate(log_2_ratio = log2(cell_frac/wt_cell_frac)) 

### creating color plate ###
color_plate_read = read.table(paste0(work_path, "/code/color_plate.txt"), as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

sub_trajectory_list = as.vector(unique(df_mutant$sub_trajectory))
for(i in 1:length(sub_trajectory_list)){
    
    xx = 'C57BL/6'
    df_sub = df_mutant %>%
        filter(sub_trajectory == sub_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x=sub_trajectory_list[i], y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p1 = p + scale_fill_manual(values = color_plate)
    
    
    xx = 'G4'
    df_sub = df_mutant %>%
        filter(sub_trajectory == sub_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x='', y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p2 = p + scale_fill_manual(values = color_plate)
    
    
    xx = 'FVB'
    df_sub = df_mutant %>%
        filter(sub_trajectory == sub_trajectory_list[i], Background == xx)
    
    p = ggplot(df_sub, aes(Mutant_id, log_2_ratio, fill = Mutant_id))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x='', y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 9) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p3 = p + scale_fill_manual(values = color_plate)
    
    name = sub_trajectory_list[i]
    name = gsub('[/]','_',name)
    name = gsub('[(|)]','',name)
    name = gsub('[-]','_',name)
    name = gsub(' ','_',name)
    pdf(paste0(work_path, '/cell_composition/sub_trajectory/', name, '.pdf'), 10, 5)
    grid.arrange(p1, p2, p3, nrow=1, ncol=3)
    dev.off()
    
}









#######################################################################
### 2 - second plot
### for each mutant type, we make the log(mutant% / wt%) for each major/sub-trajectories
### (wt with the same background are used)
### generally, seperated by mutant types for each PDF
#######################################################################




work_path = "~/work/mutant"
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(dplyr)

pd = readRDS(paste0(work_path, "/data/final/pd.rds"))

### main trajectory ###############
df = pd %>%
    group_by(Mutant_id, main_trajectory, Background) %>%
    tally() %>%
    dplyr::rename(cell_num = n)

df_sub1 = pd %>%
    group_by(Mutant_id, Background) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)

df = df %>%
    left_join(df_sub1, by = c("Mutant_id", "Background")) %>%
    mutate(cell_frac = cell_num/cell_num_total)

df_wt = df %>%
    filter(Mutant_id == "WT") %>%
    select(Mutant_id, main_trajectory, Background, wt_cell_frac = cell_frac)

df_mutant_1 = df %>%
    filter(Mutant_id != "WT") %>%
    left_join(df_wt[,c("main_trajectory", "wt_cell_frac", "Background")], by = c("main_trajectory", "Background")) %>%
    mutate(log_2_ratio = log2(cell_frac/wt_cell_frac)) 

pd_num = pd %>% ### identify major/sub-trajectory which have cell number per embryo < 10
    group_by(RT_group, main_trajectory) %>%
    tally() %>%
    group_by(main_trajectory) %>%
    summarize(n_mean = mean(n)) %>%
    filter(n_mean >= 20)

df_mutant_1 = df_mutant_1 %>%
    filter(main_trajectory %in% as.vector(pd_num$main_trajectory))


### sub trajectory ##################
df = pd %>%
    group_by(Mutant_id, sub_trajectory, Background) %>%
    tally() %>%
    dplyr::rename(cell_num = n)

df_sub1 = pd %>%
    group_by(Mutant_id, Background) %>%
    tally() %>%
    dplyr::rename(cell_num_total = n)

df = df %>%
    left_join(df_sub1, by = c("Mutant_id", "Background")) %>%
    mutate(cell_frac = cell_num/cell_num_total)

df_wt = df %>%
    filter(Mutant_id == "WT") %>%
    select(Mutant_id, sub_trajectory, Background, wt_cell_frac = cell_frac)

df_mutant_2 = df %>%
    filter(Mutant_id != "WT") %>%
    left_join(df_wt[,c("sub_trajectory", "wt_cell_frac", "Background")], by = c("sub_trajectory", "Background")) %>%
    mutate(log_2_ratio = log2(cell_frac/wt_cell_frac)) 

pd_num = pd %>% ### identify major/sub-trajectory which have cell number per embryo < 10
    group_by(RT_group, sub_trajectory) %>%
    tally() %>%
    group_by(sub_trajectory) %>%
    summarize(n_mean = mean(n)) %>%
    filter(n_mean >= 20)

df_mutant_2 = df_mutant_2 %>%
    filter(sub_trajectory %in% as.vector(pd_num$sub_trajectory))

### creating color plate ###
color_plate_read = read.table("~/work/scripts/mutant/color_plate.txt", as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

mutant_list = as.vector(unique(df_mutant$Mutant_id))
for(i in 1:length(mutant_list)){
    print(paste0(i, "/", length(mutant_list)))
    
    df_sub = df_mutant_1 %>%
        filter(Mutant_id == mutant_list[i])
    
    p = ggplot(df_sub, aes(main_trajectory, log_2_ratio, fill = main_trajectory))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x=mutant_list[i], y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 10) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p1 = p + scale_fill_manual(values = color_plate)
    
    df_sub = df_mutant_2 %>%
        filter(Mutant_id == mutant_list[i])
    
    p = ggplot(df_sub, aes(sub_trajectory, log_2_ratio, fill = sub_trajectory))
    p = p + geom_bar(stat='identity', width = 0.5) + coord_flip()
    p = p + labs(x=mutant_list[i], y="log2( % of Mutant / % of WT )", title=xx)
    p = p + theme_classic(base_size = 10) + theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    p = p + theme(legend.position="none")
    p2 = p + scale_fill_manual(values = color_plate)

    name = mutant_list[i]

    pdf(paste0(work_path, '/cell_composition/each_mutant/', name, '.pdf'), 10, 10)
    grid.arrange(p1, p2, nrow=1, ncol=2)
    dev.off()
    
}





#######################################################################
### 4 - forth plot
### creating a heatmap used to show the log2(ratio)
#######################################################################

res_G4 = readRDS(paste0(work_path, "/cell_composition/within_background/G4_pd_vgam_res.rds"))
res_FVB = readRDS(paste0(work_path, "/cell_composition/within_background/FVB_pd_vgam_res.rds"))
res_C57BL = readRDS(paste0(work_path, "/cell_composition/within_background/C57BL_pd_vgam_res.rds"))

res = rbind(res_G4[[1]], res_FVB[[1]], res_C57BL[[1]])
res$Mutant_id = gsub(" ", "_", res$mutant)
res$sub_trajectory = res$celltype

dat = df_mutant %>%
    left_join(res[,c("Mutant_id", "sub_trajectory", "pval")], by = c("Mutant_id", "sub_trajectory"))

new_name = read.table("~/work/mutant/figures/new_list.txt", header=F, as.is=T, sep="\t")
names(new_name) = c("Mutant_id", "New_Mutant_id")
dat_out = dat %>% 
    left_join(new_name, by = "Mutant_id") %>% 
    select(New_Mutant_id, sub_trajectory, Background, cell_num, cell_num_total, cell_frac, wt_cell_frac, log_2_ratio, pval)
write.table(dat_out, "~/work/mutant/figures/supplementary_table.txt", row.names=F, col.names=T, quote=F, sep="\t")

dat = dat[dat$Mutant_id %in% as.vector(res$Mutant_id) &
              dat$sub_trajectory %in% as.vector(res$sub_trajectory),]

hist(dat$log_2_ratio, 100)
dat$log_2_ratio[dat$log_2_ratio > 2] = 2
dat$log_2_ratio[dat$log_2_ratio < (-2)] = -2

dat_sig = dat %>%
    filter(pval < 0.05)
saveRDS(dat_sig, paste0(work_path, "/cell_composition/within_background/sig_pval_0.05.rds"))

library(reshape2)
library(scales)
library(corrplot)
library(ggplot2)
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
ggplot(data=celltype_order, aes(x=sub_trajectory, y=n)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

mutant_order = unique(df_mutant[,c("Mutant_id", "Background")]) %>%
    as.data.frame() %>%
    arrange(Background, Mutant_id)

M = M[as.vector(mutant_order$Mutant_id),as.vector(celltype_order$sub_trajectory)]
P = P[as.vector(mutant_order$Mutant_id),as.vector(celltype_order$sub_trajectory)]

corrplot(as.matrix(M), method="pie",
         p.mat = as.matrix(P), sig.level = 0.05, insig = "blank")

write.table(rownames(M), "~/work/mutant/figures/tmp/cell_composition_rownames.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.table(colnames(M), "~/work/mutant/figures/tmp/cell_composition_colnames.txt", sep="\t", quote=F, row.names=F, col.names=F)



