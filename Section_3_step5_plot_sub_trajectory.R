#13. performing knn on each mutant type to annotate their main trajectories
#14. make the plot for each WT and then including Mutant cells
#15. make the plot for each sub trajectory

####################################
### step 13. performing knn on each mutant type
####################################

library(monocle3)
library(Matrix)
library(dplyr)
library(FNN)

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

pd = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
mutant_list = as.vector(unique(pd$Mutant_id))
mutant_list = mutant_list[mutant_list != "WT"]

WT_pd = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/WT_pd.rds"))
WT_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/WT_umap_coor.rds"))
WT_pd = WT_pd[rownames(WT_coor),]
y_train = as.vector(WT_pd$sub_trajectory)

k_num = 15

pd_all = NULL
for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(paste0(i, "/", mutant_i))

    Mutant_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/", mutant_i, "_umap_coor.rds"))
    Mutant_pd = pd[rownames(Mutant_coor),]

    neighbors <- get.knnx(WT_coor, Mutant_coor, k = k_num)$nn.index
    
    y_test = rep(NA, nrow(Mutant_pd))
    y_score = rep(NA, nrow(Mutant_pd))
    for(j in 1:nrow(neighbors)){
        tmp = table(y_train[neighbors[j,]])
        y_test[j] = names(tmp)[which.max(tmp)]
        y_score[j] = max(tmp)/k_num
    }
    
    Mutant_pd$sub_trajectory = as.vector(y_test)
    Mutant_pd$sub_trajectory_score = as.vector(y_score)
    
    saveRDS(Mutant_pd, paste0(work_path, "/data/sub_trajectory/", name, "/", mutant_i, "_pd.rds"))
    
    pd_all = rbind(pd_all, Mutant_pd)
}

pd_all = rbind(pd_all, WT_pd[,colnames(pd_all)])
combined_umap_coor = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_umap_coor.rds"))
if(sum(!rownames(combined_umap_coor) %in% rownames(pd_all)) != 0){
    FUN()
}
pd_all = pd_all[rownames(combined_umap_coor),]
saveRDS(pd_all, paste0(work_path, "/data/sub_trajectory/", name, "/combined_pd.rds"))


### combined pd ###
pd_all = NULL
for(i in 1:length(main_trajectory_list)){
    kk = main_trajectory_list[i]
    name = gsub('[(|)]', '', kk)
    name = gsub(' ', '_', name)
    print(paste0(i, "/", kk))
    
    tmp = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_pd.rds"))
    pd_all = rbind(pd_all, tmp)
}
pd_combined = readRDS(paste0(work_path, "/data/main_trajectory/combined_pd.rds"))
pd_all = pd_all[rownames(pd_combined),]
pd_combined$sub_trajectory = as.vector(pd_all$sub_trajectory)
pd_combined$sub_trajectory_score = as.vector(pd_all$sub_trajectory_score)
saveRDS(pd_combined, paste0(work_path, "/data/final/pd.rds"))

### combined aligned coor ###
umap_all = NULL
for(i in 1:length(main_trajectory_list)){
    kk = main_trajectory_list[i]
    name = gsub('[(|)]', '', kk)
    name = gsub(' ', '_', name)
    print(paste0(i, "/", kk))
    
    tmp = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_umap_coor.rds"))
    umap_all = rbind(umap_all, tmp)
}
umap_all = umap_all[rownames(pd_combined),]
saveRDS(umap_all, paste0(work_path, "/data/final/sub_trajectory_umap_coor.rds"))

aligned_all = list()
x_num = 0
for(i in 1:length(main_trajectory_list)){
    kk = main_trajectory_list[i]
    name = gsub('[(|)]', '', kk)
    name = gsub(' ', '_', name)
    print(paste0(i, "/", kk))
    
    tmp = readRDS(paste0(work_path, "/data/sub_trajectory/", name, "/combined_aligned_coor.rds"))
    print(dim(tmp))
    aligned_all[[kk]] = tmp
    x_num = x_num + nrow(tmp)
}
saveRDS(aligned_all, paste0(work_path, "/data/final/sub_trajectory_aligned_pca_coor.rds"))

####################################
### step 14. making plots on main trajectory
####################################

### creating color plate ###
color_plate_read = read.table(paste0(work_path, "/data/final/color_plate.txt"), as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

library(htmlwidgets)
library(plotly)
library(dplyr)

umap_coor = readRDS(paste0(work_path, "/data/final/main_trajectory_umap_coor.rds"))
pd = readRDS(paste0(work_path, "/data/final/pd.rds"))
print(sum(rownames(umap_coor) != rownames(pd)))
umap_coor = umap_coor[rownames(pd),]
pd$UMAP_1 = umap_coor[,1]
pd$UMAP_2 = umap_coor[,2]
pd$UMAP_3 = umap_coor[,3]

### ploting WT cell
pd_sub = pd[pd$Mutant_id == "WT",]

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~main_trajectory, colors = color_plate)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_WT_main_trajectory.html"))

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~log_umi)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_WT_log_umi.html"))

fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~doublet_score)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_WT_doublet_score.html"))


fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~Background, colors = color_plate)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_WT_Background.html"))

### adding Mutant cells, but downsampling to 200K cells
set.seed(123)
pd_downsampling = pd[sample(1:nrow(pd), 200000),]

fig = plot_ly(pd_downsampling, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~main_trajectory, colors = color_plate)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_All_downsample200k_main_trajectory.html"))


fig = plot_ly(pd_downsampling[sample(1:nrow(pd_downsampling)),], x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~Mutant_id, colors = color_plate)
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
saveWidget(fig, paste0(work_path, "/data/plot/", "global_All_downsample200k_Mutant.html"))

### plot cells of each mutant type
pd_sub = pd[pd$Mutant_id != "WT",]
mutant_list = names(table(pd_sub$Mutant_id))

for(i in 1:length(mutant_list)){
    mutant_i = mutant_list[i]
    print(mutant_i)
    pd_sub = pd[pd$Mutant_id == mutant_i,]
    fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~Mutant_id, colors = color_plate)
    fig = fig %>% layout(
        scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
    saveWidget(fig, paste0(work_path, "/data/plot/", "global_", mutant_i,".html"))
    
}


####################################
### step 15. making plots on sub trajectory
####################################



### creating color plate ###
color_plate_read = read.table(paste0(work_path, "/data/final/color_plate.txt"), as.is=T, sep="\t", comment.char="")
color_plate = as.vector(color_plate_read$V3)
names(color_plate) = as.vector(color_plate_read$V2)

library(htmlwidgets)
library(plotly)
library(dplyr)

umap_coor = readRDS(paste0(work_path, "/data/final/sub_trajectory_umap_coor.rds"))
pd = readRDS(paste0(work_path, "/data/final/pd.rds"))
print(sum(rownames(umap_coor) != rownames(pd)))
umap_coor = umap_coor[rownames(pd),]
pd$UMAP_1 = umap_coor[,1]
pd$UMAP_2 = umap_coor[,2]
pd$UMAP_3 = umap_coor[,3]

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

for(i in c(3,4,7)){
    kk = main_trajectory_list[i]
    name = gsub('[(|)]', '', kk)
    name = gsub(' ', '_', name)
    
    ### ploting WT cell
    pd_sub = pd[pd$Mutant_id == "WT" & pd$main_trajectory == kk,]
    
    fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~sub_trajectory, colors = color_plate)
    fig = fig %>% layout(
        scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
    saveWidget(fig, paste0(work_path, "/data/plot/", name, "_WT_sub_trajectory.html"))
    
    ### ploting all
    pd_sub = pd[pd$main_trajectory == kk,]
    if(nrow(pd_sub) > 200000){
        set.seed(123)
        pd_sub = pd_sub[sample(1:nrow(pd_sub), 200000),]
    }
    
    fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~sub_trajectory, colors = color_plate)
    fig = fig %>% layout(
        scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
    saveWidget(fig, paste0(work_path, "/data/plot/", name, "_All_sub_trajectory.html"))
    
    fig = plot_ly(pd_sub, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~Mutant_id, colors = color_plate)
    fig = fig %>% layout(
        scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                     zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))
    saveWidget(fig, paste0(work_path, "/data/plot/", name, "_All_Mutant.html"))
    
}




