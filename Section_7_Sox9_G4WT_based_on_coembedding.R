# Call from conda-environment MegaBundle
"Create density plot and h5ad object for Sox9_Regulatroy_KO + G4 WT based on the co-embedding of all mutants

Usage: Sox9_G4WT_based_on_coembedding.R --file_umap_dataframe=<file> --file_counts_matrix=<file>  --file_intron_matrix=<file> --file_exon_matrix=<file> --file_functions=<file>

Options:
  -h --help			Show this screen.
  --file_umap_dataframe=<file> 	tsv file containing umap coords and the metadata
  --file_counts_matrix=<file>		 Count Matrix from CX as RDS
  --file_intron_matrix=<file>   Spliced count matrix from CX as RDS
  --file_exon_matrix=<file>   Unspliced count matrix from CX as RDS
  --file_functions=<file>   The .R file containing all functions

" -> doc

# --- Load libraries
suppressMessages(library(dplyr))
suppressMessages(library(docopt))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Matrix))

# --- Read all the arguments passed
arguments <- docopt(doc, quoted_args=TRUE)

# --- Parameters: Read in the input parameters
file_umap_dataframe <- arguments$file_umap_dataframe
file_counts_matrix <- arguments$file_counts_matrix
file_intron_matrix <- arguments$file_intron_matrix
file_exon_matrix <- arguments$file_exon_matrix
file_functions <- arguments$file_functions

message("file_umap_dataframe: ", file_umap_dataframe)
message("file_counts_matrix: ", file_counts_matrix)
message("file_intron_matrix: ", file_intron_matrix)
message("file_exon_matrix: ", file_exon_matrix)
message("file_functions: ", file_functions)

# --- load all functions defined in the r-script
source(file_functions)

debug <- FALSE
if(debug){
  #debug mode
  file_umap_dataframe <- "mmca/analysis/figuresNextflow/density_facet_coembedding/umap_coords_all.tsv"
  file_counts_matrix="mmca/all_data/gene_count.rds"
  file_intron_matrix="mmca/all_data/MMCA_count_intron.rds"
  file_exon_matrix="mmca/all_data/MMCA_count_exon.rds"
}

# --- Read input data
umap_coords_all <- read.table(file_umap_dataframe, header=TRUE, sep="\t")
counts_matrix <- readRDS(file_counts_matrix)
intron_matrix <- readRDS(file_intron_matrix)
exon_matrix <- readRDS(file_exon_matrix)

# Select gene names common to all matrices
rownames(counts_matrix) <- gsub(rownames(counts_matrix),pattern="\\.\\d+",replacement="")
genes <- intersect(rownames(counts_matrix), rownames(exon_matrix)) %>% intersect(rownames(intron_matrix))
counts_matrix <- counts_matrix[genes, ]
intron_matrix <- intron_matrix[genes, ]
exon_matrix <- exon_matrix[genes, ]

# Select the mesenchymal cells from Sox9_Regulatroy_KO and G4 WT only based on the coembedding based umap
cells <- umap_coords_all %>%
  filter(Mutant_id %in% c("WT", "Wildtype", "Sox9_Regulatroy_KO", "Sox9 Regulatory INV")) %>%
  filter(Background == "G4") %>%
  filter(main_trajectory=="Mesenchymal trajectory") %>%
  rownames()

# subset the matrices only to have the selected cells
umap_coords <- umap_coords_all[cells,]
counts_matrix <- counts_matrix[,cells]
intron_matrix <- intron_matrix[,cells]
exon_matrix <- exon_matrix[,cells]

# Make density plot
my_cols <- RColorBrewer::brewer.pal(9,name="Blues")
my_cols[1] <- "#FFFFFF"

plot <- ggplot(umap_coords, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=0.01, col="blue", alpha=1) +
  stat_density_2d_filled(bins=9, alpha=0.9, show.legend=TRUE) +
  # scale_fill_brewer(aes(palette=Mutant_id)) +
  scale_fill_manual(values=my_cols) +
  theme_void() +
  facet_grid(.~Mutant_id) +
  theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"),
    text = element_text(size=20, family="sans"),
    strip.text.y.right = element_text(angle = 0, hjust=0),
    strip.text.x = element_text(vjust=1),
    aspect.ratio=1)

plot_filename <- "Sox9_WT_based_on_coembedding_density.pdf"
width=5*length(unique(umap_coords$Mutant_id))
ggsave(plot, filename=plot_filename, width=width+2, height=5)

# convert to h5ad and save
file_h5ad <- "Sox9_WT_based_on_coembedding.h5ad"
counts_to_anndata(counts_matrix=counts_matrix,
  spliced_matrix=exon_matrix,
  unspliced_matrix=intron_matrix,
  metadata=umap_coords,
  umap_coords=as.matrix(umap_coords[, c("UMAP_1", "UMAP_2")]),
  cleanup_genenames=TRUE, file_h5ad=file_h5ad)
