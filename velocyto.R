library(monocle)
library(devtools)
library(DDRTree)
library(L1Graph)
library(reticulate)
library(monocle)
library(BiocGenerics)
library(SummarizedExperiment)
library(Seurat)
library(readr)
library(ggplot2)
library(plyr)
library(reshape2)
library(RColorBrewer)

#getting count matrix
expression_matrix <- Read10X(data.dir = "/Users/gleiria/Desktop/Projects/Vera_Martins")
filt_genes <- read_tsv("/Users/gleiria/Desktop/Projects/Vera_Martins/updated_genes.txt", col_names = FALSE)
#count matrix should have dimensions of (number of transcripts X number of barcodes)
dim(expression_matrix)
# Seurat Object
Seur_object <- CreateSeuratObject(counts = expression_matrix, project = "Vera's_project", 
                                  min.cells =3,
                                  min.features = 200)
Seur_object
# genes as a vector
vect_genes <- as.vector(filt_genes$X1)
# new filtered Seurat object (with genes of interest)
filt_Seur_object <- subset(Seur_object,
                           features = vect_genes)
# ************************************************
# ************************************************
# QC
grep("mt-", rownames(filt_Seur_object), # looking for mitochondrial genes
     value = TRUE)
#other way
# The [[]] operator can add columns to the object metadata
filt_Seur_object[["percent.mt"]] <- PercentageFeatureSet(filt_Seur_object, pattern = "mt-")
meta_data <- filt_Seur_object[[]] # refreshing meta data to add new colllumn 
#Normalizing
filt_Seur_object <- NormalizeData(filt_Seur_object,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
# Feature selection
# subset of genes showing high cell-to-cell variation
# returns 20000 features/genes by default (used for downstream PCA)
filt_Seur_object <- FindVariableFeatures(filt_Seur_object,
                                         selection.method = "vst",
                                         nfeatures = 64)
#scalling data
all.genes <- rownames(filt_Seur_object)
filt_Seur_object <- ScaleData(filt_Seur_object, features = all.genes)
# ************************************************
# ************************************************
# Linear dimensionality reduction (PCA on scalled data)
filt_Seur_object <- RunPCA(filt_Seur_object,
                           features = VariableFeatures(object = filt_Seur_object))
# Determining dimensionality
filt_Seur_object <- JackStraw(filt_Seur_object, num.replicate = 100)
filt_Seur_object <- ScoreJackStraw(filt_Seur_object, dims = 1:20)
filt_Seur_object <- FindNeighbors(filt_Seur_object, dims = 1:10)
filt_Seur_object <- FindClusters(filt_Seur_object,
                                 resolution = 1.247,
                                 algorithm = 3)
#UMAP/tSNE
filt_Seur_object <- RunUMAP(filt_Seur_object, dims = 1:10)
# ************************************************
# ***********************************************
# To subset and remove single cluster and keep the remaining clusters for new analysi
filt2_SeuratObj<- subset(filt_Seur_object, idents = c(11,13), invert = TRUE)
meta_data2 <- as.data.frame(filt2_SeuratObj@meta.data)


# preparing monocle object (getting info from Seurat)
#genes 
gene_annotation <- as.data.frame(rownames(filt2_SeuratObj@reductions[["pca"]]@feature.loadings), row.names = rownames(filt2_SeuratObj@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
#cells 
cell_metadata <- as.data.frame(filt2_SeuratObj@assays[["RNA"]]@counts@Dimnames[[2]], row.names = filt2_SeuratObj@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
#matrix 
New_matrix <- filt2_SeuratObj@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(filt2_SeuratObj@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
# creating annotated data frames
# creating barcode object
new_df <- data.frame(col1 = colnames(New_matrix),
                     labelDescription = colnames(New_matrix))
row.names(new_df) <- new_df$col1
new_df$col1 <- NULL
labelsObj <- c("rowNames", "columnNames")
feature_data <- new("AnnotatedDataFrame", data= new_df)
gene_data <- new("AnnotatedDataFrame", data=gene_annotation)
# monocle object
cds_from_seurat <- newCellDataSet(New_matrix,
                                  phenoData = feature_data,
                                  featureData = gene_data,
                                  expressionFamily = negbinomial.size())
cds_from_seurat
#normalization and variance estimation
cds_from_seurat <- estimateSizeFactors(cds_from_seurat) 
cds_from_seurat <- estimateDispersions(cds_from_seurat) 


#genes with extreme low extression are not very informative
head(dispersionTable(cds_from_seurat))

# get rid of genes with very low expression (only 2 out)
expressed_genes <- row.names(subset(fData(cds_from_seurat),
                                    num_cells_expressed >= 10))

expressed_genes_orign <- row.names(subset(fData(cds_from_seurat)))


Flt3_id <- row.names(subset(fData(cds_from_seurat), gene_short_name == "Flt3")) #marks early etp

Kit_id <- row.names(subset(fData(cds_from_seurat), gene_short_name == "Kit")) #marks ETP and DN2


Il2r_id <- row.names(subset(fData(cds_from_seurat), gene_short_name == "Il2ra")) #marks all cells DN2 onwards

Rag1_id <- row.names(subset(fData(cds_from_seurat), gene_short_name =="Rag1")) #marks DN3
Ptcra_id <- row.names(subset(fData(cds_from_seurat), gene_short_name =="Ptcra")) #marks DN3

cth <- newCellTypeHierarchy()

cth <- addCellType(cth, "early_ETP", classify_func = function(x) {x[Flt3_id,] >= 0.8})
cth <- addCellType(cth, "ETP+DN2", classify_func = function(x) {x[Kit_id,] >= 1},parent_cell_type_name = "early_ETP")

cth <- addCellType(cth, "DN2_onwards", classify_func = function(x) {x[Il2r_id,] >= 0.8})

cth <- addCellType(cth, "DN3_ptcra", classify_func = function(x) {x[Ptcra_id,] > 0.8}, parent_cell_type_name = "DN2_onwards" )
cth <- addCellType(cth, "DN3_Rag1", classify_func = function(x) {x[Rag1_id,] > 0.8},parent_cell_type_name = "DN2_onwards")



cds_from_seurat <- classifyCells(cds_from_seurat, cth, 0.1)


pie <- ggplot(pData(cds_from_seurat),
              aes(x= factor(1), fill= factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())


# clustering cells using marker genes
# this function classifies each cell into types passed to previous function
marker_diff <- markerDiffTable(cds_from_seurat,
                               cth,
                               residualModelFormulaStr = "~1",
                               cores = 4)
# ranking genes (we want top 10 or 20 for each cell type)
candidate_cluster_genes <- row.names(subset(marker_diff, qval < 0.01))
#marker_spec is a df with info of how good is a gene at classifying a group
# the closer to 1 the specificity score is the more restrictive the gene is to that cell type
marker_spec <- calculateMarkerSpecificity(cds_from_seurat[candidate_cluster_genes,],cth)
marker_spec # this feature is really cool to define new markers in new cell types for ex. 
semi_clustering_genes <- unique(selectTopMarkers(marker_spec, 64)$gene_id) #64 before
cds_from_seurat <- setOrderingFilter(cds_from_seurat, semi_clustering_genes)
plot_ordering_genes(cds_from_seurat) #black genes are used for clustering
plot_pc_variance_explained(cds_from_seurat, return_all = F)
cds_from_seurat <- reduceDimension(cds_from_seurat, max_components = 2, num_dim=3, #10 before
                                   norm_method = "log",
                                   reduction_method = "tSNE",
                                   verbose = TRUE)
cds_from_seurat <- clusterCells(cds_from_seurat, num_clusters = 12) #3 before
#store cluster info for later
#my_cluster_dim_5 <- pData(cds_from_seurat)$Cluster
myColors <- c("grey3", "blue3", "darkorchid2", "green1", "Grey","deeppink")



# get rid of genes with very low expression (only 2 out)
#expressed_genes <- row.names(subset(fData(cds_from_seurat),
                                    #num_cells_expressed >= 10))
# trajectories
# 1) choosing genes that will define progress
cds_from_seurat <- setOrderingFilter(cds_from_seurat, ordering_genes = expressed_genes)



# 2) reduce dimensionality of data

cds_from_seurat <- reduceDimension(cds_from_seurat, max_components = 2,
                                   method= "DDRTree")


# 3) ordering cells in pseudo-time
cds_from_seurat <- orderCells(cds_from_seurat,
                              reverse = TRUE)

cds_from_seurat <- orderCells(cds_from_seurat,
                              root_state = 4,
                              reverse = TRUE)

blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Il7r","Flt3")))

plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)



#**** il7r
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Il7r")))
ggplotObj1 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Il7r <- ggplotObj1[["data"]]

Il7r$expression > 0.3
Il7rPos <- Il7r[Il7r$expression > 0.3,]
Il7r_vectorPos <- as.vector(Il7rPos$Cell)

Il7r$expression < 0.3
Il7rNeg <- Il7r[Il7r$expression < 0.3,]
Il7r_vectorNeg <- as.vector(Il7rNeg$Cell)


#**** Il2ra
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Il2ra")))
ggplotObj2 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Il2ra <- ggplotObj2[["data"]]
Il2ra$expression < 0.5
Il2ra_x <- Il2ra[Il2ra$expression < 0.5,]
Il2ra_Neg_vector <- as.vector(Il2ra_x$Cell)

Il2ra <- ggplotObj2[["data"]]
Il2ra$expression > 0.5
Il2ra_Pos <- Il2ra[Il2ra$expression > 0.5,]
Il2ra_Pos_vector <- as.vector(Il2ra_Pos$Cell)

#**** Ptcra
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Ptcra")))
ggplotObj3 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Ptcra <- ggplotObj3[["data"]]
Ptcra$expression > 0.1
Ptcra_x <- Ptcra[Ptcra$expression > 0.1,]
Ptcra_vector <- as.vector(Ptcra_x$Cell)

#**** Rag
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Rag1")))
ggplotObj4 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Rag <- ggplotObj4[["data"]]
Rag$expression > 0.1
Rag_x <- Rag[Rag$expression > 0.1,]
Rag_vector <- as.vector(Rag_x$Cell)

#**** Flt3
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Flt3")))
ggplotObj5 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Flt3 <- ggplotObj5[["data"]]
Flt3$expression > 0.2
Flt3_pos <- Flt3[Flt3$expression > 0.2,]
Flt3_pos_vector <- as.vector(Flt3_pos$Cell)


Flt3$expression < 0.2
Flt3_neg <- Flt3[Flt3$expression < 0.2,]
Flt3_neg_vector <- as.vector(Flt3_neg$Cell)

#**** Bcl11b
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Bcl11b")))
ggplotObj6 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Bcl11b <- ggplotObj6[["data"]]
Bcl11b$expression >= 1
Bcl11b_pos <- Bcl11b[Bcl11b$expression >= 1,]
Bcl11b_pos_vector <- as.vector(Bcl11b_pos$Cell)


Bcl11b$expression < 1
Bcl11b_neg <- Bcl11b[Bcl11b$expression < 1,]
Bcl11b_neg_vector <- as.vector(Bcl11b_neg$Cell)

#*****
blast_genes <- row.names(subset(fData(cds_from_seurat),
                                gene_short_name %in% c("Gata3")))
ggplotObj6 <-plot_genes_jitter(cds_from_seurat[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)

Gata3 <- ggplotObj6[["data"]]
Gata3$expression > 1
Gata3_pos<- Gata3[Gata3$expression > 1,]
Gata3_pos_vec <- as.vector(Gata3_pos$Cell)

pData(cds_from_seurat)$new_cell <- "Others"


pData(cds_from_seurat)$new_cell[pData(cds_from_seurat)$State == 4 & rownames(pData(cds_from_seurat)) %in% Il2ra_Neg_vector & rownames(pData(cds_from_seurat)) %in% Bcl11b_neg_vector] <- "ETP"

barcode_df <- pData(cds_from_seurat)

ETP_DGE_vector <- as.vector(row.names(barcode_df)[which(barcode_df$new_cell== "ETP")])
Il7rPosDGE <- Il7r_vectorPos[Il7r_vectorPos %in% ETP_DGE_vector]
Flt3PosDGE <- Flt3_pos_vector[Flt3_pos_vector %in% ETP_DGE_vector]
Il7rNegDGE <- Il7r_vectorNeg[Il7r_vectorNeg %in% ETP_DGE_vector]
Flt3NegDGE <- Flt3_neg_vector[Flt3_neg_vector %in% ETP_DGE_vector]
Flt3_pos_Il7rPos <- Flt3PosDGE[Flt3PosDGE %in% Il7rPosDGE]
Flt3_pos_Il7rNeg <- Flt3PosDGE[Flt3PosDGE %in% Il7rNegDGE]
Flt3_neg_Il7rPos <- Flt3NegDGE[Flt3NegDGE %in% Il7rPosDGE]
Flt3_neg_Il7rNeg <- Flt3NegDGE[Flt3NegDGE %in% Il7rNegDGE]






