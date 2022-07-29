library(BiocManager)
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(Matrix)
library(ggplot2)
library(loomR)
library(hdf5r)
library(R6)
library(ggpubr)
library(grid)
library(gridExtra)
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


GetAssayData(object = Seur_object, slot = "data")

 
# genes as a vector
vect_genes <- as.vector(filt_genes$X1)

# new filtered Seurat object (with genes of interest)
filt_Seur_object <- subset(Seur_object,
                           features = vect_genes)

# ************************************************
# ************************************************

grep("mt-", rownames(filt_Seur_object), # looking for mitochondrial genes
     value = TRUE)

grep("Ccn", rownames(Seur_object), # looking for mitochondrial genes
     value = TRUE)
#Ccnj 
#Ccna2

#other way
# The [[]] operator can add columns to the object metadata

filt_Seur_object[["percent.mt"]] <- PercentageFeatureSet(filt_Seur_object, pattern = "mt-")
meta_data <- filt_Seur_object[[]] # refreshing meta data to add new colllumn 

# Visualize QC metrics as violins
VlnPlot(filt_Seur_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)

#Normalizing
filt_Seur_object <- NormalizeData(filt_Seur_object,
                                  normalization.method = "LogNormalize",
                                  scale.factor = 10000)
# checking normalized values
XXX <- as.data.frame(filt_Seur_object[["RNA"]]@data)
XXX <- t(XXX)
XXX <- as.data.frame(XXX)
hist(XXX$Il7r)
# Feature selection
# subset of genes showing high cell-to-cell variation
# returns 20000 features/genes by default (used for downstream PCA)
filt_Seur_object <- FindVariableFeatures(filt_Seur_object,
                                         selection.method = "vst",
                                         nfeatures = 64)

# genes showing highest cell to cell variation 
top10 <- head(VariableFeatures(filt_Seur_object),10)
top10
top20 <- head(VariableFeatures(filt_Seur_object),20)
top20
top50 <- head(VariableFeatures(filt_Seur_object),40)
top50

#scalling data
all.genes <- rownames(filt_Seur_object)
filt_Seur_object <- ScaleData(filt_Seur_object, features = all.genes)

# ************************************************
# ************************************************

# Linear dimensionality reduction (PCA on scalled data)
filt_Seur_object <- RunPCA(filt_Seur_object,
                           features = VariableFeatures(object = filt_Seur_object))

# Visualizing
DimPlot(filt_Seur_object, reduction = "pca")
VizDimLoadings(filt_Seur_object, dims = 1:2, reduction = "pca")

# Determining dimensionality
filt_Seur_object <- JackStraw(filt_Seur_object, num.replicate = 100)
filt_Seur_object <- ScoreJackStraw(filt_Seur_object, dims = 1:20)

# Variance that each PC accounts for
ElbowPlot(filt_Seur_object)


filt_Seur_object <- FindNeighbors(filt_Seur_object, dims = 1:10)

filt_Seur_object <- FindClusters(filt_Seur_object,
                                 resolution = 0.5,
                                 algorithm = 3)

# cluster IDs of the first 5 cells
head(Idents(filt_Seur_object),5)

#UMAP/tSNE
filt_Seur_object <- RunUMAP(filt_Seur_object, dims = 1:10)

# label = TRUE can be setted or use LabelClusters 
# individual clusters
plot_dim <- DimPlot(filt_Seur_object, reduction = "umap",
                    pt.size = 1.2,
                    label = FALSE)
plot_dim


# How many cells are in each cluster
table(Idents(filt_Seur_object))
#There are 1333 cells in cluster 1
# What are the cell names of all NK cells?
cluster0_vector <- WhichCells(filt_Seur_object, idents = 0)
cluster1_vector <- WhichCells(filt_Seur_object, idents = 1)
cluster2_vector <- WhichCells(filt_Seur_object, idents = 2)
cluster3_vector <- WhichCells(filt_Seur_object, idents = 3)
cluster4_vector <- WhichCells(filt_Seur_object, idents = 4)
cluster5_vector <- WhichCells(filt_Seur_object, idents = 5)
cluster6_vector <- WhichCells(filt_Seur_object, idents = 6)
cluster7_vector <- WhichCells(filt_Seur_object, idents = 7)
#ggsave("umap_1.pdf", plot_dim)
ggsave("umap_7_plus_1.pdf", plot_dim, width = 3, height = 3, units = "in", scale = 3, path = "/Users/gleiria/Desktop/Plots_vera")

`%!in%` <- Negate(`%in%`)

title1=text_grob("Expression of Flt3 and Il7r in State 4", size = 25, face = "bold")   #### this worked for me
grid.arrange(plot_dim3, plot_dim2, plot_dim5, plot_dim4 ,nrow = 2,
             top = title1)
             






plot_dim2 <- DimPlot(filt_Seur_object, reduction = "umap",
                     cells.highlight = overlap_barcodes,
                     cols.highlight = "blue",cols = "grey", order = TRUE,
                     sizes.highlight = 0.5) +
  theme(legend.position="none") +
  ggtitle("Flt3(+) Il7r(+)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=15))
  
                  
plot(plot_dim2)

plot_dim3 <- DimPlot(filt_Seur_object, reduction = "umap",
                     cells.highlight = Flt3Yes_Il7rNo_vector,
                     cols.highlight = "blue",cols = "grey", order = TRUE,
                     sizes.highlight = 0.5) +
  theme(legend.position="none") +
  ggtitle("Flt3(+) Il7r(-)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=15))

plot(plot_dim3)

plot_dim4 <- DimPlot(filt_Seur_object, reduction = "umap",
                     cells.highlight = Il7rYes_Flt3Novector,
                     cols.highlight = "blue",cols = "grey", order = TRUE,
                     sizes.highlight = 0.5) +
  theme(legend.position="none") +
  ggtitle("Flt3(-) Il7r(+)") +
  theme(plot.title = element_text(hjust = 0.3)) +
  theme(plot.title = element_text(size=15))
plot_dim4

plot_dim5 <- DimPlot(filt_Seur_object, reduction = "umap",
                     cells.highlight = none_vector,
                     cols.highlight = "blue",cols = "grey", order = TRUE,
                     sizes.highlight = 0.5) +
  theme(legend.position="none") +
  ggtitle("Flt3(-) Il7r(-)") +
  theme(plot.title = element_text(hjust = 0.3)) +
  theme(plot.title = element_text(size=15))
plot_dim5



ggp <- ggplot(NULL) +
  plot_dim2 + 
  plot_dim3 +
  plot_dim4 +
  plot_dim5
ggp



# Feature plot
plot_x <- FeaturePlot(filt_Seur_object,
                      features = c("Flt3", "Il7r"),
                      pt.size = 0.2,
                      cells = overlap_barcodes)
plot_x

filt_Seur_object %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top15
DoHeatmap(filt_Seur_object, features = top20) + NoLegend()

DimPlot(filt_Seur_object,
        reduction = "umap",
        label = TRUE,
        pt.size = 0.1) + NoLegend()
# ************************************************
# ************************************************

filt_Seur_object %>%
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) -> top15
DoHeatmap(filt_Seur_object, features = top50) + NoLegend()

#===========
new.cluster.ids <- c("ETP", "DN2-DN3", "NotSure", "NotSure", "NotSure", "ETP",
                     "DN3", "NotSure", "DN3","NotSure","OutGroup","NotSure","OutGroup")
names(new.cluster.ids) <- levels(filt_Seur_object)
filt_Seur_object <- RenameIdents(filt_Seur_object, new.cluster.ids)
DimPlot(filt_Seur_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# ************************************************
# ************************************************

# find all markers of cluster 0
cluster0.markers <- FindMarkers(filt_Seur_object, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n=10)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(filt_Seur_object, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n=10)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(filt_Seur_object, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n=10)

# find all markers of cluster 3
cluster3.markers <- FindMarkers(filt_Seur_object, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n=10)

# find all markers of cluster 4
cluster4.markers <- FindMarkers(filt_Seur_object, ident.1 = 4, min.pct = 0.25)
head(cluster6.markers, n=10)

# find all markers of cluster 5
cluster5.markers <- FindMarkers(filt_Seur_object, ident.1 = 5, min.pct = 0.25)
head(cluster3.markers, n=10)

# find all markers of cluster 6
cluster6.markers <- FindMarkers(filt_Seur_object, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n=10)

# find all markers of cluster 7
cluster7.markers <- FindMarkers(filt_Seur_object, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n=10)














# fins all markers distinguishing cluster 1 from cluster 4 and 2
clusterXXX.markers <- FindMarkers(filt_Seur_object,
                                  ident.1 = 1,
                                  ident.2 = c(2,4),
                                  min.pct = 0.25)
head(clusterXXX.markers, n = 5)

# expression probability
VlnPlot(filt_Seur_object,
        features = c("Il7r", "Bcl11b"),
        slot = "counts",
        log = TRUE)

VlnPlot(filt_Seur_object, 
        features = c("Il7r", "Bcl11b"))

##***********************************************************************************************************########
##***********************************************************************************************************#######


# To subset and remove single cluster and keep the remaining clusters for new analysi

filt2_SeuratObj<- subset(filt_Seur_object, idents = c(11,13), invert = TRUE)
meta_data2 <- as.data.frame(filt2_SeuratObj@meta.data)



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








plot_x <- FeaturePlot(filt2_SeuratObj,
                      features = c("Flt3", "Il7r","Bcl11b","Il2ra"),
                      pt.size = 0.05,
                      cols = c("grey", "green"),
                      combine = TRUE) & NoLegend() &  NoAxes() 



plot_x 


  

#df_1 <- filt_Seur_object@meta.data

#df_2 <- data.frame(col1 = row.names(df_1),
                   #col2 = df_1$seurat_clusters)

#df_2$col2 == 11 & df_2$col2 == 13 

#df_3 <- df_2[df_2$col2 == 13,]
#df_4 <- df_2[df_2$col2 == 11,]

#cell_vec1 <- as.vector(df_3$col1)
#cell_vec2 <- as.vector(df_4$col1)
#cell_vecX <- c(cell_vec2,cell_vec1)

#cell_vec_keep <- as.vector(row.names(df_1))

#final_vec <- setdiff(cell_vec_keep, cell_vecX)
#write.csv(final_vec, file = "final_vector.csv")




#library(sceasy)
#library(reticulate)
#library(anndata)


#sceasy::convertFormat(filt2_SeuratObj, from="seurat", to="anndata",
                      #outFile='filt2_SeuratObj.h5ad')



barcode_df <- pData(cds_from_seurat)






