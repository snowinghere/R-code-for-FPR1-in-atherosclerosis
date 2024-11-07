# R-code-for-FPR1-in-atherosclerosis
AS1<-read.table( file = "RPE004.txt", header=T, row.names=1 )
AS2<-read.table( file = "RPE005.txt", header=T, row.names=1)
AS2<-read.table( file = "RPE006.txt", header=T, row.names=1,comment.char = "#",check.names=F )
AS11 <- CreateSeuratObject(counts = AS1, project = "Agedbra", min.cells = 3, min.features = 200)
AS21 <- CreateSeuratObject(counts = AS2, project = "Agedbra", min.cells = 3, min.features = 200)
AS21 <- CreateSeuratObject(counts = AS2, project = "Agedbra", min.cells = 3, min.features = 200)
merged<-merge(AS11,c(AS21,AS31))
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged,vars.to.regress ="percent.mt")
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
Idents(merged)="orig.ident"
Idents(merged)
merged <- IntegrateLayers(object = merged, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca", 
                                 new.reduction = "harmony",
                                 verbose = FALSE)
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

ElbowPlot(merged, ndims = 50)
merged <- FindNeighbors(merged,reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
merged <- RunUMAP(merged, dims = 1:30, reduction = "harmony")
merged <- RunTSNE(merged, dims = 1:30, reduction = "harmony")
library(clustree)
clustree(merged)

DimPlot(merged, reduction = "umap")+ggtitle("Harmony") 



