# singleLiverCells
Scripts for processing single liver cells RNA-Seq 10X Genomics Data

The 20 Clusters of total Liver cells were determined by the package "Seurat" in R using the following parameters:
FindClusters(Data,pc.use=1:29,print.output=F,save.SNN=T,resolution=0.8)

The 9 Clusters of Immune Cell Subset of Liver were determined by the package "Seurat" in R using the following parameters:
FindClusters(Data,pc.use=1:9,print.output=F,save.SNN=T,resolution=0.4)
