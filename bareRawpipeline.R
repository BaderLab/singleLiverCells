library(MASS)
library(scran)
library(org.Mm.eg.db)
library(Matrix)
library(scales)
library(Seurat)
library(cluster)
library(viridis)
forText_dropNorm <- forText_libSizeP <- NA
alldata<-list.files(pattern="matrix")
temp_genes1 <- read.table(paste0(unlist(strsplit(alldata[1],"_"))[1],"_genes.tsv"),sep="\t",header = FALSE)
anno <- temp_genes1
temp_cells<-NULL
eb_raw<-matrix(nrow = nrow(temp_genes1))
for (i in alldata){
  dataname<-unlist(strsplit(i,"_"))[1]
  temp_cells1 <- scan(paste0(dataname,"_barcodes.tsv"),character(),sep="\t")
  temp_cells1 <-paste0(rep(paste0(dataname,"_"),length(temp_cells1)),temp_cells1)
  temp_cells1 <- gsub("-","_",temp_cells1)
  temp_cells <- c(temp_cells,temp_cells1)
  raw1<-readMM(paste0(dataname,"_matrix.mtx"))
  eb_raw <- cbind(eb_raw,raw1)
}
eb_raw<-eb_raw[,2:ncol(eb_raw)]
colnames(eb_raw) <- temp_cells
rownames(eb_raw) <- temp_genes1[,2]

# Remove duplicated gene names (a couple genes are in under their MGI and HGNC symbols)
temp_r <- rownames(eb_raw)[which(duplicated(toupper(rownames(eb_raw))))]
temp_r <- lapply(temp_r,function(X) grep(paste0("^",X,"$"),rownames(eb_raw),ignore.case=T))
temp_r <- which(rownames(eb_raw) %in% 
                  names(sapply(temp_r,function(X) which.min(apply(eb_raw[X,],1,function(Y) sum(Y>0))))))
if (length(temp_r) > 0) { eb_raw <- eb_raw[-temp_r,] }
eb_raw <- eb_raw[rowSums(eb_raw) > 0,]
eb_raw <- eb_raw[,colSums(eb_raw) > 500]
rm(list=ls()[grepl("temp",ls())])
cS <- data.frame(libSize=colSums(eb_raw),
                 geneDetect=colSums(eb_raw>0))
temp_mito <- eb_raw[grepl("^MT-",rownames(eb_raw)),]
cS$mitoPct <- Matrix::colSums(temp_mito)/cS$libSize
cS2<-cS[order(cS$libSize,decreasing = T),]
eb_raw <- eb_raw[,cS$libSize > 1500]
cS <- cS[cS$libSize > 1500,]
tM <- cS$mitoPct < 0.5
eb_raw <- eb_raw[,tM]
cS <- cS[tM,]
fitLibSize <- fitdistr(cS$libSize,"negative binomial")
forText_libSizeP <- 1e-4
tB <- with(cS,geneDetect < 0)
tD <- cS$libSize < 0
eb_rawF1 <- eb_raw[,!tD]
doublets <- colnames(eb_raw)[tD]
eb_rawF1 <- eb_rawF1[rowSums(eb_rawF1) > 0,]
rm(list=ls()[grepl("temp",ls())])
rm(eb_raw)
tB1 <- tB[!tD]
eb_rawF <- eb_rawF1[,!tB1]
lowcomplex <- colnames(eb_rawF1)[tB1]
eb_rawF <- eb_rawF[rowSums(eb_rawF) > 0,]
rm(eb_rawF1)
pdf("FilterConditions.pdf")
plot(mitoPct~libSize,data=cS2,log="x",
     pch=21,cex=0.5,col=viridis(nrow(cS2),0.2),bg=viridis(nrow(cS2),0.1),
     xlab="Library Size",ylab="Mitochondrial Transcript Percent",main="Library Size and Mitochondrial Transcript Cutoffs")
abline(h=0.5,lwd=2,lty=2,col=alpha("red",0.5))
abline(v=1500,lwd=2,lty=2,col=alpha("blue",0.5))
points(mitoPct~libSize,data = cS2[rownames(cS2)%in%doublets,],pch="X",cex=0.5,col="red")
points(mitoPct~libSize,data = cS2[rownames(cS2)%in%lowcomplex,],pch="X",cex=0.5,col="green")
dev.off()
rm(cS2)
eb1 <- newSCESet(countData=eb_rawF)
eb1 <- calculateQCMetrics(eb1,feature_controls=list(Mt=grepl("^MT-",rownames(eb1))))
rm(list=ls()[grepl("raw",ls())])
eb1 <- eb1[rowSums(counts(eb1)) > 0,]
if (!file.exists("e15_cycScores.RData")) {
  colnames(anno) <- c("ENSEMBL","SYMBOL")
  cycScores <- cyclone(eb1,gene.names=anno$ENSEMBL[match(rownames(eb1), anno$SYMBOL)],
                       pairs=readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran")))
  cycScores$phases <- as.factor(cycScores$phases)
  cycScores$cellID <- rownames(eb1@phenoData)
  save(cycScores,file="e15_cycScores.RData")
} else {
  load("e15_cycScores.RData")
}
cycDlibSize <- tapply(pData(eb1)$total_counts,cycScores$phases,function(X) density(X))
cycDgeneDetect <- tapply(pData(eb1)$total_features,cycScores$phases,function(X) density(X))
geneStatsR <- with(fData(eb1),data.frame(DR=n_cells_exprs/ncol(eb1),
                                         MDTC=total_feature_counts/n_cells_exprs,
                                         MTC=total_feature_counts/ncol(eb1),
                                         sumTC=total_feature_counts))
geneStatsR$cellMax <- apply(counts(eb1),1,max)
lowCellNum <- 3
DRcut <- lowCellNum/ncol(eb1)
drop_lowCell <- geneStatsR$DR < DRcut
eb1F <- eb1[!drop_lowCell,]
rm(eb1)
if (!file.exists("e15_eb1Fnorm.RData")) {
  qClust <- quickCluster(eb1F,min.size=200) #default
  names(qClust) <- colnames(eb1F)
  forText_numQClust <- length(levels(qClust))
  eb1F <- computeSumFactors(eb1F,clusters=qClust,positive=T)
  forText_dropNorm <- sum(sizeFactors(eb1F) <= 0)
  eb1F <- eb1F[,!sizeFactors(eb1F) <= 0]
  eb1F <- normalize(eb1F)
  naCells <- apply(exprs(eb1F),2,function(X) any(is.na(X)))
  if (any(naCells)) {
    exprs(eb1F)[,naCells] <- min(apply(exprs(eb1F),1,function(X) min(X,na.rm = T)))
  }
  drop_lowGene <- rowSums(exprs(eb1F))==0
  eb1F <- eb1F[!drop_lowGene,]
  save(forText_dropNorm,forText_numQClust,qClust,eb1F,file="e15_eb1Fnorm.RData")
} else {
  load("e15_eb1Fnorm.RData")
}
geneStatsN <- data.frame(DR=apply(exprs(eb1F),1,function(X) sum(X > 0))/ncol(eb1F),
                         MDTC=apply(exprs(eb1F),1,function(X) mean(X[X > 0])),
                         MTC=rowMeans(exprs(eb1F)),sumTC=rowSums(exprs(eb1F)),
                         cellMax=apply(exprs(eb1F),1,max))
var.fit <- trendVar(eb1F,trend="semiloess",span=10,use.spikes=F)
var.out <- decomposeVar(eb1F,var.fit)
bioCut <- 0
bioCutFDR <- 1e-3
hvg.out <- var.out[which(var.out$FDR <= bioCutFDR & var.out$bio >= bioCut),]
hvg.out <- hvg.out[order(hvg.out$bio,decreasing=T),]
rm(list=ls()[!(ls() %in% c("cycScores","eb1F","hvg.out") | grepl("forText",ls()))])
if (!file.exists("e15_eb1Sjack.RData")) {
  eb1S <- new("seurat",raw.data=exprs(eb1F))
  eb1S <- Setup(eb1S,"e15",min.cells=0,min.genes=0,do.logNormalize=F,save.raw=T)
  eb1S@var.genes <- rownames(hvg.out) 
  eb1S <- PCA(eb1S,pc.genes=eb1S@var.genes,do.print=F)
  eb1S <- JackStraw(eb1S)
  save(eb1S,file="e15_eb1Sjack.RData")
} else {
  load("e15_eb1Sjack.RData")
}

if (!file.exists("e15_eb1S.RData")) {
  temp <- as.numeric(sub("PC[0-9]+ ","",levels(JackStrawPlot(eb1S,PCs=1:30)[["data"]]$PC.Score)))
  maxPCt <- min(which(temp > 1e-4)) - 1
  if(maxPCt > 30){maxPCt<-29}
  temp_PCs <- seq(1,maxPCt) #based on JackStraw output
  temp2<- RunTSNE(eb1S,dims.use=temp_PCs,do.fast=T,dim_embed = 3)
  tsne3d<-temp2@tsne.rot
  rm(temp2)
  eb1S <- RunTSNE(eb1S,dims.use=temp_PCs,do.fast=T)
  minDEgenes <- uniqDE <- list()
  k <- 0; minDEgenesToNeighbour <- 100
  rn <- 0
  while (minDEgenesToNeighbour > 1 & rn < 11) {
    rn <- rn+1
    if (minDEgenesToNeighbour <= 30) {
      k <- k + 0.2
    } else {
      k <- k + 0.4
    }
    print(paste0("~~~~~~~~~~~~ res.",k," ~~~~~~~~~~~~"))
    if (!any(grepl("res",colnames(eb1S@data.info)))) {
      eb1S <- FindClusters(eb1S,pc.use=temp_PCs,print.output=F,save.SNN=T,resolution=k)
    } else {
      eb1S <- FindClusters(eb1S,pc.use=temp_PCs,print.output=F,reuse.SNN=T,resolution=k)
    }
    res <- colnames(eb1S@data.info)[length(colnames(eb1S@data.info))]
    eb1S@ident <- eb1S@data.info[,res] <- as.factor(as.integer(eb1S@data.info[,res]) + 1)
    names(eb1S@ident) <- rownames(eb1S@data.info)
    CL <- makeCluster(detectCores() - 4)
    tempCL <- clusterEvalQ(CL,library(Seurat))
    clusterExport(CL,"eb1S")
    
    deGall <- parLapplyLB(cl=CL,levels(eb1S@ident),function(X) {
      temp <- FindMarkers(eb1S,X)
      temp$fwer <- p.adjust(temp$p_val,"holm",nrow(eb1S@data))
      return(temp[temp$fwer <= 1e-2 & temp$avg_diff > 0,]) # postive DE at FWER of 1%
    })
    
    deGvs <- parApply(CL,combn(levels(eb1S@ident),2),2,function(cl) {
      tempOut <- FindMarkers(eb1S,ident.1=cl[1],ident.2=cl[2])
      tempOut$fwer <- p.adjust(tempOut$p_val,"holm",n=nrow(eb1S@data))
      tempOut$clustPair <- factor(paste(cl,collapse="~"))
      tempOut$gene <- rownames(tempOut)
      return(tempOut[tempOut$fwer <= 5e-2,]) # FWER at 5%
    })
    stopCluster(CL)
    names(deGvs) <- apply(combn(levels(eb1S@ident),2),2,function(X) paste(X,collapse="~"))
    deGvs <- do.call(rbind,deGvs[order(names(deGvs))])
    deGvs$posClust <- as.factor(mapply(function(a,b) strsplit(as.character(a),"~")[[1]][(b < 0)+1],
                                       deGvs$clustPair,deGvs$avg_diff))
    uniqDE[[res]] <- tapply(deGvs$gene,deGvs$posClust,function(X) 
      names(table(X))[table(X) == length(levels(deGvs$posClust))-1])
    
    tempDist <- dist(eb1S@pca.rot[,seq(1,maxPCt)],method="euclidean")
    tempSil <- silhouette(as.integer(eb1S@ident),dist=tempDist)
    pw <- list()
    cc <- t(apply(tempSil[,1:2],1,sort))
    for (x in sort(unique(cc[,1]))) {
      for (y in sort(unique(cc[cc[,1] == x,2]))) {
        pw[[length(pw)+1]] <- c(x,y)
      }
    }
    pw <- sapply(pw,function(X) paste(X,collapse="~"))
    minDEgenesToNeighbour <- with(deGvs,tapply(fwer,clustPair,function(X) sum(X <= 1e-2)))
    minDEgenesToNeighbour <- min(minDEgenesToNeighbour[pw])
    if (is.na(minDEgenesToNeighbour)){minDEgenesToNeighbour<-0}
    minDEgenes[[res]] <- minDEgenesToNeighbour
    
    DR <- apply(eb1S@data,1,function(X) tapply(X,eb1S@ident,function(Y) sum(Y>0)/length(Y)))
    MDTC <- apply(eb1S@data,1,
                  function(X) tapply(X,eb1S@ident,
                                     function(Y) {
                                       temp <- mean(Y[Y>0])
                                       if (is.na(temp)) { temp <- 0 }
                                       return(temp)
                                     }))
    MTC <- apply(eb1S@data,1,function(X) tapply(X,eb1S@ident,mean))
    CGS <- lapply(levels(eb1S@ident), function(X) data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]))
    
    save(deGvs,file=paste0("e15_precalc_",gsub(".","",res,fixed=T),"_deGvs.RData"))
    save(deGall,file=paste0("e15_precalc_",gsub(".","",res,fixed=T),"_deGall.RData"))
    save(CGS,file=paste0("e15_precalc_",gsub(".","",res,fixed=T),"_CGS.RData"))
  }
  pDat <- pData(eb1F)[,c("total_counts","total_features")]
  save(eb1S,maxPCt,minDEgenes,uniqDE,pDat,tsne3d,file="e15_eb1S.RData")
  rm(list=ls()[!(ls() %in% c("cycScores","eb1F","hvg.out","pwGenes",
                             "tempConds","i","CL") | grepl("forText",ls()))])
  gc()
} else {
  load("e15_eb1S.RData")
}
