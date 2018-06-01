#Output Brendan's Pipeline results as pdf
library(Seurat)
library(cluster)
library(gplots)
library(viridis)
library(RColorBrewer)
library(scales)
library(TeachingDemos)
library(vioplot)
library(parallel)
require(rgl)

#load results
timePoint <- "e15"
if (!exists("eb1S")) {
  load(paste0(timePoint,"_eb1S.RData"))
}

if (!exists("cycScores")) {
  load(paste0(timePoint,"_cycScores.RData"))
}

#color gradient
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#define signature genes for potential clusters
m1<-readline(prompt="Please enter the full name of the file (including .txt or .csv) listing Cell Type and Marker Genes (Cell Type 1st column, gene 2nd column): ")
markers<-read.table(m1, header = TRUE,sep = "\t")
markers2<-names(table(markers[,1]))
cellMarkers <- list()
for (i in 1:length(markers2)){
  temp1<-as.character(markers[which(markers[,1]==markers2[i]),2])
  cellMarkers[[i]]<-temp1[!duplicated(temp1)]
}
names(cellMarkers)<-markers2

if (file.exists(paste0(timePoint,"_savedRes.RData"))) {
  load(paste0(timePoint,"_savedRes.RData"))
} else {
  savedRes <- NULL
}

#Select the resolution for later analysis
cat("Please choose the resolution for clustering (#): \n")
for(i in 1:floor(length(minDEgenes)/2)){
  cat("#",2*i-1,"\t :",names(minDEgenes)[2*i-1],"\t \t","#",2*i,"\t :",names(minDEgenes)[2*i],"\n")
}
f3<-as.numeric(readline(prompt="Please enter the # of the resolution used for clustering (#): "))

#Show the resolution, number of clusters, and number of DE genes among the clusters
numClust <- apply(eb1S@data.info[,grepl("res",colnames(eb1S@data.info))],2,function(Y) length(unique(Y)))
plot(x=numClust,y=unlist(minDEgenes),type="b",cex=1.2,
     xlab="Number of clusters",
     ylab="DE genes (@FWER <= 1e-2) b/w most similar clusters")
abline(h=seq(0,max(unlist(minDEgenes)),10),lty=3,col=alpha(1,0.3))
points(x=numClust[f3],y=unlist(minDEgenes)[f3],
       pch=16,cex=1.5,col="red")
readline(prompt="Press Enter to continue: ")

#Set parameters
clusts <- eb1S@data.info[,names(minDEgenes)[f3]]
CGS <- list()
  load(paste0(timePoint,"_precalc_",gsub(".","",names(minDEgenes)[f3],fixed=T),"_CGS.RData"))
  for (i in seq_along(CGS)) {
    CGS[[i]]$MTCrank <- rank(CGS[[i]]$MTC,ties.method="min")/nrow(CGS[[i]])
    CGS[[i]]$overCut <- CGS[[i]]$MTC > mean(CGS[[i]]$MTC)
    CGS[[i]]$genes <- rownames(CGS[[i]])
  }
clusterID <- names(sapply(CGS,function(X) which.max(sapply(cellMarkers,function(Y) median(X$MTC[X$genes %in% Y])))))
if (length(levels(clusts)) <= 8) {
  clustCols<-brewer.pal(length(levels(clusts)),"Dark2")
} else {
  clustCols<-gg_colour_hue(length(levels(clusts)))
}

#silhouette plot
tempDist <- dist(eb1S@pca.rot[,seq(1,maxPCt)],method="euclidean")
tempSil <- silhouette(as.integer(clusts),dist=tempDist)
plot(tempSil,beside=T,border=NA,main=NA,col=clustCols,do.n.k=T)
readline(prompt="Press Enter to continue: ")

#tSNE Plot with Cluster#
transp <- data.frame(bg=rep(0.5,ncol(eb1S@data)),col=rep(1,ncol(eb1S@data)))
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste(paste0("Labelled by Cluster#: ", nrow(eb1S@tsne.rot)," Cells"),
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=alpha(clustCols[clusts],transp$col),
       bg=alpha(clustCols[clusts],transp$bg))
text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],INDEX=clusts),
     y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],INDEX=clusts),
     labels=levels(clusts),col="black",font=2,cex=1.5)
readline(prompt="Press Enter to continue: ")

#tSNE Plot with Potential cluster names
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste("Labelled by Potential Cell Type",
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=alpha(clustCols[clusts],transp$col),
       bg=alpha(clustCols[clusts],transp$bg))
text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],INDEX=clusts),
     y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],INDEX=clusts),
     labels=clusterID,col="black",font=2,cex=0.8)
readline(prompt="Press Enter to continue: ")

#3D plots
plot3d(tsne3d,type="s", size=1,col=alpha(clustCols[clusts],transp$col),
       bg=alpha(clustCols[clusts],transp$bg))
text3d(x=(tapply(FUN=mean,X=tsne3d[,1],INDEX=clusts)+3),
     y=tapply(FUN=mean,X=tsne3d[,2],INDEX=clusts),
     z=(tapply(X=tsne3d[,3],FUN=function(X) {quantile(X,0.99)},INDEX=clusts)+3),
     texts=levels(clusts),col="black",font=2,cex=1)
readline(prompt="Please adjust 3D orientation and press Enter to take a snapshot")
rgl.snapshot("3DtSNE.png","png")

#tSNE with cell cycle states
tempCyc<-cbind(cycScores$cellID,cycScores$phases)
key5<-unlist(lapply(rownames(eb1S@tsne.rot), function(x) tempCyc[which(x==tempCyc[,1]),2]))
layout(cbind(2:1),heights=c(1,7))
par(mar=c(3,3,0,1),mgp=2:0)
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=viridis(3,.8)[c(3,1,2)][as.factor(key5)],
       bg=viridis(3,0.4)[c(3,1,2)][as.factor(key5)])
par(mar=c(0,3,0,1))
plot.new()
legend("top",bty="n",horiz=T,pch=c(NA,21,21,21),
       legend=c("Cell cycle\nphase:",levels(cycScores$phases)),
       col=c(NA,viridis(3)[c(3,1,2)]),pt.bg=c(NA,viridis(3,0.5)[c(3,1,2)]))
readline(prompt="Press Enter to continue: ")

#tSNE Plot with Sample#
transp <- data.frame(bg=rep(0.5,ncol(eb1S@data)),col=rep(1,ncol(eb1S@data)))
ptID<-rownames(eb1S@tsne.rot)
ptID1<-ptID2<-NULL
for (i in 1:length(ptID)){
  ptID1<-c(ptID1,unlist(strsplit(ptID[i],"_"))[1])
  ptID2<-c(ptID2,unlist(strsplit(ptID[i],"_"))[2])
}
ptNames<-names(table(ptID1))
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste("Labelled by Sample#",
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=16,
       col=as.factor(ptID1))
legend("topleft",bty="n",pch=c(16,16,16,16),
       legend=ptNames,
       col=as.factor(ptNames))
IDOut<-cbind(rownames(eb1S@tsne.rot),ptID1,ptID2,clusts,names(table(cycScores$phases))[as.numeric(key5)])
colnames(IDOut)<-c("CellName","Sample","Cell#","Cluster#","CyclePhase")
write.table(IDOut,"Cell_clusterID_cycle.txt",sep = "\t",row.names = FALSE)
readline(prompt="Press Enter to continue: ")

#heatmap of top 15 genes per cluster compared to all others
deG<-get(load(paste0(timePoint,"_precalc_",gsub(".","",names(minDEgenes)[f3],fixed=T),"_deGall.RData")))
heatGenes <- lapply(deG,function(X) rownames(X[order(X$fwer),])[1:15])
clustMeans <- sapply(CGS,function(X) X$MTC[X$genes %in% unique(unlist(heatGenes))])
rownames(clustMeans) <- CGS[[1]]$genes[CGS[[1]]$genes %in% unique(unlist(heatGenes))]
hG <- hclust(dist(clustMeans),"complete")
hC <- hclust(dist(t(clustMeans)),"single")
tempLabCol <- paste(paste0("#",seq_along(deG)),
                    paste(sapply(deG,nrow),"DE"),sep=": ")
clustMeans2 <- sapply(CGS,function(X) X$MTC[X$genes %in% unique(markers[,2])])
rownames(clustMeans2) <- CGS[[1]]$genes[CGS[[1]]$genes %in% unique(markers[,2])]
colnames(clustMeans2) <- paste0("Clust #",seq_along(deG))
ord1<-hclust(dist(t(clustMeans2)))
clustMeans3<-clustMeans2[,ord1$order]
Labels3<-cbind(rep("Col Labels",ncol(clustMeans3)),colnames(clustMeans3))
heatmap.2(clustMeans,Rowv=as.dendrogram(hG),Colv=as.dendrogram(hC),scale="row",
          col="viridis",trace="none",
          ColSideColors=clustCols,labCol=tempLabCol)
readline(prompt="Press Enter to continue: ")

#output a list for top 20 genes per comparison
CL <- makeCluster(detectCores() - 2)
tempCL <- clusterEvalQ(CL,library(Seurat))
eb2S<-eb1S
eb2S@ident <- as.factor(as.integer(eb1S@data.info[,f3+3]))
names(eb2S@ident) <- rownames(eb2S@data.info)
clusterExport(CL,"eb2S")
deGall2 <- parLapplyLB(cl=CL,levels(eb2S@ident),function(X) {
  temp <- FindMarkers(eb2S,X)
  temp$fwer <- p.adjust(temp$p_val,"holm",nrow(eb2S@data))
  return(temp) 
})
deGvs2 <- parApply(CL,combn(levels(eb2S@ident),2),2,function(cl) {
  tempOut <- FindMarkers(eb2S,ident.1=cl[1],ident.2=cl[2],thresh.use = -1, min.pct = -1,min.diff.pct = -1)
  tempOut$fwer <- p.adjust(tempOut$p_val,"holm",n=nrow(eb2S@data))
  tempOut$clustPair <- factor(paste(cl,collapse="vs"))
  tempOut$gene <- rownames(tempOut)
  return(tempOut)
})
stopCluster(CL)
for (i in 1:length(deGall2)){
  tmpData<-deGall2[[i]][order(deGall2[[i]]$avg_diff,decreasing=TRUE),]
  tmpData<-cbind(rownames(tmpData),tmpData)
  tmpData<-tmpData[tmpData$p_val<=0.05,] #only significant changes
  Mean<-apply(tmpData,1, function(X) CGS[[i]]$MTC[which(CGS[[i]]$genes == X[1])][1])
  tmpData<-cbind(tmpData,Mean)
  colnames(tmpData)[1]<-"Genes"
  write.table(tmpData,paste0("SigGenes_Cluster",i,".txt"), sep = "\t",row.names = FALSE)
  }

#List of top genes most unique to the cluster
sigGvs<-lapply(deGvs2,function(x){return(x[x$p_val<0.05,])})
deGvs<-NULL
for (i in 1:length(sigGvs)){
  deGvs<-rbind(deGvs,sigGvs[[i]])
}
cls<-table(clusts)
groups <- matrix(unlist(strsplit(as.character(deGvs[,6]),"vs")),nrow = 2)
for (i in 1:length(cls)){
  tmp1a<-deGvs[which(groups[1,]==names(cls)[i]),]
  tmp1b<-deGvs[which(groups[2,]==names(cls)[i]),]
  tmp1<-rbind(tmp1a,tmp1b)
  idx<-c(rep(1,nrow(tmp1a)),rep(-1,nrow(tmp1b)))
  tmp1<-cbind(tmp1,idx)
  tmp2<-table(tmp1$gene)
  tmp2<-tmp2[order(tmp2,decreasing = T)]
  tmp2<-tmp2[tmp2>0.7*(length(cls)-1)]
  tmp3<-tmp1[tmp1$gene%in%names(tmp2),]
  tmp3$avg_diff<-tmp3$avg_diff*tmp3$idx
  tmp4<-aggregate(avg_diff~gene,data=tmp3,function(x) {return(round(100*sum(x>0)/length(x),0))})
  tmp5<-aggregate(avg_diff~gene,data=tmp3,mean)
  tmp6<-aggregate(p_val~gene,data=tmp3,max)
  tmp7<-aggregate(fwer~gene,data=tmp3,max)
  clustMeans2 <- sapply(CGS,function(X) {apply(tmp7,1,function(t) return(X$MTC[X$gene==t[1]]))})
  rGvs<-cbind(tmp5,tmp6[,2],tmp7[,2],tmp4[,2],clustMeans2)
  colnames(rGvs)<-c("Gene","Ave_diff","maxP_val","maxFWER","Percent>0",names(cls))
  write.table(rGvs,paste0("SigGenes_uniqCluster_Aggr",i,".txt"), sep = "\t",row.names = FALSE)
}

#Context plot with Cluster#
densGene<-list()
densUMI<-list()
densGene[[1]] <- density(pDat$total_features)
densUMI[[1]] <- density(pDat$total_counts)
for (i in 1:length(levels(clusts))){
  densGene[[i+1]] <- density(pDat[which(clusts==i),"total_features"])
  densUMI[[i+1]] <- density(pDat[which(clusts==i),"total_counts"])
}

for (i in 1:length(levels(clusts))){
   layout(matrix(c(2,1,4,3),2),c(6,2),c(2,6))
   par(mar=c(3,3,0,0),mgp=2:0)
   plot(total_features~total_counts,data=pDat,
               pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
               xlab="Library Size",ylab="Genes Detected")
   points(total_features~total_counts,data=pDat[which(clusts==i),],
                   pch=21,col=alpha("red",0.4),bg=alpha("red",0.2),cex=1.2)
   par(mar=c(0,3,1,0))
   plot(x=NULL,y=NULL,ylab="Density",xaxt="n",
               xlim=range(pDat$total_counts),
               ylim=c(min(sapply(densUMI,function(X) min(X$y))),
                                    max(sapply(densUMI,function(X) max(X$y)))))
   lines(densUMI[[1]],col="black",lwd=3)
   lines(densUMI[[i+1]],col="red",lwd=3)
   par(mar=c(3,0,0,1))
   plot(x=NULL,y=NULL,xlab="Density",yaxt="n",
               xlim=c(min(sapply(densGene,function(X) min(X$y))),
                                    max(sapply(densGene,function(X) max(X$y)))),
               ylim=range(pDat$total_features))
   lines(x=densGene[[1]]$y,y=densGene[[1]]$x,col="black",lwd=3)
   lines(x=densGene[[i+1]]$y,y=densGene[[i+1]]$x,col="red",lwd=3)
   par(mar=c(1,1,1,1))
   plot.new()
   legend("center",bty="n",horiz=F,pch=c("-","-"),
          legend=c("Total",paste0("Cluster",i)),
          col=c("black","red"),pt.bg=c("black","red"), lwd = 3)
   readline(prompt="Press Enter to continue: ")
   }

#Gene expression distribution
lg<-"y"
r1<-NULL
layout(cbind(1:2),heights=c(8,4))
while(lg=="y"){
goi1 <- readline(prompt="Please enter the name of the gene for examination: ")
  if (!goi1%in%rownames(eb1S@data)) {
   cat("Gene not found \n")
  } else {
    par(mar=c(2,3,2,3))
    iH <- 101 - cut(eb1S@data[goi1,],breaks = 100,labels=F)
    iH2 <- 5 - cut(eb1S@data[goi1,],breaks = 4,labels=F)
    plot(eb1S@tsne.rot[order(iH,decreasing=T),],xlab=NA,ylab=NA,pch=21,
         col=viridis(100,0.5)[sort(iH,decreasing=T)],
         bg=viridis(100,0.3)[sort(iH,decreasing=T)],
         main=paste0("Distribution for ",goi1," Gene"))
    readline(prompt="Press Enter to continue: ")
    par(mar=c(3,3,0,3))
    plot(x=NULL,y=NULL,xlim=c(0,length(levels(clusts))+1),
         ylim=range(eb1S@data[goi1,]),
         xlab=paste(goi1,"gene Expression by Cluster"),ylab=NA,xaxt="n")
    mtext(levels(clusts),side=1,line=0,at=1:length(levels(clusts)))
    for (i in 1:length(levels(clusts))) {
      if (sum(eb1S@data[goi1,clusts == i])!=0) {
        vioplot(eb1S@data[goi1,clusts == i], at = i,add=T,col=clustCols[i])
        } else { 
          text(x=i,y=0,labels ="_",col=clustCols[i])
        }
      }
    readline(prompt="Press Enter to continue: ")
    par(mar=c(2,3,2,3))
    plot(eb1S@tsne.rot[order(iH2,decreasing=T),],xlab=NA,ylab=NA,pch=21,
         col=viridis(4,0.5)[sort(iH2,decreasing=T)],
         bg=viridis(4,0.3)[sort(iH2,decreasing=T)],
         main=paste0("Relative Distribution for ",goi1," Gene"))
    plot.new()
    legend("top",bty="n",horiz=T,pch=c(21,21,21,21),
           legend=c("High","Med","Low","Bottom"),
    col=viridis(4,0.5)[1:4],
    pt.bg=viridis(4,0.3)[1:4])
    q1 <- readline(prompt=paste0("Save the gene ",goi1," for report? (y/n) "))
    if(q1=="y"){r1<-c(r1,goi1)}
  }
lg<-readline(prompt="Examine another gene? (y/n) ")
}

#output to pdf
pdf(paste0("VisResult",timePoint,".pdf"))
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste("Labelled by Sample#",
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=16,
       col=as.factor(ptID1))
legend("topleft",bty="n",pch=c(16,16,16,16),
       legend=ptNames,
       col=as.factor(ptNames))
plot(x=numClust,y=unlist(minDEgenes),type="b",cex=1.2,
     xlab="Number of clusters",
     ylab="DE genes (@FWER <= 1e-2) b/w most similar clusters")
abline(h=seq(0,max(unlist(minDEgenes)),10),lty=3,col=alpha(1,0.3))
points(x=numClust[f3],y=unlist(minDEgenes)[f3],
       pch=16,cex=1.5,col="red")
plot(tempSil,beside=T,border=NA,main=NA,col=clustCols,do.n.k=T)
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste("Labelled by Cluster#",
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=alpha(clustCols[clusts],transp$col),
       bg=alpha(clustCols[clusts],transp$bg))
text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],INDEX=clusts),
     y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],INDEX=clusts),
     labels=levels(clusts),col="black",font=2,cex=1.5)
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     main=paste("Labelled by Potential Cell Type",
                paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                sep="\n"),
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=alpha(clustCols[clusts],transp$col),
       bg=alpha(clustCols[clusts],transp$bg))
text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],INDEX=clusts),
     y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],INDEX=clusts),
     labels=clusterID,col="black",font=2,cex=0.8)
layout(cbind(2:1),heights=c(1,7))
par(mar=c(3,3,0,1),mgp=2:0)
plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
     xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
points(eb1S@tsne.rot,pch=21,
       col=viridis(3,.8)[c(3,1,2)][as.factor(key5)],
       bg=viridis(3,0.4)[c(3,1,2)][as.factor(key5)])
par(mar=c(0,3,0,1))
plot.new()
legend("top",bty="n",horiz=T,pch=c(NA,21,21,21),
       legend=c("Cell cycle\nphase:",levels(cycScores$phases)),
       col=c(NA,viridis(3)[c(3,1,2)]),pt.bg=c(NA,viridis(3,0.5)[c(3,1,2)]))
heatmap.2(clustMeans,Rowv=as.dendrogram(hG),Colv=as.dendrogram(hC),scale="row",
          col="viridis",trace="none", cexCol = 0.5,
          ColSideColors=clustCols,labCol=tempLabCol)
for (i in markers2){
  tempG<-markers[which(markers[,1]==i),2]
  if(length(tempG)>1){
    tempH<-clustMeans3[rownames(clustMeans3)%in%tempG,]
    ord2<-hclust(dist(tempH))
    tempH<-tempH[ord2$order,]
    heatmap.2(tempH,dendrogram = "none",Colv=FALSE,Rowv=FALSE,scale="row",col="viridis",trace="none",labCol=colnames(clustMeans3),main = paste0("Identifier Genes for ",i))
    tempLab3<-cbind(rep(i,nrow(tempH)),rownames(tempH))
    Labels3<-rbind(Labels3,tempLab3)
  }
}
write.table(Labels3,"CellMarkersHeatMapLabel.txt",sep = "\t",row.names = FALSE)
for (i in 1:length(levels(clusts))){
  layout(matrix(c(2,1,4,3),2),c(6,2),c(2,6))
  par(mar=c(3,3,0,0),mgp=2:0)
  plot(total_features~total_counts,data=pDat,
       pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
       xlab="Library Size",ylab="Genes Detected")
  points(total_features~total_counts,data=pDat[which(clusts==i),],
         pch=21,col=alpha("red",0.4),bg=alpha("red",0.2),cex=1.2)
  par(mar=c(0,3,1,0))
  plot(x=NULL,y=NULL,ylab="Density",xaxt="n",
       xlim=range(pDat$total_counts),
       ylim=c(min(sapply(densUMI,function(X) min(X$y))),
              max(sapply(densUMI,function(X) max(X$y)))))
  lines(densUMI[[1]],col="black",lwd=3)
  lines(densUMI[[i+1]],col="red",lwd=3)
  par(mar=c(3,0,0,1))
  plot(x=NULL,y=NULL,xlab="Density",yaxt="n",
       xlim=c(min(sapply(densGene,function(X) min(X$y))),
              max(sapply(densGene,function(X) max(X$y)))),
       ylim=range(pDat$total_features))
  lines(x=densGene[[1]]$y,y=densGene[[1]]$x,col="black",lwd=3)
  lines(x=densGene[[i+1]]$y,y=densGene[[i+1]]$x,col="red",lwd=3)
  par(mar=c(1,1,1,1))
  plot.new()
  legend("center",bty="n",horiz=F,pch=c("-","-"),
         legend=c("Total",paste0("Cluster",i)),
         col=c("black","red"),pt.bg=c("black","red"), lwd = 3)
}
par(mfrow=c(1,1))
par(mar=c(3,3,3,1))
for (goi1 in r1){
iH <- 101 - cut(eb1S@data[goi1,],breaks = 100,labels=F)
iH2 <- 5 - cut(eb1S@data[goi1,],breaks = 4,labels=F)
plot(eb1S@tsne.rot[order(iH,decreasing=T),],xlab=NA,ylab=NA,pch=21,
     col=viridis(100,0.5)[sort(iH,decreasing=T)],
     bg=viridis(100,0.3)[sort(iH,decreasing=T)],
     main=paste0("Distribution for ",goi1," Gene"))
plot(eb1S@tsne.rot[order(iH2,decreasing=T),],xlab=NA,ylab=NA,pch=21,
     col=viridis(4,0.5)[sort(iH2,decreasing=T)],
     bg=viridis(4,0.3)[sort(iH2,decreasing=T)],
     main=paste0("Relative Distribution for ",goi1," Gene"))
legend("topleft",bty="n",pch=c(21,21,21,21),
       legend=c("High","Med","Low","Bottom"),
       col=viridis(4,0.5)[1:4],
       pt.bg=viridis(4,0.3)[1:4])
plot(x=NULL,y=NULL,xlim=c(0,length(levels(clusts))+1),
     ylim=range(eb1S@data[goi1,]),
     xlab=paste(goi1,"gene Expression by Cluster"),ylab=NA,xaxt="n")
mtext(levels(clusts),side=1,line=0,at=1:length(levels(clusts)))
for (i in 1:length(levels(clusts))) {
  if (sum(eb1S@data[goi1,clusts == i])!=0) {
    vioplot(eb1S@data[goi1,clusts == i], at = i,add=T,col=clustCols[i])
  } else { 
    text(x=i,y=0,labels ="_",col=clustCols[i])
  }
}}

genes1<-unlist(cellMarkers)
for (i in names(table(ptID1))){
  data1<-eb1S@data[,ptID1==i]
  TSNE1<-eb1S@tsne.rot[ptID1==i,]
  for (goi1 in genes1){
    if (goi1%in%rownames(eb1S@data)){
      iH3<-5-cut(data1[goi1,],breaks = 4,labels=F)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           main=paste0("Relative Distribution for ",goi1," Gene in ",i),
           xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
      points(TSNE1[order(iH3,decreasing=T),],pch=21,
             col=viridis(4,0.5)[sort(iH3,decreasing=T)],
             bg=viridis(4,0.3)[sort(iH3,decreasing=T)]
      )
      legend("topleft",bty="n",horiz=T,pch=c(21,21,21,21),
             legend=c("High","Med","Low","Bottom"),
             col=viridis(4,0.5)[1:4],
             pt.bg=viridis(4,0.3)[1:4])
    }
  }
}

n<-1
keep1<-genes1[genes1%in%rownames(eb1S@data)]
data3<-eb1S@data[match(keep1,rownames(eb1S@data)),]
for (j in names(table(clusts))){
  data2<-data3[,clusts==j]
  geneCols<-gg_colour_hue(nrow(data2))
  ncex<-1
  if(nrow(data2)>40){ncex<-40/nrow(data2) }
  plot(x=NULL,y=NULL,xlim=c(0,nrow(data2)+1),
       ylim=c(-2,max(data2)+0.5),
       main = paste("Gene Expressions for Cluster #",j," Predicted as: ",clusterID[n]),ylab=NA,xaxt="n")
  text(1:nrow(data2),-1,labels = rownames(data2),srt=90,pos = 1,offset = 0,cex = ncex)
  for (i in 1:nrow(data2)) {
    if (sum(data2[i,])!=0) {
      vioplot(data2[i,], at = i,add=T,col=geneCols[i])
    } else { 
      text(x=i,y=0,labels ="_",col=geneCols[i])
    }
  }
  boxplot(t(data2),col=geneCols,las=2,main=paste("Gene Expressions for Cluster #",j," Predicted as: ",clusterID[n]))
  n<-n+1
  }

dev.off()

#Output rnk files
q2<-readline(prompt="Export rnk files for all cluster comparisons? (y/n): ")
if(q2=="y"){
  for (i in 1:length(deGvs2)){
    tmp8<-cbind(deGvs2[[i]]$gene,-sign(deGvs2[[i]]$avg_diff),deGvs2[[i]]$p_val)
    tmpMin<-min(as.numeric(tmp8[tmp8[,3]!=0,3]))
    tmp8[which(tmp8[,3]==0),3]<-tmpMin
    tmp8<-cbind(tmp8[,1],round(as.numeric(tmp8[,2])*log10(as.numeric(tmp8[,3])),4))
    write.table(tmp8,paste0(deGvs2[[i]]$clustPair[1],".rnk"),sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  }
  for (i in 1:length(deGall2)){
    tmp9<-cbind(rownames(deGall2[[i]]),-sign(deGall2[[i]]$avg_diff),deGall2[[i]]$p_val)
    tmpMin<-min(as.numeric(tmp9[tmp9[,3]!=0,3]))
    tmp9[which(tmp9[,3]==0),3]<-tmpMin
    tmp9<-cbind(tmp9[,1],round(as.numeric(tmp9[,2])*log10(as.numeric(tmp9[,3])),4))
    write.table(tmp9,paste0("Cluster",i,"vsAll.rnk"),sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  }
}