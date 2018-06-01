#GSVA for samples, please keep the .gmt file in the same folder
require(GSVA)
require(parallel)
require(GSA)
require(qvalue)

#Input data with gene symbol as the first column (name ID) and sample names as the first row
f1<-readline(prompt="Please enter the full name (including .txt, .csv, etc) of the gene expression file (gene ID 1st column, sample name 1st row): ")
datam <- read.table( f1, header = TRUE, sep = "\t", quote="\"",  stringsAsFactors = FALSE)
f3<-readline(prompt="Is the Data: 1. normalized microarray/RNA-seq; or 2. RNA-seq raw counts (Input 1 or 2): ")
if (f3=="2") {rseq=TRUE} else {rseq=FALSE}
f4<-readline(prompt="Please enter the full name (including .txt, .csv, etc) of the Gene Type file (gene ID 1st column, gene symbol 2nd column, and gene type in the 3rd column): ")
genetype <-  read.table( f4, header = TRUE, sep = "\t", quote="\"",  stringsAsFactors = FALSE)
proteingene<-genetype[which(genetype[,3]=="protein_coding"),1:2]
#remove all duplicated genes
dupgene<-proteingene[duplicated(proteingene[,2]),2]
proteingene<-proteingene[which(!proteingene[,2]%in%dupgene),]

#Only protein coding and unique genes used
datau<-datam[which(datam[,1]%in%proteingene[,1]),]
datau<-datau[!duplicated(datau[,1]),]
rownames(datau)<-proteingene[match(datau[,1],proteingene[,1]),2]
datau<-as.matrix(datau[,2:ncol(datau)])

#Input the file with subtype information for each sample, first column contains sample name, second column records group name 
f2<-list.files(pattern = ".gmt")
gmt1<-GSA.read.gmt(f2)
gmt2<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

#GSVA and output
result1<-gsva(datau, gmt2, min.sz=10, max.sz=200, mx.diff=TRUE, verbose=F, rnaseq=rseq, parallel.sz=0)
result2<-NULL
for (i in 1:nrow(result1)){
tmp1<-unlist(strsplit(gmt1[[2]][match(rownames(result1)[i],gmt1[[3]])],"%"))[3]
tmp2<-paste(unlist(gmt2[[match(rownames(result1)[i],names(gmt2))]]),collapse = ",")
result2<-rbind(result2,c(tmp1,tmp2))
}

dir.create("GSVA_Gprofile")
for (j in 1:ncol(result1)){
test1<-pnorm(-abs(scale(result1[,j])))
tmp3<-qvalue(test1,lambda=seq(0.05,0.45,0.01))$lfdr
tmp4<-cbind(result2[,1],rownames(result1),test1,tmp3,sign(result1[,j]),result2[,2])
colnames(tmp4)<-c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")
write.table(tmp4,paste0("./GSVA_Gprofile/",colnames(result1)[j],"_gProfile.txt"),sep = "\t",row.names = F,quote = F)
write.table(tmp4[tmp4[,5]>0,],paste0("./GSVA_Gprofile/",colnames(result1)[j],"_gProfilePos.txt"),sep = "\t",row.names = F,quote = F)
}
write.csv(result1,file=paste0(f1,"_GSVAprotein_result.csv"))