check_log<-function(x) {
	ex <- exprs(x)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if (LogC) { 
		ex[which(ex <= 0)] <- NaN
	  	exprs(x) <- log2(ex)
	}
	x
}

library(Biobase)
library(GEOquery)
library(limma)

gset<-getGEO("GSE68889", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE68889",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE59445", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE59445",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE64250", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
ann_extr=function(x) {
    unlist(strsplit(x," \\/\\/ "))[2]
}
ann=apply(data.frame(fData(data)[,10]),1,FUN=ann_extr)
agg=aggregate(exprs(data),by=list(ann),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE64250",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE26529", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
ann_extr=function(x) {
    unlist(strsplit(x," \\/\\/ "))[2]
}
ann=apply(data.frame(fData(data)[,10]),1,FUN=ann_extr)
agg=aggregate(exprs(data),by=list(ann),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE26529",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE115406", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE115406",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE43090", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,7])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE43090",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE10927", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE10927",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE18271", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE18271",cor.all,PHOX2B.val,agg)


gset<-getGEO("GSE3678", GSEMatrix=TRUE, getGPL=TRUE)
data <- gset[[1]]
data<-check_log(data)
glist=gsub(" \\/\\/\\/.+","",perl=T,fData(data)[,11])
agg=aggregate(exprs(data),by=list(glist),FUN=mean)
rownames(agg)=agg[,1]
agg=agg[,-1]
agg=agg[-1,]
agg=normalizeBetweenArrays(agg)
PHOX2B<-grep("PHOX2B", rownames(agg))
PHOX2B.val<-as.numeric(agg[PHOX2B,])
mycor<- function(x) {
	c<-cor.test(PHOX2B.val,x,method="spearman")
	c(c$estimate,c$p.value)
}
tmp<-data.frame(t((apply(agg,1,FUN=mycor))))
cor.all<-data.frame("COR"=tmp[,1],"PVAL"=tmp[,2],"ADJP"=p.adjust(tmp[,2],method="BH"))
rownames(cor.all)<-rownames(agg)
save(file="GSE3678",cor.all,PHOX2B.val,agg)

