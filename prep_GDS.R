############## DATASET PREPARATION ###################
library(Biobase)
library(limma)
library(GEOquery)

df.multi<-function(x,y,z) {
	d<-aggregate(exprs(x),by=list(fData(x)[,y]),mean)
	d<-d[-1,]
	newx<-NULL
	rows<-NULL
	for (i in 1:dim(d)[1]) {
		k<-d[i,1]
		l<-unlist(strsplit(k,"\\/\\/\\/"))
		for (j in 1:length(l)) {
			rows<-c(rows,l[j])
			newx<-rbind(newx,d[i,-1])
		}
		if (i %% 1000 == 0) {
			cat(i,dim(newx)[1],"\n")
		}
	}
	newx<-aggregate(newx,by=list(rows),mean)
	rownames(newx)<-newx[,1]
	newx[,-1]
}

gds<-c("GDS1713","GDS2609","GDS5367","GDS2860","GDS3898","GDS3787","GDS3501","GDS4011","GDS3557","GDS3513","GDS1384","GDS2705","GDS1813","GDS4362","GDS1869","GDS4006","GDS4104","GDS5408","GDS5615","GDS2763","GDS5298","GDS2414","GDS1972","GDS1646","GDS3029","GDS5263","GDS2852","GDS1926")

for(ds in 1:length(gds)) {
#	if(file.exists(gds[ds])) {
#		next
#	}
	if(file.exists(paste0("tmp/",gds[ds]))) {
		load(paste0("tmp/",gds[ds]))
	} else {
		cat("Downloading",gds[ds],"\n") 
		gset <- getGEO(gds[ds], GSEMatrix =TRUE,getGPL=TRUE)
		g<-NULL
		if (grepl("GDS",gds[ds])==T) { 
			g<-GDS2eSet(gset,do.log2=FALSE)
		} else {
			g<-gset[[1]]
		}
		#this is taken from GEO2R, needed to understand is a dataset has to be logged or not
		cat("Preprocessing",gds[ds],"\n") 
		ex <- exprs(g)
		qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
		LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
		if (LogC) { 
			ex[which(ex <= 0)] <- NaN
		  	exprs(g) <- log2(ex)
		}
		cat("Preaggregating probes\n") 
		m<-df.multi(g,3)
	}
	cat("Normalizing probes\n") 
	m<-normalizeBetweenArrays(m)
	v<-match("PHOX2B",rownames(m))
	if(is.na(v)) {
		cat(ds,gds[ds],"cannot find PHOX2B annotation\n",file="PHOX2B.log",append=T);
		next
	}
	cat(ds,gds[ds],"good\n",file="PHOX2B.log",append=T);
	val<-as.numeric(m[v,])
	cors<-NULL
	ps<-NULL
	genes<-NULL
	cat("Computing correlations\n") 
	for(i in 1:dim(m)[1]) {
		if(sum(ifelse(is.finite(as.numeric(m[i, ])),1,0))>=length(m[i,])*0.5) {
			t<-cor.test(as.numeric(m[i,]),val,method="spearman")
			ps<-c(ps,as.numeric(t$p.value))
			cors<-c(cors,as.numeric(t$estimate))
		} else {
			ps<-c(ps,NA)
			cors<-c(cors,NA)	
		}
		genes<-c(genes,rownames(m)[i])
	}
	res<-data.frame("COR"=cors,"PVAL"=ps,"ADJP"=p.adjust(ps,method="BH"),row.names=genes)
	res.opt<-res[as.numeric(res[,3])<=0.05,]
	save(file=gds[ds],m,val,res,res.opt)
}
