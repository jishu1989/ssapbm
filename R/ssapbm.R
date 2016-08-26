#Subsample analysis of pathways-based methods
ssapbm <- function(data=NULL, pathways=NULL, ref=NULL, minp=2, maxp=9999,method= c("sumoftsq","ht2", "GSEA", "GAGE","GSA"), perm=T, sampling=c("sample.labels","gene.labels"),B=100, sample.size=NULL, steps=10, fdr=NULL,thr=.05,dc=T,rep=F,cr=NULL){
	d <- processdata(data, pathways)
	data <- d[[1]]
	pathways <- d[[2]]
	cls1 <- ref
	cls2 <- c(1:ncol(data))[-ref]
	TP <- pvl <- fdr.tp <- NULL
	if(is.null(fdr)){
		fdr <- TRUE
		fdr.method <- "BH"
	}
	else{
		fdr.method <- fdr
		fdr <- TRUE
	}
	print(fdr.method)
	if(method[1]=="sumoftsq"){
		TP <- p.squared.t.test(data, cls1, cls2,steps=B, sampling=sampling[1], pathways=pathways)
	 	if(!is.null(sample.size)){
			
			pvl <- subsampanalysis(data, cls1, cls2, pathways, B=100, steps, sample.size=sample.size, method="sumoftsq", sampling=sampling[1],rep=rep, rn=rownames(TP),dc=dc,thr=thr)
			pvtmp <- c()
			for(m in 1:length(pvl)){
				pvtmp <- cbind(pvtmp, pvl[[m]])
			}
			pvl <- pvtmp
			rm(pvtmp)
		}	 
		if(fdr){
			print("hi")
			print(fdr.method)
			fdr.tp <- p.adjust(TP[,"pval"], method=fdr.method)
		}
	
	}
	if(method[1]=="ht2"){
		TP <- pval.htsq(data, cls1, cls2,pathways=pathways, shrinkage=shrinkage, perm=perm, sampling=sampling[1], steps=B)
	 	TP <- cbind(TP[[1]],TP[[2]])
		colnames(TP) <- c("ht2","pval")
		if(!is.null(sample.size)){
			
			pvl <- subsampanalysis(data, cls1, cls2, pathways, B=100, steps, sample.size=sample.size, method="sumoftsq", sampling=sampling[1],rep=rep, shrinkage=shrinkage, rn=rownames(TP), dc=dc,thr=thr)
			pvtmp <- c()
			for(m in 1:length(pvl)){
				pvtmp <- cbind(pvtmp, pvl[[m]])
			}
			pvl <- pvtmp
			rm(pvtmp)
			#pvl <- as.matrix(as.data.frame(pvl))
		}	 
		if(fdr){
			fdr.tp <- p.adjust(TP[,2], method=fdr.method)
		}
	
	}
	if(method[1]=="GSEA"){
		clsn <- rep(0, length(cls1)+length(cls2))
		clsn[cls2] <- 1
		lp<- sapply(pathways, length)
		lp <- c(min(lp), max(lp))
		TP <- GSEAx(data, ref=clsn, pathways, reshuffling.type=sampling[1], nperm=B, gs.size.threshold.min = lp[1], gs.size.threshold.max=lp[2])
		TP <- as.matrix(rbind(TP[[1]][,c("ES","NOM p-val")], TP[[2]][,c("ES","NOM p-val")]))
		print(TP)	
		tmp1 <- cbind(as.numeric(as.vector(TP[,1])), as.numeric(as.vector(TP[,2])))
		rownames(tmp1) <- rownames(TP)
		colnames(tmp1) <- colnames(TP)
		TP <- tmp1
		TP <- TP[names(pathways),]
		colnames(TP) <- c("GS", "pval")

		if(!is.null(sample.size)){
			pvl <- subsampanalysis(data, cls1, cls2, pathways, B=100, steps, sample.size=sample.size, method="GSEA", sampling=sampling[1],rep=rep, rn=rownames(TP),dc=dc,thr=thr)
		}	 
		if(fdr){
			fdr.tp <- p.adjust(TP[,"pval"], method=fdr.method)
		}
		 
	}
	if(method[1]=="GAGE"){
		lp<- sapply(pathways, length)
		lp <- c(min(lp), max(lp))
		print(head(data))
		pvl1 <- gage(data, gsets = pathways, ref = cls1, samp = cls2, compare = "unpaired", set.size=c(lp[1],lp[2]))
		pvl1 <- lapply(pvl1[1:2], function(x)x[,1:5])
		pvx <- lapply(pvl1[1:2], function(x)x[,"p.val"])
		pvx <- cbind(pvx[[1]], pvx[[2]][names(pvx[[1]])], names(pvx[[1]]))
		#tx <- apply(pvx,1,which.min)
		tx <- t(apply(pvx,1,function(x)c(which.min(as.numeric(x[1:2])),x[3])))		
		TP <- t(apply(tx,1, function(x)pvl1[[as.numeric(x[1])]][x[2],]))
		TP <- TP[,2:4]
		colnames(TP) <- c("stat.mean", "pval", "fdr.tp")
		if(!is.null(sample.size)){

                	pvl <- subsampanalysis(data, cls1, cls2, pathways, B=100, steps, sample.size=sample.size, method="GAGE",dc=dc,thr=thr)
                }
	}
	if(method[1]=="GSA"){
		clsn <- rep(1, length(cls1)+length(cls2))
                clsn[cls2] <- 2
		
		lp<- sapply(pathways, length)
		lp <- c(min(lp), max(lp))
		TP<-GSA(data, clsn, genenames=rownames(data), genesets=pathways,  resp.type="Two class unpaired", nperms=B, minsize=lp[1], maxsize=lp[2])
		TP <- cbind(TP$GSA.scores, TP$pvalues.lo, TP$pvalues.hi)
		rownames(TP) <- names(pathways)
		tx <- (apply(TP,1,function(x)c(which.min(as.numeric(x[2:3])))))
		TP <- cbind(TP,tx)
		print(TP)
		TP <- t(apply(TP,1, function(x)x[c(1,(x[4]+1))]))
		colnames(TP) <- c("GSA_score", "pval")
		if(fdr){
			fdr.tp <- p.adjust(TP[,"pval"], method=fdr.method)
		}
		if(!is.null(sample.size)){
                	pvl <- subsampanalysis(data, cls1, cls2, pathways, B=100, steps, sample.size=sample.size, method="GSA",dc=dc,thr=thr)
                }
				
	}
	dcoverall <- corpaths <- NULL
	if(!is.null(cr)){
		crall <- avgcorpath(data, ref=c(cls1,cls2), pathways, method=cr)
		crctrl <- avgcorpath(data, ref=c(cls1), pathways, method=cr)
		crtreat <- avgcorpath(data, ref=c(cls2), pathways, method=cr)
		corpaths <- cbind(crall, crctrl, crtreat)	
	
	}
	if(dc){
		clsn <- rep(1, length(cls1)+length(cls2))
                clsn[cls2] <- 2
		dcoverall <-  pcdetectioncall(data, ref=clsn, pathways=pathways, fdr=fdr.method, thr=thr)
	}
	if(fdr){
		qval<- fdr.tp
		
		sigpath <- names(which(qval<=thr))
		list(cbind(TP, qval), pvl, dcoverall=dcoverall, corpaths=corpaths, thr=thr, pathways=names(pathways),sigpath=sigpath,fdr=fdr.method)	
	}
	#else{
		
	#}
}
subsampanalysis <- function(data, cls1, cls2, pathways, B=100, steps, sample.size=NULL, method=NULL, rep=F, sampling=NULL, shrinkage=T, rn=NULL, fdr="BH", dc=T, thr=0.05){
	plist <- list()
	dcall <- list()
	for(i in 1:length(sample.size)){
		pvec <- c()
		dc_p <- c()
		for(j in 1:steps){
			print("*****")
			print(cls1)
			smp1 <- sample(cls1, sample.size[i], rep=rep)
			smp2 <- sample(cls2, sample.size[i], rep=rep)
			if(method=="sumoftsq"){
				tmp <- p.squared.t.test(data[,c(smp1,smp2)], c(1:length(smp1)), c((length(smp1)+1):(length(smp2)+length(smp1))), steps=B, sampling=sampling[1], pathways=pathways) 
				pvec <- cbind(pvec, tmp)
			}
			if(method=="ht2"){
				tmp <- pval.htsq(data[,c(smp1,smp2)], c(1:length(smp1)), c((length(smp1)+1):(length(smp2)+length(smp1))), steps=B, sampling=sampling, perm=perm, shrinkage=shrinkage)
				pvec <- cbind(pvec, tmp[[2]])
			}
			if(method=="GSEA"){
				clsn <- rep(0, length(smp1)+length(smp2))
                		clsn <- c(rep(0, length(smp1)), rep(1, length(smp2)))
                		lp<- sapply(pathways, length)
                		lp <- c(min(lp), max(lp))
				print("xxxxxxsdfsdfsdf")
				print(smp1)
				
                		tmp <- GSEAx(data[,c(smp1, smp2)], ref=clsn, pathways, reshuffling.type=sampling[1], nperm=B, gs.size.threshold.min = lp[1], gs.size.threshold.max=lp[2])
				ptmp1 <- as.numeric(as.vector(tmp[[1]][,"NOM p-val"]))
				ptmp2 <- as.numeric(as.vector(tmp[[2]][,"NOM p-val"]))
				tmp1 <- c(ptmp1,ptmp2)
				tmp <- tmp1[rn]	
				pvec <- cbind(pvec, tmp)
			}
			if(method=="GAGE"){
				print(c(1:length(smp1)))
				print(head(data[,c(smp1,smp2)]))
				pvl <- gage(data[,c(smp1,smp2)], gsets = pathways, ref = c(1:length(smp1)), samp = c((length(smp1)+1):(length(smp1)+length(smp2))), compare = "unpaired", set.size=c(2, 9999))
                		pvl <- lapply(pvl[1:2], function(x)x[,1:5])
                		pvx <- lapply(pvl[1:2], function(x)x[,"p.val"])
                		pvx <- cbind(pvx[[1]], pvx[[2]][names(pvx[[1]])], names(pvx[[1]]))
                		tx <- t(apply(pvx,1,function(x)c(which.min(as.numeric(x[1:2])),x[3])))
                		pvn <- t(apply(tx,1, function(x)pvl[[as.numeric(x[1])]][x[2],]))
				pvec <- cbind(pvec, pvn[names(pathways),"p.val"])
				
			}
			if(method=="GSA"){
				clsn <- rep(1, length(smp1)+length(smp2))
                		clsn <- c(rep(1, length(smp1)), rep(2, length(smp2)))
			 	TPx<-GSA(data[,c(smp1, smp2)], clsn, genenames=rownames(data), genesets=pathways,  resp.type="Two class unpaired", nperms=B,minsize=2, maxsize=9999)                     
				TPx <- cbind(TPx$GSA.scores, TPx$pvalues.lo, TPx$pvalues.hi)
                		rownames(TPx) <- names(pathways)
                		tx <- (apply(TPx,1,function(x)c(which.min(as.numeric(x[2:3])))))       
                		TPx <- cbind(TPx,tx) 
				TPx <- (apply(TPx,1, function(x)x[c((x[4]+1))]))
				pvec <- cbind(pvec, TPx[names(pathways)])
			}
			if(dc){
			 dctmp <- pcdetectioncall(data[,c(smp1, smp2)], ref=c(rep(1,length(smp1)), rep(2,length(smp2))), pathways=pathways, fdr=fdr, thr=thr)	
			dc_p <- cbind(dc_p, dctmp[names(pathways)])
			}
			
		}
		plist[[i]] <- pvec	
		if(dc){
			dcall[[i]] <- dc_p		
			
		}
	}
	
	names(plist) <- sample.size
	if(dc){
	names(dcall) <- sample.size
	}
	plist$dcall <- dcall
	invisible(plist)
}

powerpath <- function(x, type="pval"){
	tp <- x[[1]]
	thr <- x$thr
	fdr <- x$fdr
	k <- x$sigpath
	if(length(k)==0){
		stop("no p-value is below threshold")
	}
	resampval <- length(x[[2]][-length(x[[2]])])
	rml <- list()
	tprx <- fprx <- fdrx <- list()
	tpx <- k
	tnx <- setdiff(rownames(tp), k) 
	#if(type=="qval"){
		#bh <- lapply(ff, function(x)apply(x,2,function(y)p.adjust(y,method=xx$fdr)))
	#}
	for(i in 1:resampval){
		p <- x[[2]][[i]]
		p <- apply(p,2,function(x)p.adjust(x,method=fdr))
		ptmp <- apply(p,2,function(x){m<- rep(0, length(x));m[which(x<=thr)]<-1;m})
		rownames(ptmp) <- rownames(p)
		rml[[i]] <- (rowMeans(ptmp[k,]))
		fdrtmp <- fprtmp <- tprtmp <- c()
		for(j in 1:ncol(ptmp)){
			n1 <- rownames(ptmp[which(ptmp[,j]==1),])
			tprtmp <- c(tprtmp, length(intersect(tpx, n1))/length(tpx))
			fprtmp <- c(fprtmp, length(intersect(tnx, n1))/length(tnx))
			fdrtmp <- c(fdrtmp, length(intersect(tnx, n1))/length(n1))
		}	
		tprx[[i]] <- tprtmp
		fprx[[i]] <- fprtmp
		fdrx[[i]] <- fdrtmp
	}
	names(tprx) <- names(fprx) <- names(fdrx) <- names(rml) <- names(xx[[2]][1:resampval])
	list(power=rml, tpr=tprx, fpr=fprx, fdr=fdrx, sigpath=k)
}
pcdetectioncall <- function(data,pvl=NULL, ref=NULL, pathways=NULL, fdr="BH", thr=.05){
	if(is.null(pvl)){
		tb <-unique(ref)
		k1 <- which(ref==tb[1])
		k2 <- which(ref==tb[2])
		print(k1)
		print(k2)
		pvl <- apply(data,1, function(x){t.test(x[k1],x[k2])$p.value})
	}
	cnt <- which(pvl<=thr)
	nm <- numeric(length(pvl))
	names(nm) <- names(pvl)
	fdr.method <- fdr
	fdr <- TRUE
	if(fdr){
		pvl <- p.adjust(pvl, method="BH")
	}
	nm[cnt] <- 1
	diff <- sapply(pathways, function(x){count <- length(which(nm[x]==1));100*(count/length(x))})
	names(diff) <- names(pathways)
	diff
}
avgcorpath <- function(data, ref=NULL, pathways, method="pearson"){
	if(is.null(ref)){ref <- 1:nrow(data)}
	crx <- sapply(pathways, function(x){cr <- t(cor(m[x,ref], method=method)); mean(cr[upper.tri(cr)])})	
	crx
}
plot.ssapbm <- function(xx, data, type=c("power","fdr","fpr","tpr","dc","cor"), val="pval", pathways=NULL,nc=1,nr=1,pos="bottom"){
		print(type)
		pth <- NULL
		if(length(pathways)==0){
			pth <- xx$pathways
		}
		else{pth <- pathways}
		type <- type[1]
		typex=c("power","fdr","fpr","tpr","dc","cor")
		nm <- c("Power", "FDR", "FPR", "TPR", "Detection call", "Correlation")
		names(nm) <- typex
		if((type!="dc")&&(type!="cor")){
		pwr <- powerpath(xx, type=val)
		dat <- pwr[[type[1]]]
		dat <- getdf(dat)
		colnames(dat) <-  c("sample", "stat")
		}
		dftrl <- list()
		if(type=="dc"){
			dcv <- xx[[2]][[length(xx[[2]])]]
			for(i in 1:(length(dcv))){
				#if(i==(length(dcv)+1)){
				#dctmp <- dcv[[i-1]][pth,]
				#}
				#else{
				dctmp <- dcv[[i]][pth,]
				#}
				dftr <- data.frame(rep(rownames(dctmp), ncol(dctmp)), as.vector(dctmp) )
				colnames(dftr) <- c("Pathways", "DC")
				#if(i!=(length(dcv)+1)){	
				xtst <- cbind(xx$dcoverall[pth],pth)
				colnames(xtst) <- c("dcall", "Pathways")
				xtst <- data.frame(xtst, type=factor(rep(1,length(pth))))
				xtst$dcall <- as.numeric(as.vector(xtst$dcall))
				print(xtst)	
				gdc <- ggplot(dftr, aes(x=Pathways, y=DC,fill=Pathways)) + geom_boxplot()+xlab(paste("Pathways","(sample: ",names(dcv)[i],")",sep=""))+ylab("Detection Call")+theme(legend.position="none") +  theme(axis.text.x=element_blank())
				gdc <- gdc +geom_point(data=xtst,aes(x=Pathways, y=dcall, colour = factor(type,labels="Overall detection call")),shape=5, size=2, fill="red")+ ylim(min(c(xtst$dcall,dftr$DC)), max(c(xtst$dcall,dftr$DC)))+theme(legend.title=element_blank())
				#}
				#else{
				#gdc <- ggplot(dftr, aes(x=Pathways, y=DC,fill=Pathways)) + geom_boxplot()
				#gdc <- g_legend(gdc)
				#}
				dftrl[[i]] <- gdc
			}
		}
		if(type=="cor"){
			sig <- rep(0, length(xx$pathways))
			crdat <- data.frame(xx$corpaths,sig)
			crdat[xx$sigpath, "sig"] <- 1
			crdat$sig <- as.factor(crdat$sig)
			p1 <-  ggplot(crdat, aes(x=rank(crall), y=crall, color=factor(sig, labels=c("Differentially expressed", "Not Differentially expressed")) ))+geom_point()+theme(legend.position="none")+ labs(color = "",x="Ordered pathways", y="Average correlation of individual pathways")
			p2 <- ggplot(crdat, aes(x=crctrl, y=crtreat, color=sig))+geom_point()+theme(legend.position="none")
			crd <- getdf(list(ctrl=crdat[xx$sigpath,2], crtreat=crdat[xx$sigpath,3]))
			p3 <- ggplot(crd, aes(x=sample, y=stat, fill=sample))+geom_boxplot()+theme(legend.position="none")
			grid_arrange_shared_legend(list(p1,p2,p3),nrow=nr, ncol=nc,position=pos)
		}	
		if(type=="dc"){
			#multiplot(dftrl, cols=length(dftrl))
			print(length(dftrl))
			grid_arrange_shared_legend(dftrl,nrow=nr, ncol=nc,position=pos)
		}
		if((type!="dc")&&(type!="cor")){
		ggplot(dat, aes(y=stat, x=factor(sort(as.numeric(sample))),fill=sample)) + geom_boxplot()+xlab("Sample-size")+ylab(nm[type[1]])		
		 }
}
getdf <- function(lst=NULL){
	if(!is.null(lst)){
		nm <- names(lst)
		len <- sapply(lst, length)
		lab <- rep(nm, len)
		dfx <- data.frame(sample=lab, stat=unlist(lst),stringsAsFactors=FALSE )
	
	}
}



library(ggplot2)
library(gridExtra)
library(grid)


grid_arrange_shared_legend <- function(plots, ncol = 1, nrow = 1, position = c("bottom", "right")) {
   if(is.null(ncol)||((ncol*nrow)<length(plots))){
   	ncol <- length(plots) 
   	nrow <- 1
    }
   if(is.null(ncol)||((ncol*nrow)<length(plots))){
        nrow <- length(plots)
   	ncol <- 1
   }		

  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

