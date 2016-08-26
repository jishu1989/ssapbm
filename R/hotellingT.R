htsq <- function(data, cls1, cls2, shrinkage=F){
	meanS1 <- apply(data[,cls1],1, mean)
	meanS2 <- apply(data[,cls2],1, mean)
	covS1 <- covS2 <- NULL
	if((nrow(data)-1)>(ncol(data))){
		shrinkage=T
	}
	if(shrinkage){
		covS1 <- cov.shrink(t(data[,cls1]),verbose=F)
		covS2 <- cov.shrink(t(data[,cls2]), verbose=F)
    	}
	else{
		covS1 <- cov(t(data[,cls1]))
		covS2 <- cov(t(data[,cls2]))
    	}
	n1 = ncol(data[,cls1])
    	n2 = ncol(data[,cls2])
	p <- nrow(data)
	
	
  	covS = ((n1 - 1)*covS1 + (n2 - 1)*covS2)/(n1 + n2 - 2)
	covSinv <- solve(covS)
	htsqx <- (t(meanS1-meanS2) %*% covSinv %*%(meanS1-meanS2))*(n1*n2/(n1+n2))

	list(statistic=htsqx, df1=p, df2=n1+n2-p-1, n1=n1,n2=n2)
}
pval.htsq <- function(data, cls1, cls2, pathways=NULL, shrinkage=F, perm=F, sampling="sample.labels", steps=100){
	pval <-  tvec <- permtvec <- c()
	tvec <- c()
		if(!perm){
			for(i in 1:length(pathways)){
				datn <- data[pathways[[i]], ]
				tstat <- htsq(datn, cls1, cls2, shrinkage=shrinkage)
				m <- (tstat$n1 + tstat$n2 - tstat$p - 1)/(tstat$p*(tstat$n1 + tstat$n2 -2))
				ptmp <- 1 - pf(m*tstat$statistic, df1, df2)
				pval <- c(pval, ptmp)
				tvec <- c(tvec, tstat$statistic)
			}
			names(tvec) <- names(pval) <- names(pathways)
			
		}
		else{
			permtemp <- c()
			for(i in 1:length(pathways)){
				tvec <- c(tvec, htsq(data[pathways[[i]],], cls1, cls2, shrinkage=shrinkage)$statistic)
			}
			names(tvec) <- names(pathways)
			for(i in 1:steps){
				tmp <- c()
				if(sampling=="sample.labels"){
					smp <-  sample(1:ncol(data))
					#cls1n <- smp[1:length(cls1)]
					#cls2n <- smp[(length(cls1)+1):ncol(data)]
					for(j in 1:length(pathways)){
						tmp <- c(tmp, htsq(data[pathways[[j]],smp], cls1, cls2, shrinkage=shrinkage)$statistic)
					}
					names(tmp) <- names(pathways)
				}
				if(sampling=="gene.labels"){
					datn <- data[sample(nrow(data)),]
					rownames(datn) <- rownames(data)
			
					for(j in 1:length(pathways)){
						tmp <- c(tmp, htsq(datn[pathways[[j]],], cls1, cls2, shrinkage=shrinkage)$statistic)
					}
					names(tmp) <- names(pathways)

				}	
				permtemp <- cbind(permtemp, abs(tmp))
			}
			htval <- cbind(abs(tvec), permtemp)
        		pval <- apply(htval, 1, function(x){sum(x[2:(steps+1)]>x[1])/steps})
		
		}
		list((tvec), pval)
		
	
	}

