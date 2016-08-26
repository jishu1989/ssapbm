
mat <- matrix(rnorm(1000*40), 1000,40)
rownames(mat) <- paste("g", c(1:nrow(mat)),sep="")
pl <- list()
s=25
b <- 0
for(i in 1:40){
       
	a <- b+1
	b <- b+s
	pl[[i]] <- paste("g", c(a:b),sep="")
	print(a)
	print(b)

}
names(pl) <- paste("pth", c(1:length(pl)), sep="")
mat[1:10,21:40] <- mat[1:10,21:40] + 1 
cls1 <- c(1:20)
cls2 <- c(21:40)
#xx <- sscmp(mat, pathways=pl, ref=1:20, method="sumoftsq", perm=T,sampling="sample.labels", B=1000, sample.size=c(5,10),steps=10)

processdata <- function(data, pathways, min=5, max=500){
	if(is.null(names(pathways))){
		names(pathways) <- paste("pathways", c(1:length(pathways)), sep="")
	}
	if(is.null(colnames(data))){
		colnames(data) <- paste("cl",1:ncol(data), sep="")
	}
	rn <- rownames(data)
	cnt <- 1
	nm <- names(pathways)
	pth <- list()
	nmptmp <- c()
	for(i in 1:length(pathways)){
		tmp <- intersect(rn, pathways[[i]])
		if((length(tmp) >=min)&& (length(tmp) <=max)){
			pth[[cnt]] <- tmp
			nmptmp[cnt] <- nm[i]
			cnt <- cnt + 1
		}
	}
	names(pth) <- nmptmp
	
	
	list(data[unique(unlist(pth)),],pth)
	
}
ttest <- function(data, cls1, cls2, pathways=NULL) { apply(data,1, function(x){t.test(x[cls1], x[cls2])$statistic})}
squared.t.test <- function(data, cls1, cls2, pathways=NULL){
        tval <- ptval <-  c()
        tval <- ttest(data,cls1,cls2)
	for(i in 1:length(pathways)){
		tmp <- c()
		nm1 <- pathways[[i]]
		nm2 <- names(pathways)[i]
		tmp <- tval[nm1]
		tmp <- tmp[!is.na(tmp)]
		tmp <- (sqrt(sum(tmp^2)))/length(tmp)
        	names(tmp) <- nm2
		ptval <- c(ptval, tmp)
		
	}
        return(ptval)
}
resampling.squared.t.test <- function(data,cls1,cls2,steps=100,pathways=NULL, sampling="sample.labels"){
        tvec <- c()
	if(sampling=="sample.labels"){
        	for(i in 1:steps){
                	smp <- sample(ncol(data))
			cls1n <- smp[1:length(cls1)]
			cls2n <- smp[(length(cls1)+1):(length(cls1)+length(cls2))]
                	tmp <- squared.t.test(data, cls1n, cls2n, pathways=pathways)
			tvec <- cbind(tvec, tmp)
        	}
	}
	if(sampling=="gene.labels"){
		for(i in 1:steps){
			smp <- sample(1:nrow(data))
			datax <- data[smp, ]
			rownames(datax) <- rownames(data)
			tmp <- squared.t.test(datax, cls1, cls2, pathways=pathways)
			tvec <- cbind(tvec, tmp)
		}
		    	
	}
	
        return(tvec)
}
p.squared.t.test <- function(data,cls1, cls2, steps=100, sampling="sample.labels", pathways=NULL){
        tv <- squared.t.test(data, cls1, cls2, pathways=pathways)
	tvec <- resampling.squared.t.test(data, cls1, cls2, steps=steps, pathways=pathways, sampling=sampling)
        tvec <- cbind(tv, tvec)
	pval <- apply(tvec, 1, function(x){sum(x[2:(steps+1)]>x[1])/steps})
        cbind(tv,pval)

}


