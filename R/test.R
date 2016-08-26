p <- sample(unlist(x), 500)

p <- unique(unlist(x))
pnew <- intersect(p1,p)
pnew <- p
nm <- c()
ptog <- list()

for (i in 1:length(pnew)){
	
	 px <-sapply(x, function(y){intersect(y,pnew[i])})
	 pl <- sapply(px, length)
	 k <- which(pl>0)
	ptog[[i]] <- names(pl[k])

}



simdiffx <- function(d, pc=1, diff=.1){
    nr <- nrow(d)*pc/100
    nr <- round(nr)
    #diff <- seq(.1, .4, by=.3/(nr-1))
    print(diff)
    for(i in 1:nr){
        d[i,] <- d[i, ]- diff
        k <- which(d[i,]<0)
        print(length(k))
        if(length(k)>0){
            tmp <- rbeta(length(k),.1,200000)
            #d[i,k] <- 1+d[i,k]
            d[i,k] <- tmp
        }
    }
    d
    
}

m <- matrix(rnorm(80000), 2000, 40 )
m[1:1000, 21:40] <- m[1:1000, 21:40] + .75
rownames(m) <- paste("g", c(1:2000),sep="")
pathways <- list()
s1 <- 1
s2 <- 50
for(i in 1:40){
	
	pathways[[i]] <- paste("g",c(s1:s2), sep="")
	s1 <- s2+1
	s2 <- s2+50
	print(s1)
	print(s2)	
}
gg <- GSEAx(m, ref=c(rep(1,20), rep(2,20)), pathways)
xx <- sscmp(m, pathways, ref=c(rep(1:15)), sample.size=c(5,10,15), method="GSEA", fdr="BH")
tmp <- powerpath(xx, pathname=paste("p",c(1:20),sep=""))
kk <- detectioncall(m, ref=c(rep(1,20), rep(2,20)), pathways=pathways)


