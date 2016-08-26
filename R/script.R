cl <- 10;set.seed(3)
for(smp in 1:10){
print(cl)
rc <- 10;pth=1
pvlls <- list();pvgls <- list()
for(ind in 1:20){
print(rc);r=10;
pvl <- pvgl <- c()
for(i in 1:25){
bg <- c()
m <- matrix(rnorm(rc*cl),rc,cl)
m[1:(rc), ((cl/2)+1):cl] <- m[1:(rc), ((cl/2)+1):cl] + abs(rnorm(rc,0,.12))
pv <- pvg <- c()
bg <- c();r=10
for(j in 1:25){
	pt <- paste("g", c(1:rc), sep="")
	rownames(m) <- pt
	pt <- list(pt)
	names(pt) <- "pt1"
	colnames(m) <- paste("s",1:ncol(m),sep="")
	bg <- rbind(bg, matrix(rnorm(r*cl),r,cl))
	rownames(bg) <- paste("gx", 1:nrow(bg),sep="")
	mn <- rbind(m,bg)
	#gg <- gage(mn, gsets = pt, ref=c(1:(cl/2)),samp=c(((cl/2)+1):cl), compare= "unpaired")
	gg <- GSEAx(mn, ref=c(rep(0,cl/2),rep(1,cl/2)),pt,adjust.FDR.q.val=FALSE, nperm = 200)
	pv <- c(pv, as.numeric(as.vector(gg$report1[,"NOM p-val"])))
	pvg <- c(pvg, as.numeric(as.vector(gg$report1[,"glob.p.val"])))
	#pv <- c(pv, as.numeric(as.vector(gg$greater[,"p.val"])))
	r <- r+20
}
	pvl <- rbind(pvl, pv)
	pvgl <- rbind(pvgl, pvg)
} 
	pvlls[[ind]] <- pvl
	pvgls[[ind]] <- pvgl
	rc <- rc + 10
}
save(pvlls, file=paste("gsea",cl,"_",rc,sep=""))
cl <- cl+10
}

r=10;cl <- 30;rc <- 50
pvl <- pvgl <- c()
set.seed(3)
for(i in 1:25){

bg <- c()
m <- matrix(rnorm(600),rc,cl)
m[1:(rc), ((cl/2)+1):cl] <- m[1:(rc), ((cl/2)+1):cl] + abs(rnorm(rc,0,.12))
pv <- pvg <- c()
bg <- c()
r=10
for(j in 1:50){
	bg <- c()	
	pt <- paste("g", c(1:rc), sep="")
	rownames(m) <- pt
	pt <- list(pt)
	names(pt) <- "pt1"
	colnames(m) <- paste("s",1:ncol(m),sep="")
	bg <- rbind(bg, matrix(rnorm(r*cl),r,cl))
	rownames(bg) <- paste("gx", 1:nrow(bg),sep="")
	mn <- rbind(m,bg)
	#gg <- gage(mn, gsets = pt, ref=c(1:(cl/2)),samp=c(((cl/2)+1):cl), compare= "unpaired")
	gg <- GSEAx(mn, ref=c(rep(0,cl/2),rep(1,cl/2)),pt,adjust.FDR.q.val=FALSE, nperm = 200)
	pv <- c(pv, as.numeric(as.vector(gg$report1[,"NOM p-val"])))
	pvg <- c(pvg, as.numeric(as.vector(gg$report1[,"glob.p.val"])))
	#pv <- c(pv, as.numeric(as.vector(gg$greater[,"p.val"])))
	r <- r+30
}
	print(i)
	print(nrow(bg))
	pvl <- rbind(pvl, pv)
	pvgl <- rbind(pvgl, pvg)
} 

pvl1 <- cbind(pvl1, pvl)
boxplot(pvl1[ ,c(1:50,401,402,403)])
