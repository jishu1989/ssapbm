\name{ssapbm}
\alias{ssapbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssapbm(data = NULL, pathways = NULL, ref = NULL, target = NULL, minp = 2, maxp = 9999, method = c("sumoftsq", "ht2", "GSEA", "GAGE", "GSA"), perm = T, sampling = c("sample.labels", "gene.labels"), B = 100, sample.size = NULL, steps = 10, fdr = NULL, thr = 0.05, dc = T, rep = F, cr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{pathways}{
%%     ~~Describe \code{pathways} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{minp}{
%%     ~~Describe \code{minp} here~~
}
  \item{maxp}{
%%     ~~Describe \code{maxp} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{perm}{
%%     ~~Describe \code{perm} here~~
}
  \item{sampling}{
%%     ~~Describe \code{sampling} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{sample.size}{
%%     ~~Describe \code{sample.size} here~~
}
  \item{steps}{
%%     ~~Describe \code{steps} here~~
}
  \item{fdr}{
%%     ~~Describe \code{fdr} here~~
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
}
  \item{dc}{
%%     ~~Describe \code{dc} here~~
}
  \item{rep}{
%%     ~~Describe \code{rep} here~~
}
  \item{cr}{
%%     ~~Describe \code{cr} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data = NULL, pathways = NULL, ref = NULL, target = NULL, 
    minp = 2, maxp = 9999, method = c("sumoftsq", "ht2", "GSEA", 
        "GAGE", "GSA"), perm = T, sampling = c("sample.labels", 
        "gene.labels"), B = 100, sample.size = NULL, steps = 10, 
    fdr = NULL, thr = 0.05, dc = T, rep = F, cr = NULL) 
{
    if (is.null(target)) {
        cls2 <- c(1:ncol(data))[-ref]
    }
    else {
        cls2 <- target
    }
    d <- .processdata(data, pathways)
    data <- d[[1]]
    pathways <- d[[2]]
    cls1 <- ref
    TP <- pvl <- fdr.tp <- NULL
    if (is.null(fdr)) {
        fdr <- TRUE
        fdr.method <- "BH"
    }
    else {
        fdr.method <- fdr
        fdr <- TRUE
    }
    print(fdr.method)
    if (method[1] == "sumoftsq") {
        TP <- p.squared.t.test(data, cls1, cls2, steps = B, sampling = sampling[1], 
            pathways = pathways)
        if (!is.null(sample.size)) {
            pvl <- subsampanalysis(data, cls1, cls2, pathways, 
                B = 100, steps, sample.size = sample.size, method = "sumoftsq", 
                sampling = sampling[1], rep = rep, rn = rownames(TP), 
                dc = dc, thr = thr)
        }
        if (fdr) {
            fdr.tp <- p.adjust(TP[, "pval"], method = fdr.method)
        }
    }
    if (method[1] == "ht2") {
        TP <- pval.htsq(data, cls1, cls2, pathways = pathways, 
            perm = perm, sampling = sampling[1], steps = B)
        TP <- cbind(TP[[1]], TP[[2]])
        colnames(TP) <- c("ht2", "pval")
        if (!is.null(sample.size)) {
            pvl <- subsampanalysis(data, cls1, cls2, pathways, 
                B = 100, steps, sample.size = sample.size, method = "sumoftsq", 
                sampling = sampling[1], rep = rep, rn = rownames(TP), 
                dc = dc, thr = thr)
        }
        if (fdr) {
            fdr.tp <- p.adjust(TP[, 2], method = fdr.method)
        }
    }
    if (method[1] == "GSEA") {
        clsn <- rep(0, length(cls1) + length(cls2))
        clsn[cls2] <- 1
        lp <- sapply(pathways, length)
        lp <- c(min(lp), max(lp))
        TP <- GSEAx(data, ref = clsn, pathways, reshuffling.type = sampling[1], 
            nperm = B, gs.size.threshold.min = lp[1], gs.size.threshold.max = lp[2])
        TP <- as.matrix(rbind(TP[[1]][, c("ES", "NOM p-val")], 
            TP[[2]][, c("ES", "NOM p-val")]))
        print(TP)
        tmp1 <- cbind(as.numeric(as.vector(TP[, 1])), as.numeric(as.vector(TP[, 
            2])))
        rownames(tmp1) <- rownames(TP)
        colnames(tmp1) <- colnames(TP)
        TP <- tmp1
        TP <- TP[names(pathways), ]
        colnames(TP) <- c("GS", "pval")
        if (!is.null(sample.size)) {
            pvl <- subsampanalysis(data, cls1, cls2, pathways, 
                B = 100, steps, sample.size = sample.size, method = "GSEA", 
                sampling = sampling[1], rep = rep, rn = rownames(TP), 
                dc = dc, thr = thr)
        }
        if (fdr) {
            fdr.tp <- p.adjust(TP[, "pval"], method = fdr.method)
        }
    }
    if (method[1] == "GAGE") {
        lp <- sapply(pathways, length)
        lp <- c(min(lp), max(lp))
        print(head(data))
        pvl1 <- gage(data, gsets = pathways, ref = cls1, samp = cls2, 
            compare = "unpaired", set.size = c(lp[1], lp[2]), 
            same.dir = FALSE)
        TP <- pvl1[[1]][, 2:4]
        colnames(TP) <- c("stat.mean", "pval", "fdr.tp")
        if (fdr) {
            fdr.tp <- p.adjust(TP[, "pval"], method = fdr.method)
        }
        if (!is.null(sample.size)) {
            pvl <- subsampanalysis(data, cls1, cls2, pathways, 
                B = 100, steps, sample.size = sample.size, method = "GAGE", 
                dc = dc, thr = thr)
        }
    }
    if (method[1] == "GSA") {
        clsn <- rep(1, length(cls1) + length(cls2))
        clsn[cls2] <- 2
        lp <- sapply(pathways, length)
        lp <- c(min(lp), max(lp))
        TP <- GSA(data, clsn, genenames = rownames(data), genesets = pathways, 
            resp.type = "Two class unpaired", nperms = B, minsize = lp[1], 
            maxsize = lp[2])
        TP <- cbind(TP$GSA.scores, TP$pvalues.lo, TP$pvalues.hi)
        rownames(TP) <- names(pathways)
        tx <- (apply(TP, 1, function(x) c(which.min(as.numeric(x[2:3])))))
        TP <- cbind(TP, tx)
        print(TP)
        TP <- t(apply(TP, 1, function(x) x[c(1, (x[4] + 1))]))
        colnames(TP) <- c("GSA_score", "pval")
        if (fdr) {
            fdr.tp <- p.adjust(TP[, "pval"], method = fdr.method)
        }
        if (!is.null(sample.size)) {
            pvl <- subsampanalysis(data, cls1, cls2, pathways, 
                B = 100, steps, sample.size = sample.size, method = "GSA", 
                dc = dc, thr = thr)
        }
    }
    dcoverall <- corpaths <- NULL
    if (!is.null(cr)) {
        crall <- avgcorpath(data, ref = c(cls1, cls2), pathways, 
            method = cr)
        crctrl <- avgcorpath(data, ref = c(cls1), pathways, method = cr)
        crtreat <- avgcorpath(data, ref = c(cls2), pathways, 
            method = cr)
        corpaths <- cbind(crall, crctrl, crtreat)
    }
    if (dc) {
        clsn <- rep(1, length(cls1) + length(cls2))
        clsn[cls2] <- 2
        dcoverall <- pcdetectioncall(data, ref = clsn, pathways = pathways, 
            fdr = fdr.method, thr = thr)
    }
    if (fdr) {
        qval <- fdr.tp
        sigpath <- names(which(qval <= thr))
        list(cbind(TP, qval), pvl, dcoverall = dcoverall, corpaths = corpaths, 
            thr = thr, pathways = names(pathways), sigpath = sigpath, 
            fdr = fdr.method)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
