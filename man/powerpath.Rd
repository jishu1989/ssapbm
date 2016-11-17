\name{powerpath}
\alias{powerpath}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
powerpath(x, type = "pval")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (x, type = "pval") 
{
    tp <- x[[1]]
    thr <- x$thr
    fdr <- x$fdr
    k <- x$sigpath
    if (length(k) == 0) {
        stop("no p-value is below threshold")
    }
    resampval <- length(x[[2]][-length(x[[2]])])
    rml <- list()
    tprx <- fprx <- fdrx <- list()
    tpx <- k
    tnx <- setdiff(rownames(tp), k)
    for (i in 1:resampval) {
        p <- x[[2]][[i]]
        p <- apply(p, 2, function(x) p.adjust(x, method = fdr))
        ptmp <- apply(p, 2, function(x) {
            m <- rep(0, length(x))
            m[which(x <= thr)] <- 1
            m
        })
        rownames(ptmp) <- rownames(p)
        rml[[i]] <- (rowMeans(ptmp[k, ]))
        fdrtmp <- fprtmp <- tprtmp <- c()
        for (j in 1:ncol(ptmp)) {
            n1 <- rownames(ptmp[which(ptmp[, j] == 1), ])
            tprtmp <- c(tprtmp, length(intersect(tpx, n1))/length(tpx))
            fprtmp <- c(fprtmp, length(intersect(tnx, n1))/length(tnx))
            fdrtmp <- c(fdrtmp, length(intersect(tnx, n1))/length(n1))
        }
        tprx[[i]] <- tprtmp
        fprx[[i]] <- fprtmp
        fdrx[[i]] <- fdrtmp
    }
    names(tprx) <- names(fprx) <- names(fdrx) <- names(rml) <- names(xx[[2]][1:resampval])
    list(power = rml, tpr = tprx, fpr = fprx, fdr = fdrx, sigpath = k)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
