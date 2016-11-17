\name{plot.ssapbm}
\alias{plot.ssapbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.ssapbm(xx, data, type = c("power", "fdr", "fpr", "tpr", "dc", "cor"), val = "pval", pathways = NULL, nc = 1, nr = 1, pos = "bottom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{
%%     ~~Describe \code{xx} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{val}{
%%     ~~Describe \code{val} here~~
}
  \item{pathways}{
%%     ~~Describe \code{pathways} here~~
}
  \item{nc}{
%%     ~~Describe \code{nc} here~~
}
  \item{nr}{
%%     ~~Describe \code{nr} here~~
}
  \item{pos}{
%%     ~~Describe \code{pos} here~~
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
function (xx, data, type = c("power", "fdr", "fpr", "tpr", "dc", 
    "cor"), val = "pval", pathways = NULL, nc = 1, nr = 1, pos = "bottom") 
{
    print(type)
    pth <- NULL
    if (length(pathways) == 0) {
        pth <- xx$pathways
    }
    else {
        pth <- pathways
    }
    type <- type[1]
    typex = c("power", "fdr", "fpr", "tpr", "dc", "cor")
    nm <- c("Power", "FDR", "FPR", "TPR", "Detection call", "Correlation")
    names(nm) <- typex
    if ((type != "dc") && (type != "cor")) {
        pwr <- powerpath(xx, type = val)
        dat <- pwr[[type[1]]]
        dat <- getdf(dat)
        colnames(dat) <- c("sample", "stat")
    }
    dftrl <- list()
    if (type == "dc") {
        dcv <- xx[[2]][[length(xx[[2]])]]
        for (i in 1:(length(dcv))) {
            dctmp <- dcv[[i]][pth, ]
            dftr <- data.frame(rep(rownames(dctmp), ncol(dctmp)), 
                as.vector(dctmp))
            colnames(dftr) <- c("Pathways", "DC")
            xtst <- cbind(xx$dcoverall[pth], pth)
            colnames(xtst) <- c("dcall", "Pathways")
            xtst <- data.frame(xtst, type = factor(rep(1, length(pth))))
            xtst$dcall <- as.numeric(as.vector(xtst$dcall))
            print(xtst)
            gdc <- ggplot(dftr, aes(x = Pathways, y = DC, fill = Pathways)) + 
                geom_boxplot() + xlab(paste("Pathways", "(sample: ", 
                names(dcv)[i], ")", sep = "")) + ylab("Detection Call") + 
                theme(legend.position = "none") + theme(axis.text.x = element_blank())
            gdc <- gdc + scale_x_discrete(limits = pth)
            gdc <- gdc + geom_point(data = xtst, aes(x = Pathways, 
                y = dcall, colour = factor(type, labels = "Overall detection call")), 
                shape = 5, size = 2, fill = "red") + ylim(min(c(xtst$dcall, 
                dftr$DC)), max(c(xtst$dcall, dftr$DC))) + theme(legend.title = element_blank())
            dftrl[[i]] <- gdc
        }
    }
    if (type == "cor") {
        sig <- rep(0, length(xx$pathways))
        crdat <- data.frame(xx$corpaths, sig)
        crdat[xx$sigpath, "sig"] <- 1
        crdat$sig <- as.factor(crdat$sig)
        p1 <- ggplot(crdat, aes(x = rank(crall), y = crall, color = factor(sig, 
            labels = c("Differentially expressed", "Not Differentially expressed")))) + 
            geom_point() + theme(legend.position = "none") + 
            labs(color = "", x = "Ordered pathways", y = "Average correlation of individual pathways")
        p2 <- ggplot(crdat, aes(x = crctrl, y = crtreat, color = factor(sig))) + 
            geom_point() + theme(legend.position = "none") + 
            geom_line(data = crdat, color = rgb(0, 0, 1, alpha = 0.4), 
                size = 1.5, aes(x = crctrl, y = crctrl)) + labs(color = "", 
            x = "Average correlation:\n Individual pathways (control data)", 
            y = "Average correlation: \n Individual pathways (treatment data)")
        crd <- getdf(list(ctrl = crdat[xx$sigpath, 2], crtreat = crdat[xx$sigpath, 
            3]))
        print(crd)
        p3 <- ggplot(crd, aes(x = sample, y = stat, fill = sample)) + 
            geom_boxplot() + scale_x_discrete(limits = c("ctrl", 
            "crtreat"), breaks = c("ctrl", "crtreat"), labels = c("Control", 
            "Treatment")) + theme(legend.position = "none")
        grid_arrange_shared_legend(list(p1, p2, p3), nrow = nr, 
            ncol = nc, position = pos)
    }
    if (type == "dc") {
        print(length(dftrl))
        grid_arrange_shared_legend(dftrl, nrow = nr, ncol = nc, 
            position = pos)
    }
    if ((type != "dc") && (type != "cor")) {
        ggplot(dat, aes(y = stat, x = factor(sort(as.numeric(sample))), 
            fill = sample)) + geom_boxplot() + xlab("Sample-size") + 
            ylab(nm[type[1]])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
