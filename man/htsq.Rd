\name{htsq}
\alias{htsq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
htsq(data, cls1, cls2, shrinkage = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{cls1}{
%%     ~~Describe \code{cls1} here~~
}
  \item{cls2}{
%%     ~~Describe \code{cls2} here~~
}
  \item{shrinkage}{
%%     ~~Describe \code{shrinkage} here~~
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
function (data, cls1, cls2, shrinkage = F) 
{
    meanS1 <- apply(data[, cls1], 1, mean)
    meanS2 <- apply(data[, cls2], 1, mean)
    covS1 <- covS2 <- NULL
    if ((nrow(data) - 1) > (ncol(data))) {
        shrinkage = T
    }
    if (shrinkage) {
        covS1 <- cov.shrink(t(data[, cls1]), verbose = F)
        covS2 <- cov.shrink(t(data[, cls2]), verbose = F)
    }
    else {
        covS1 <- cov(t(data[, cls1]))
        covS2 <- cov(t(data[, cls2]))
    }
    n1 = ncol(data[, cls1])
    n2 = ncol(data[, cls2])
    p <- nrow(data)
    covS = ((n1 - 1) * covS1 + (n2 - 1) * covS2)/(n1 + n2 - 2)
    covSinv <- solve(covS)
    htsqx <- (t(meanS1 - meanS2) \%*\% covSinv \%*\% (meanS1 - meanS2)) * 
        (n1 * n2/(n1 + n2))
    list(statistic = htsqx, df1 = p, df2 = n1 + n2 - p - 1, n1 = n1, 
        n2 = n2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
