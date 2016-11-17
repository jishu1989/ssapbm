\name{squared.t.test}
\alias{squared.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
squared.t.test(data, cls1, cls2, pathways = NULL)
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
  \item{pathways}{
%%     ~~Describe \code{pathways} here~~
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
function (data, cls1, cls2, pathways = NULL) 
{
    tval <- ptval <- c()
    tval <- ttest(data, cls1, cls2)
    for (i in 1:length(pathways)) {
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
