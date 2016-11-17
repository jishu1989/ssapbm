\name{p.squared.t.test}
\alias{p.squared.t.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p.squared.t.test(data, cls1, cls2, steps = 100, sampling = "sample.labels", pathways = NULL)
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
  \item{steps}{
%%     ~~Describe \code{steps} here~~
}
  \item{sampling}{
%%     ~~Describe \code{sampling} here~~
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
function (data, cls1, cls2, steps = 100, sampling = "sample.labels", 
    pathways = NULL) 
{
    tv <- squared.t.test(data, cls1, cls2, pathways = pathways)
    tvec <- resampling.squared.t.test(data, cls1, cls2, steps = steps, 
        pathways = pathways, sampling = sampling)
    tvec <- cbind(tv, tvec)
    pval <- apply(tvec, 1, function(x) {
        sum(x[2:(steps + 1)] > x[1])/steps
    })
    cbind(tv, pval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
