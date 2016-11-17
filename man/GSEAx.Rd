\name{GSEAx}
\alias{GSEAx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GSEAx(dataset, ref = NULL, gsets, reshuffling.type = "sample.labels", nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, topgs = 10, adjust.FDR.q.val = F, reverse.sign = F, preproc.type = 0, random.seed = NULL, perm.type = 0, fraction = 1, replace = F, gs.size.threshold.min = 5, gs.size.threshold.max = 500, use.fast.enrichment.routine = T, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{gsets}{
%%     ~~Describe \code{gsets} here~~
}
  \item{reshuffling.type}{
%%     ~~Describe \code{reshuffling.type} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{weighted.score.type}{
%%     ~~Describe \code{weighted.score.type} here~~
}
  \item{nom.p.val.threshold}{
%%     ~~Describe \code{nom.p.val.threshold} here~~
}
  \item{fwer.p.val.threshold}{
%%     ~~Describe \code{fwer.p.val.threshold} here~~
}
  \item{fdr.q.val.threshold}{
%%     ~~Describe \code{fdr.q.val.threshold} here~~
}
  \item{topgs}{
%%     ~~Describe \code{topgs} here~~
}
  \item{adjust.FDR.q.val}{
%%     ~~Describe \code{adjust.FDR.q.val} here~~
}
  \item{reverse.sign}{
%%     ~~Describe \code{reverse.sign} here~~
}
  \item{preproc.type}{
%%     ~~Describe \code{preproc.type} here~~
}
  \item{random.seed}{
%%     ~~Describe \code{random.seed} here~~
}
  \item{perm.type}{
%%     ~~Describe \code{perm.type} here~~
}
  \item{fraction}{
%%     ~~Describe \code{fraction} here~~
}
  \item{replace}{
%%     ~~Describe \code{replace} here~~
}
  \item{gs.size.threshold.min}{
%%     ~~Describe \code{gs.size.threshold.min} here~~
}
  \item{gs.size.threshold.max}{
%%     ~~Describe \code{gs.size.threshold.max} here~~
}
  \item{use.fast.enrichment.routine}{
%%     ~~Describe \code{use.fast.enrichment.routine} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (dataset, ref = NULL, gsets, reshuffling.type = "sample.labels", 
    nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, 
    fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, topgs = 10, 
    adjust.FDR.q.val = F, reverse.sign = F, preproc.type = 0, 
    random.seed = NULL, perm.type = 0, fraction = 1, replace = F, 
    gs.size.threshold.min = 5, gs.size.threshold.max = 500, use.fast.enrichment.routine = T, 
    verbose = FALSE) 
{
    if (.Platform$OS.type == "windows") {
        memory.limit(6e+09)
        memory.limit()
    }
    if (!is.null(random.seed)) {
        set.seed(seed = random.seed, kind = NULL)
    }
    adjust.param <- 0.5
    gc()
    time1 <- proc.time()
    gene.labels <- row.names(dataset)
    gene.ann = ""
    gs.ann = ""
    if (class(dataset) == "data.frame") {
        sample.names <- names(dataset)
    }
    else {
        sample.names = colnames(dataset)
    }
    A <- data.matrix(dataset)
    dim(A)
    cols <- length(A[1, ])
    rows <- length(A[, 1])
    CLS <- ref
    if (class(CLS) != "list") {
        if (length(CLS) == ncol(dataset)) {
            tmp1 <- names(table(CLS))
            tmp2 <- rep(0, ncol(dataset))
            tmp2[which(CLS == tmp1[1])] <- 1
            tmp2[which(CLS == tmp1[2])] <- 0
            CLS <- list(phen = tmp1, class.v = tmp2)
        }
    }
    class.labels <- CLS$class.v
    class.phen <- CLS$phen
    if (reverse.sign == T) {
        phen1 <- class.phen[2]
        phen2 <- class.phen[1]
    }
    else {
        phen1 <- class.phen[1]
        phen2 <- class.phen[2]
    }
    col.index <- order(class.labels, decreasing = F)
    class.labels <- class.labels[col.index]
    sample.names <- sample.names[col.index]
    for (j in 1:rows) {
        A[j, ] <- A[j, col.index]
    }
    names(A) <- sample.names
    temp <- gsets
    max.Ng <- length(temp)
    temp.size.G <- vector(length = max.Ng, mode = "numeric")
    for (i in 1:max.Ng) {
        temp.size.G[i] <- length(temp[[i]])
    }
    max.size.G <- max(temp.size.G)
    gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, 
        ncol = max.size.G)
    temp.names <- vector(length = max.Ng, mode = "character")
    temp.desc <- vector(length = max.Ng, mode = "character")
    gs.count <- 1
    for (i in 1:max.Ng) {
        gene.set.size <- length(temp[[i]])
        gs.line <- noquote(temp[[i]])
        gene.set.name <- names(temp)[i]
        gene.set.desc <- noquote(" ")
        gene.set.tags <- gs.line
        existing.set <- is.element(gene.set.tags, gene.labels)
        set.size <- length(existing.set[existing.set == T])
        if ((set.size < gs.size.threshold.min) || (set.size > 
            gs.size.threshold.max)) 
            next
        temp.size.G[gs.count] <- set.size
        gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", 
            max.size.G - temp.size.G[gs.count]))
        temp.names[gs.count] <- gene.set.name
        temp.desc[gs.count] <- gene.set.desc
        gs.count <- gs.count + 1
    }
    Ng <- gs.count - 1
    gs.names <- vector(length = Ng, mode = "character")
    gs.desc <- vector(length = Ng, mode = "character")
    size.G <- vector(length = Ng, mode = "numeric")
    gs.names <- temp.names[1:Ng]
    gs.desc <- temp.desc[1:Ng]
    size.G <- temp.size.G[1:Ng]
    N <- length(A[, 1])
    Ns <- length(A[1, ])
    all.gene.descs <- gene.labels
    all.gene.symbols <- gene.labels
    all.gs.descs <- gs.desc
    if (verbose) {
        print(c("Number of genes:", N))
        print(c("Number of Gene Sets:", Ng))
        print(c("Number of samples:", Ns))
        print(c("Original number of Gene Sets:", max.Ng))
        print(c("Maximum gene set size:", max.size.G))
    }
    all.gene.descs <- vector(length = N, mode = "character")
    all.gene.symbols <- vector(length = N, mode = "character")
    all.gs.descs <- vector(length = Ng, mode = "character")
    if (is.data.frame(gene.ann)) {
        temp <- gene.ann
        a.size <- length(temp[, 1])
        print(c("Number of gene annotation file entries:", a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gene.labels, accs)
        all.gene.descs <- as.character(temp[locs, "Gene.Title"])
        all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
        rm(temp)
    }
    else if (gene.ann == "") {
        for (i in 1:N) {
            all.gene.descs[i] <- gene.labels[i]
            all.gene.symbols[i] <- gene.labels[i]
        }
    }
    else {
        temp <- read.delim(gene.ann, header = T, sep = ",", comment.char = "", 
            as.is = T)
        a.size <- length(temp[, 1])
        print(c("Number of gene annotation file entries:", a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gene.labels, accs)
        all.gene.descs <- as.character(temp[locs, "Gene.Title"])
        all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
        rm(temp)
    }
    if (is.data.frame(gs.ann)) {
        temp <- gs.ann
        a.size <- length(temp[, 1])
        print(c("Number of gene set annotation file entries:", 
            a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gs.names, accs)
        all.gs.descs <- as.character(temp[locs, "SOURCE"])
        rm(temp)
    }
    else if (gs.ann == "") {
        for (i in 1:Ng) {
            all.gs.descs[i] <- gs.desc[i]
        }
    }
    else {
        temp <- read.delim(gs.ann, header = T, sep = "\t", comment.char = "", 
            as.is = T)
        a.size <- length(temp[, 1])
        print(c("Number of gene set annotation file entries:", 
            a.size))
        accs <- as.character(temp[, 1])
        locs <- match(gs.names, accs)
        all.gs.descs <- as.character(temp[locs, "SOURCE"])
        rm(temp)
    }
    Obs.indicator <- matrix(nrow = Ng, ncol = N)
    Obs.RES <- matrix(nrow = Ng, ncol = N)
    Obs.ES <- vector(length = Ng, mode = "numeric")
    Obs.arg.ES <- vector(length = Ng, mode = "numeric")
    Obs.ES.norm <- vector(length = Ng, mode = "numeric")
    time2 <- proc.time()
    obs.s2n <- vector(length = N, mode = "numeric")
    signal.strength <- vector(length = Ng, mode = "numeric")
    tag.frac <- vector(length = Ng, mode = "numeric")
    gene.frac <- vector(length = Ng, mode = "numeric")
    coherence.ratio <- vector(length = Ng, mode = "numeric")
    obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
    correl.matrix <- matrix(nrow = N, ncol = nperm)
    obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
    order.matrix <- matrix(nrow = N, ncol = nperm)
    obs.order.matrix <- matrix(nrow = N, ncol = nperm)
    nperm.per.call <- 100
    n.groups <- nperm\%/\%nperm.per.call
    n.rem <- nperm\%\%nperm.per.call
    n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
    n.ends <- cumsum(n.perms)
    n.starts <- n.ends - n.perms + 1
    if (n.rem == 0) {
        n.tot <- n.groups
    }
    else {
        n.tot <- n.groups + 1
    }
    for (nk in 1:n.tot) {
        call.nperm <- n.perms[nk]
        if (verbose) {
            print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", 
                n.starts[nk], "--", n.ends[nk], sep = " "))
        }
        O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, 
            permutation.type = perm.type, sigma.correction = "GeneCluster", 
            fraction = fraction, replace = replace, reverse.sign = reverse.sign)
        gc()
        order.matrix[, n.starts[nk]:n.ends[nk]] <- O$order.matrix
        obs.order.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
        correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
        obs.correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
        rm(O)
    }
    obs.s2n <- apply(obs.correl.matrix, 1, median)
    obs.index <- order(obs.s2n, decreasing = T)
    obs.s2n <- sort(obs.s2n, decreasing = T)
    obs.gene.labels <- gene.labels[obs.index]
    obs.gene.descs <- all.gene.descs[obs.index]
    obs.gene.symbols <- all.gene.symbols[obs.index]
    for (r in 1:nperm) {
        correl.matrix[, r] <- correl.matrix[order.matrix[, r], 
            r]
    }
    for (r in 1:nperm) {
        obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[, 
            r], r]
    }
    gene.list2 <- obs.index
    for (i in 1:Ng) {
        if (verbose) {
            print(paste("Computing observed enrichment for gene set:", 
                i, gs.names[i], sep = " "))
        }
        gene.set <- gs[i, gs[i, ] != "null"]
        gene.set2 <- vector(length = length(gene.set), mode = "numeric")
        gene.set2 <- match(gene.set, gene.labels)
        GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, 
            gene.set = gene.set2, weighted.score.type = weighted.score.type, 
            correl.vector = obs.s2n)
        Obs.ES[i] <- GSEA.results$ES
        Obs.arg.ES[i] <- GSEA.results$arg.ES
        Obs.RES[i, ] <- GSEA.results$RES
        Obs.indicator[i, ] <- GSEA.results$indicator
        if (Obs.ES[i] >= 0) {
            tag.frac[i] <- sum(Obs.indicator[i, 1:Obs.arg.ES[i]])/size.G[i]
            gene.frac[i] <- Obs.arg.ES[i]/N
        }
        else {
            tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
            gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
        }
        signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * 
            (N/(N - size.G[i]))
    }
    phi <- matrix(nrow = Ng, ncol = nperm)
    phi.norm <- matrix(nrow = Ng, ncol = nperm)
    obs.phi <- matrix(nrow = Ng, ncol = nperm)
    if (reshuffling.type == "sample.labels") {
        for (i in 1:Ng) {
            if (verbose) {
                print(paste("Computing random permutations' enrichment for gene set:", 
                  i, gs.names[i], sep = " "))
            }
            gene.set <- gs[i, gs[i, ] != "null"]
            gene.set2 <- vector(length = length(gene.set), mode = "numeric")
            gene.set2 <- match(gene.set, gene.labels)
            for (r in 1:nperm) {
                gene.list2 <- order.matrix[, r]
                if (use.fast.enrichment.routine == F) {
                  GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = correl.matrix[, r])
                }
                else {
                  GSEA.results <- GSEA.EnrichmentScore2(gene.list = gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = correl.matrix[, r])
                }
                phi[i, r] <- GSEA.results$ES
            }
            if (fraction < 1) {
                for (r in 1:nperm) {
                  obs.gene.list2 <- obs.order.matrix[, r]
                  if (use.fast.enrichment.routine == F) {
                    GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                      gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                      correl.vector = obs.correl.matrix[, r])
                  }
                  else {
                    GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
                      gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                      correl.vector = obs.correl.matrix[, r])
                  }
                  obs.phi[i, r] <- GSEA.results$ES
                }
            }
            else {
                obs.gene.list2 <- obs.order.matrix[, 1]
                if (use.fast.enrichment.routine == F) {
                  GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.correl.matrix[, r])
                }
                else {
                  GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.correl.matrix[, r])
                }
                obs.phi[i, 1] <- GSEA.results$ES
                for (r in 2:nperm) {
                  obs.phi[i, r] <- obs.phi[i, 1]
                }
            }
            gc()
        }
    }
    else if (reshuffling.type == "gene.labels") {
        for (i in 1:Ng) {
            gene.set <- gs[i, gs[i, ] != "null"]
            gene.set2 <- vector(length = length(gene.set), mode = "numeric")
            gene.set2 <- match(gene.set, gene.labels)
            for (r in 1:nperm) {
                reshuffled.gene.labels <- sample(1:rows)
                if (use.fast.enrichment.routine == F) {
                  GSEA.results <- GSEA.EnrichmentScore(gene.list = reshuffled.gene.labels, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.s2n)
                }
                else {
                  GSEA.results <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.s2n)
                }
                phi[i, r] <- GSEA.results$ES
            }
            if (fraction < 1) {
                for (r in 1:nperm) {
                  obs.gene.list2 <- obs.order.matrix[, r]
                  if (use.fast.enrichment.routine == F) {
                    GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                      gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                      correl.vector = obs.correl.matrix[, r])
                  }
                  else {
                    GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
                      gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                      correl.vector = obs.correl.matrix[, r])
                  }
                  obs.phi[i, r] <- GSEA.results$ES
                }
            }
            else {
                obs.gene.list2 <- obs.order.matrix[, 1]
                if (use.fast.enrichment.routine == F) {
                  GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.correl.matrix[, r])
                }
                else {
                  GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
                    gene.set = gene.set2, weighted.score.type = weighted.score.type, 
                    correl.vector = obs.correl.matrix[, r])
                }
                obs.phi[i, 1] <- GSEA.results$ES
                for (r in 2:nperm) {
                  obs.phi[i, r] <- obs.phi[i, 1]
                }
            }
            gc()
        }
    }
    if (verbose) {
        print("Computing nominal p-values...")
    }
    p.vals <- matrix(0, nrow = Ng, ncol = 2)
    for (i in 1:Ng) {
        pos.phi <- NULL
        neg.phi <- NULL
        for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                pos.phi <- c(pos.phi, phi[i, j])
            }
            else {
                neg.phi <- c(neg.phi, phi[i, j])
            }
        }
        ES.value <- Obs.ES[i]
        if (ES.value >= 0) {
            p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), 
                digits = 5)
        }
        else {
            p.vals[i, 1] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), 
                digits = 5)
        }
    }
    erf <- function(x) {
        2 * pnorm(sqrt(2) * x)
    }
    KS.mean <- function(N) {
        S <- 0
        for (k in -100:100) {
            if (k == 0) 
                next
            S <- S + 4 * (-1)^(k + 1) * (0.25 * exp(-2 * k * 
                k * N) - sqrt(2 * pi) * erf(sqrt(2 * N) * k)/(16 * 
                k * sqrt(N)))
        }
        return(abs(S))
    }
    if (verbose) {
        print("Computing rescaling normalization for each gene set null...")
    }
    for (i in 1:Ng) {
        pos.phi <- NULL
        neg.phi <- NULL
        for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                pos.phi <- c(pos.phi, phi[i, j])
            }
            else {
                neg.phi <- c(neg.phi, phi[i, j])
            }
        }
        pos.m <- mean(pos.phi)
        neg.m <- mean(abs(neg.phi))
        pos.phi <- pos.phi/pos.m
        neg.phi <- neg.phi/neg.m
        for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                phi.norm[i, j] <- phi[i, j]/pos.m
            }
            else {
                phi.norm[i, j] <- phi[i, j]/neg.m
            }
        }
        for (j in 1:nperm) {
            if (obs.phi[i, j] >= 0) {
                obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
            }
            else {
                obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
            }
        }
        if (Obs.ES[i] >= 0) {
            Obs.ES.norm[i] <- Obs.ES[i]/pos.m
        }
        else {
            Obs.ES.norm[i] <- Obs.ES[i]/neg.m
        }
    }
    if (verbose) {
        print("Computing FWER p-values...")
    }
    max.ES.vals.p <- NULL
    max.ES.vals.n <- NULL
    for (j in 1:nperm) {
        pos.phi <- NULL
        neg.phi <- NULL
        for (i in 1:Ng) {
            if (phi.norm[i, j] >= 0) {
                pos.phi <- c(pos.phi, phi.norm[i, j])
            }
            else {
                neg.phi <- c(neg.phi, phi.norm[i, j])
            }
        }
        if (length(pos.phi) > 0) {
            max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
        }
        if (length(neg.phi) > 0) {
            max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
        }
    }
    for (i in 1:Ng) {
        ES.value <- Obs.ES.norm[i]
        if (Obs.ES.norm[i] >= 0) {
            p.vals[i, 2] <- signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), 
                digits = 5)
        }
        else {
            p.vals[i, 2] <- signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), 
                digits = 5)
        }
    }
    if (verbose) {
        print("Computing FDR q-values...")
    }
    NES <- vector(length = Ng, mode = "numeric")
    phi.norm.mean <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.mean <- vector(length = Ng, mode = "numeric")
    phi.norm.median <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.median <- vector(length = Ng, mode = "numeric")
    phi.norm.mean <- vector(length = Ng, mode = "numeric")
    obs.phi.mean <- vector(length = Ng, mode = "numeric")
    FDR.mean <- vector(length = Ng, mode = "numeric")
    FDR.median <- vector(length = Ng, mode = "numeric")
    phi.norm.median.d <- vector(length = Ng, mode = "numeric")
    obs.phi.norm.median.d <- vector(length = Ng, mode = "numeric")
    Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
    Orig.index <- seq(1, Ng)
    Orig.index <- Orig.index[Obs.ES.index]
    Orig.index <- order(Orig.index, decreasing = F)
    Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
    gs.names.sorted <- gs.names[Obs.ES.index]
    for (k in 1:Ng) {
        NES[k] <- Obs.ES.norm.sorted[k]
        ES.value <- NES[k]
        count.col <- vector(length = nperm, mode = "numeric")
        obs.count.col <- vector(length = nperm, mode = "numeric")
        for (i in 1:nperm) {
            phi.vec <- phi.norm[, i]
            obs.phi.vec <- obs.phi.norm[, i]
            if (ES.value >= 0) {
                count.col.norm <- sum(phi.vec >= 0)
                obs.count.col.norm <- sum(obs.phi.vec >= 0)
                count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= 
                  ES.value)/count.col.norm, 0)
                obs.count.col[i] <- ifelse(obs.count.col.norm > 
                  0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 
                  0)
            }
            else {
                count.col.norm <- sum(phi.vec < 0)
                obs.count.col.norm <- sum(obs.phi.vec < 0)
                count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= 
                  ES.value)/count.col.norm, 0)
                obs.count.col[i] <- ifelse(obs.count.col.norm > 
                  0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 
                  0)
            }
        }
        phi.norm.mean[k] <- mean(count.col)
        obs.phi.norm.mean[k] <- mean(obs.count.col)
        phi.norm.median[k] <- median(count.col)
        obs.phi.norm.median[k] <- median(obs.count.col)
        FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 
            1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
        FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 
            1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
    }
    if (adjust.FDR.q.val == T) {
        pos.nes <- length(NES[NES >= 0])
        min.FDR.mean <- FDR.mean[pos.nes]
        min.FDR.median <- FDR.median[pos.nes]
        for (k in seq(pos.nes - 1, 1, -1)) {
            if (FDR.mean[k] < min.FDR.mean) {
                min.FDR.mean <- FDR.mean[k]
            }
            if (min.FDR.mean < FDR.mean[k]) {
                FDR.mean[k] <- min.FDR.mean
            }
        }
        neg.nes <- pos.nes + 1
        min.FDR.mean <- FDR.mean[neg.nes]
        min.FDR.median <- FDR.median[neg.nes]
        for (k in seq(neg.nes + 1, Ng)) {
            if (FDR.mean[k] < min.FDR.mean) {
                min.FDR.mean <- FDR.mean[k]
            }
            if (min.FDR.mean < FDR.mean[k]) {
                FDR.mean[k] <- min.FDR.mean
            }
        }
    }
    obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
    phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
    FDR.mean.sorted <- FDR.mean[Orig.index]
    FDR.median.sorted <- FDR.median[Orig.index]
    glob.p.vals <- vector(length = Ng, mode = "numeric")
    NULL.pass <- vector(length = nperm, mode = "numeric")
    OBS.pass <- vector(length = nperm, mode = "numeric")
    for (k in 1:Ng) {
        NES[k] <- Obs.ES.norm.sorted[k]
        if (NES[k] >= 0) {
            for (i in 1:nperm) {
                NULL.pos <- sum(phi.norm[, i] >= 0)
                NULL.pass[i] <- ifelse(NULL.pos > 0, sum(phi.norm[, 
                  i] >= NES[k])/NULL.pos, 0)
                OBS.pos <- sum(obs.phi.norm[, i] >= 0)
                OBS.pass[i] <- ifelse(OBS.pos > 0, sum(obs.phi.norm[, 
                  i] >= NES[k])/OBS.pos, 0)
            }
        }
        else {
            for (i in 1:nperm) {
                NULL.neg <- sum(phi.norm[, i] < 0)
                NULL.pass[i] <- ifelse(NULL.neg > 0, sum(phi.norm[, 
                  i] <= NES[k])/NULL.neg, 0)
                OBS.neg <- sum(obs.phi.norm[, i] < 0)
                OBS.pass[i] <- ifelse(OBS.neg > 0, sum(obs.phi.norm[, 
                  i] <= NES[k])/OBS.neg, 0)
            }
        }
        glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
    }
    glob.p.vals.sorted <- glob.p.vals[Orig.index]
    Obs.ES <- signif(Obs.ES, digits = 5)
    Obs.ES.norm <- signif(Obs.ES.norm, digits = 5)
    p.vals <- signif(p.vals, digits = 4)
    signal.strength <- signif(signal.strength, digits = 3)
    tag.frac <- signif(tag.frac, digits = 3)
    gene.frac <- signif(gene.frac, digits = 3)
    FDR.mean.sorted <- signif(FDR.mean.sorted, digits = 5)
    FDR.median.sorted <- signif(FDR.median.sorted, digits = 5)
    glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits = 5)
    report <- data.frame(cbind(gs.names, size.G, all.gs.descs, 
        Obs.ES, Obs.ES.norm, p.vals[, 1], FDR.mean.sorted, p.vals[, 
            2], tag.frac, gene.frac, signal.strength, FDR.median.sorted, 
        glob.p.vals.sorted))
    names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", 
        "FDR q-val", "FWER p-val", "Tag \\\%", "Gene \\\%", "Signal", 
        "FDR (median)", "glob.p.val")
    report2 <- report
    report.index2 <- order(Obs.ES.norm, decreasing = T)
    for (i in 1:Ng) {
        report2[i, ] <- report[report.index2[i], ]
    }
    report3 <- report
    report.index3 <- order(Obs.ES.norm, decreasing = F)
    for (i in 1:Ng) {
        report3[i, ] <- report[report.index3[i], ]
    }
    phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
    phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
    report.phen1 <- report2[1:phen1.rows, ]
    report.phen2 <- report3[1:phen2.rows, ]
    rownames(report.phen1) <- as.vector(report.phen1[, 1])
    rownames(report.phen2) <- as.vector(report.phen2[, 1])
    return(list(report1 = report.phen1, report2 = report.phen2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
