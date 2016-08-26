GSEAx <- function(dataset, ref=NULL, gsets, reshuffling.type = "sample.labels", 
nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, topgs = 10, adjust.FDR.q.val = F, reverse.sign = F, preproc.type = 0, random.seed = NULL, perm.type = 0, fraction = 1.0, replace = F, gs.size.threshold.min = 5, gs.size.threshold.max = 500 ,use.fast.enrichment.routine = T, verbose=FALSE) {

  if (.Platform$OS.type == "windows") {
      memory.limit(6000000000)
      memory.limit()
#      print(c("Start memory size=",  memory.size()))
  }

  # Read input data matrix

 if(!is.null(random.seed)){ 
 	set.seed(seed=random.seed, kind = NULL)
 }	
  adjust.param <- 0.5

  gc()

  time1 <- proc.time()

  gene.labels <- row.names(dataset)
  gene.ann=""
  gs.ann =""  
  if(class(dataset)=="data.frame"){sample.names <- names(dataset)}
  else{sample.names= colnames(dataset)} 
  A <- data.matrix(dataset)
  dim(A) 
  cols <- length(A[1,])
  rows <- length(A[,1])

  # Read input class vector
  CLS <- ref
  if(class(CLS)!="list"){
  	if(length(CLS)==ncol(dataset)){
	  tmp1 <- names(table(CLS))
	  tmp2 <- rep(0, ncol(dataset))
	  tmp2[which(CLS==tmp1[1])] <- 1
	  tmp2[which(CLS==tmp1[2])] <- 0
	  CLS <- list(phen=tmp1, class.v=tmp2)	 	 		
  		
  	}
  }
  class.labels <- CLS$class.v
  class.phen <- CLS$phen

  if (reverse.sign == T) {
     phen1 <- class.phen[2]
     phen2 <- class.phen[1]
  } else {
     phen1 <- class.phen[1]
     phen2 <- class.phen[2]
  }

  # sort samples according to phenotype
 
 col.index <- order(class.labels, decreasing=F)
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
      gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
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
          if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
          temp.size.G[gs.count] <- set.size
          gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
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

  	N <- length(A[,1])
  	Ns <- length(A[1,])
 
        all.gene.descs <- gene.labels
        all.gene.symbols <- gene.labels
        all.gs.descs <- gs.desc
	if(verbose){	
		print(c("Number of genes:", N))
  		print(c("Number of Gene Sets:", Ng))
  		print(c("Number of samples:", Ns))
  		print(c("Original number of Gene Sets:", max.Ng))
  		print(c("Maximum gene set size:", max.size.G))
	}
		
# Read gene and gene set annotations if gene annotation file was provided

  all.gene.descs <- vector(length = N, mode ="character") 
  all.gene.symbols <- vector(length = N, mode ="character") 
  all.gs.descs <- vector(length = Ng, mode ="character") 

  if (is.data.frame(gene.ann)) {
     temp <- gene.ann
     a.size <- length(temp[,1])
     print(c("Number of gene annotation file entries:", a.size))  
     accs <- as.character(temp[,1])
     locs <- match(gene.labels, accs)
     all.gene.descs <- as.character(temp[locs, "Gene.Title"])
     all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
     rm(temp)
  } else  if (gene.ann == "") {
     for (i in 1:N) {
        all.gene.descs[i] <- gene.labels[i]
        all.gene.symbols[i] <- gene.labels[i]
     }
  } else {
     temp <- read.delim(gene.ann, header=T, sep=",", comment.char="", as.is=T)
     a.size <- length(temp[,1])
     print(c("Number of gene annotation file entries:", a.size))  
     accs <- as.character(temp[,1])
     locs <- match(gene.labels, accs)
     all.gene.descs <- as.character(temp[locs, "Gene.Title"])
     all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
     rm(temp)
  }

  if (is.data.frame(gs.ann)) {
     temp <- gs.ann
     a.size <- length(temp[,1])
     print(c("Number of gene set annotation file entries:", a.size))  
     accs <- as.character(temp[,1])
     locs <- match(gs.names, accs)
     all.gs.descs <- as.character(temp[locs, "SOURCE"])
     rm(temp)
  } else if (gs.ann == "") {
     for (i in 1:Ng) {
        all.gs.descs[i] <- gs.desc[i]
     }
  } else {
     temp <- read.delim(gs.ann, header=T, sep="\t", comment.char="", as.is=T)
     a.size <- length(temp[,1])
     print(c("Number of gene set annotation file entries:", a.size))  
     accs <- as.character(temp[,1])
     locs <- match(gs.names, accs)
     all.gs.descs <- as.character(temp[locs, "SOURCE"])
     rm(temp)
  }

  
  Obs.indicator <- matrix(nrow= Ng, ncol=N)
  Obs.RES <- matrix(nrow= Ng, ncol=N)

  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")
  Obs.ES.norm <- vector(length = Ng, mode = "numeric")

  time2 <- proc.time()

  # GSEA methodology

  # Compute observed and random permutation gene rankings

  obs.s2n <- vector(length=N, mode="numeric")
  signal.strength <- vector(length=Ng, mode="numeric")
  tag.frac <- vector(length=Ng, mode="numeric")
  gene.frac <- vector(length=Ng, mode="numeric")
  coherence.ratio <- vector(length=Ng, mode="numeric")
  obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
  correl.matrix <- matrix(nrow = N, ncol = nperm)
  obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
  order.matrix <- matrix(nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(nrow = N, ncol = nperm)

   nperm.per.call <- 100
   n.groups <- nperm %/% nperm.per.call
   n.rem <- nperm %% nperm.per.call
   n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
   n.ends <- cumsum(n.perms)
   n.starts <- n.ends - n.perms + 1

   if (n.rem == 0) {
     n.tot <- n.groups
   } else {
     n.tot <- n.groups + 1
   }

 for (nk in 1:n.tot) {
   call.nperm <- n.perms[nk]
   if(verbose){	
   	print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep=" "))
   }	
   O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, permutation.type = perm.type, sigma.correction = "GeneCluster", fraction=fraction, replace=replace, reverse.sign = reverse.sign)
   gc()

   order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
   obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
   correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
   obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
    rm(O)
 }

  obs.s2n <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
  obs.index <- order(obs.s2n, decreasing=T)            
  obs.s2n   <- sort(obs.s2n, decreasing=T)            

  obs.gene.labels <- gene.labels[obs.index]       
  obs.gene.descs <- all.gene.descs[obs.index]       
  obs.gene.symbols <- all.gene.symbols[obs.index]       

  for (r in 1:nperm) {
      correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
  }
  for (r in 1:nperm) {
      obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
  }

  gene.list2 <- obs.index
  for (i in 1:Ng) {
	if(verbose){
       		print(paste("Computing observed enrichment for gene set:", i, gs.names[i], sep=" ")) 
       }
	gene.set <- gs[i,gs[i,] != "null"]
       gene.set2 <- vector(length=length(gene.set), mode = "numeric")
       gene.set2 <- match(gene.set, gene.labels)
       GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
       Obs.ES[i] <- GSEA.results$ES
       Obs.arg.ES[i] <- GSEA.results$arg.ES
       Obs.RES[i,] <- GSEA.results$RES
       Obs.indicator[i,] <- GSEA.results$indicator
       if (Obs.ES[i] >= 0) {  # compute signal strength
           tag.frac[i] <- sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
           gene.frac[i] <- Obs.arg.ES[i]/N
       } else {
           tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
           gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
       }
       signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * (N / (N - size.G[i]))
   }

# Compute enrichment for random permutations 

   phi <- matrix(nrow = Ng, ncol = nperm)
   phi.norm <- matrix(nrow = Ng, ncol = nperm)
   obs.phi <- matrix(nrow = Ng, ncol = nperm)

   if (reshuffling.type == "sample.labels") { # reshuffling phenotype labels
      for (i in 1:Ng) {
	if(verbose){
        print(paste("Computing random permutations' enrichment for gene set:", i, gs.names[i], sep=" ")) 
        }
	 gene.set <- gs[i,gs[i,] != "null"]
        gene.set2 <- vector(length=length(gene.set), mode = "numeric")
        gene.set2 <- match(gene.set, gene.labels)
        for (r in 1:nperm) {
            gene.list2 <- order.matrix[,r]
            if (use.fast.enrichment.routine == F) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
            }
            phi[i, r] <- GSEA.results$ES
        }
        if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
            for (r in 1:nperm) {
                obs.gene.list2 <- obs.order.matrix[,r]
                if (use.fast.enrichment.routine == F) {
                   GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
                } else {
                   GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
                }
                obs.phi[i, r] <- GSEA.results$ES
            }
        } else { # if no resampling then compute only one column (and fill the others with the same value)
             obs.gene.list2 <- obs.order.matrix[,1]
            if (use.fast.enrichment.routine == F) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
            }
            obs.phi[i, 1] <- GSEA.results$ES
            for (r in 2:nperm) {
               obs.phi[i, r] <- obs.phi[i, 1]
            }
        }
        gc()
     }

   } else if (reshuffling.type == "gene.labels") { # reshuffling gene labels
      for (i in 1:Ng) {
        gene.set <- gs[i,gs[i,] != "null"]
        gene.set2 <- vector(length=length(gene.set), mode = "numeric")
        gene.set2 <- match(gene.set, gene.labels)
        for (r in 1:nperm) {
            reshuffled.gene.labels <- sample(1:rows)
            if (use.fast.enrichment.routine == F) {
               GSEA.results <- GSEA.EnrichmentScore(gene.list=reshuffled.gene.labels, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n)   
            } else {
               GSEA.results <- GSEA.EnrichmentScore2(gene.list=reshuffled.gene.labels, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n)   
            }
            phi[i, r] <- GSEA.results$ES
        }
        if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
           for (r in 1:nperm) {
              obs.gene.list2 <- obs.order.matrix[,r]
              if (use.fast.enrichment.routine == F) {
                 GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
              } else {
                 GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
              }
              obs.phi[i, r] <- GSEA.results$ES
           }
        } else { # if no resampling then compute only one column (and fill the others with the same value)
           obs.gene.list2 <- obs.order.matrix[,1]
           if (use.fast.enrichment.routine == F) {
              GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
           } else {
              GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
           }
           obs.phi[i, 1] <- GSEA.results$ES
           for (r in 2:nperm) {
              obs.phi[i, r] <- obs.phi[i, 1]
           }
        }
        gc()
     }
   }

# Compute 3 types of p-values

# Find nominal p-values       
if(verbose){
print("Computing nominal p-values...")
}
p.vals <- matrix(0, nrow = Ng, ncol = 2)

   for (i in 1:Ng) {
      pos.phi <- NULL
      neg.phi <- NULL
      for (j in 1:nperm) {
         if (phi[i, j] >= 0) {
            pos.phi <- c(pos.phi, phi[i, j]) 
         } else {
            neg.phi <- c(neg.phi, phi[i, j]) 
         }
      }
      ES.value <- Obs.ES[i]
      if (ES.value >= 0) {
         p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=5)
      } else {
         p.vals[i, 1] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=5)
      }
   }

# Find effective size 

 erf <- function (x) 
 {
    2 * pnorm(sqrt(2) * x)
 }

 KS.mean <- function(N) { # KS mean as a function of set size N
       S <- 0
       for (k in -100:100) {
          if (k == 0) next
          S <- S + 4 * (-1)**(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) *  erf(sqrt(2 * N) * k)/(16 * k * sqrt(N)))
       }
      return(abs(S))
 }

# KS.mean.table <- vector(length=5000, mode="numeric")

# for (i in 1:5000) {
#    KS.mean.table[i] <- KS.mean(i)
# }

# KS.size <-  vector(length=Ng, mode="numeric")

# Rescaling normalization for each gene set null
if(verbose){
print("Computing rescaling normalization for each gene set null...")
}
   for (i in 1:Ng) {
         pos.phi <- NULL
         neg.phi <- NULL
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
               pos.phi <- c(pos.phi, phi[i, j]) 
            } else {
               neg.phi <- c(neg.phi, phi[i, j]) 
            }
         }
         pos.m <- mean(pos.phi)
         neg.m <- mean(abs(neg.phi))

#         if (Obs.ES[i] >= 0) {
#            KS.size[i] <- which.min(abs(KS.mean.table - pos.m))
#         } else {
#            KS.size[i] <- which.min(abs(KS.mean.table - neg.m))
#         }

         pos.phi <- pos.phi/pos.m
         neg.phi <- neg.phi/neg.m
         for (j in 1:nperm) {
            if (phi[i, j] >= 0) {
                phi.norm[i, j] <- phi[i, j]/pos.m
            } else {
                phi.norm[i, j] <- phi[i, j]/neg.m
            }
          }
          for (j in 1:nperm) {
             if (obs.phi[i, j] >= 0) {
                obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
             } else {
                obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
             }
          }
          if (Obs.ES[i] >= 0) {
             Obs.ES.norm[i] <- Obs.ES[i]/pos.m
          } else {
             Obs.ES.norm[i] <- Obs.ES[i]/neg.m
          }
   }


# Compute FWER p-vals
	if(verbose){
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
          } else {
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
          p.vals[i, 2] <- signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), digits=5)
       } else {
          p.vals[i, 2] <- signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), digits=5)
       }
     }

# Compute FDRs 
	if(verbose){
      		print("Computing FDR q-values...")
	}
      NES <- vector(length=Ng, mode="numeric")
      phi.norm.mean  <- vector(length=Ng, mode="numeric")
      obs.phi.norm.mean  <- vector(length=Ng, mode="numeric")
      phi.norm.median  <- vector(length=Ng, mode="numeric")
      obs.phi.norm.median  <- vector(length=Ng, mode="numeric")
      phi.norm.mean  <- vector(length=Ng, mode="numeric")
      obs.phi.mean  <- vector(length=Ng, mode="numeric")
      FDR.mean <- vector(length=Ng, mode="numeric")
      FDR.median <- vector(length=Ng, mode="numeric")
      phi.norm.median.d <- vector(length=Ng, mode="numeric")
      obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")

      Obs.ES.index <- order(Obs.ES.norm, decreasing=T)
      Orig.index <- seq(1, Ng)
      Orig.index <- Orig.index[Obs.ES.index]
      Orig.index <- order(Orig.index, decreasing=F)
      Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
      gs.names.sorted <- gs.names[Obs.ES.index]

      for (k in 1:Ng) {
         NES[k] <- Obs.ES.norm.sorted[k]
         ES.value <- NES[k]
         count.col <- vector(length=nperm, mode="numeric")
         obs.count.col <- vector(length=nperm, mode="numeric")
         for (i in 1:nperm) {
            phi.vec <- phi.norm[,i]
            obs.phi.vec <- obs.phi.norm[,i]
            if (ES.value >= 0) {
               count.col.norm <- sum(phi.vec >= 0)
               obs.count.col.norm <- sum(obs.phi.vec >= 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
            } else {
               count.col.norm <- sum(phi.vec < 0)
               obs.count.col.norm <- sum(obs.phi.vec < 0)
               count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
               obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
            }
         }
        phi.norm.mean[k] <- mean(count.col)
        obs.phi.norm.mean[k] <- mean(obs.count.col)
        phi.norm.median[k] <- median(count.col)
        obs.phi.norm.median[k] <- median(obs.count.col)
        FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
        FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
      }

# adjust q-values

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

#   Compute global statistic

    glob.p.vals <- vector(length=Ng, mode="numeric")
    NULL.pass <- vector(length=nperm, mode="numeric")
    OBS.pass <- vector(length=nperm, mode="numeric")

    for (k in 1:Ng) {
      NES[k] <- Obs.ES.norm.sorted[k]
      if (NES[k] >= 0) {
         for (i in 1:nperm) {
             NULL.pos <- sum(phi.norm[,i] >= 0)            
             NULL.pass[i] <- ifelse(NULL.pos > 0, sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
             OBS.pos <- sum(obs.phi.norm[,i] >= 0)
             OBS.pass[i] <- ifelse(OBS.pos > 0, sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
         }
      } else {
         for (i in 1:nperm) {
             NULL.neg <- sum(phi.norm[,i] < 0)
             NULL.pass[i] <- ifelse(NULL.neg > 0, sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
             OBS.neg <- sum(obs.phi.norm[,i] < 0)
             OBS.pass[i] <- ifelse(OBS.neg > 0, sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
         }
      }
      glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
    }
    glob.p.vals.sorted <- glob.p.vals[Orig.index]

# Produce results report


       Obs.ES <- signif(Obs.ES, digits=5)
       Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
       p.vals <- signif(p.vals, digits=4)
       signal.strength <- signif(signal.strength, digits=3)
       tag.frac <- signif(tag.frac, digits=3)
       gene.frac <- signif(gene.frac, digits=3)
       FDR.mean.sorted <- signif(FDR.mean.sorted, digits=5)
       FDR.median.sorted <-  signif(FDR.median.sorted, digits=5)
       glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits=5)

       report <- data.frame(cbind(gs.names, size.G, all.gs.descs, Obs.ES, Obs.ES.norm, p.vals[,1], FDR.mean.sorted, p.vals[,2], tag.frac, gene.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
       names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag \\%", "Gene \\%", "Signal", "FDR (median)", "glob.p.val")
#       print(report)
       report2 <- report
       report.index2 <- order(Obs.ES.norm, decreasing=T)
       for (i in 1:Ng) {
           report2[i,] <- report[report.index2[i],]
       }   
       report3 <- report
       report.index3 <- order(Obs.ES.norm, decreasing=F)
       for (i in 1:Ng) {
           report3[i,] <- report[report.index3[i],]
       }   
       phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
       phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
       report.phen1 <- report2[1:phen1.rows,]
       report.phen2 <- report3[1:phen2.rows,]
	rownames(report.phen1) <- as.vector(report.phen1[,1])
	rownames(report.phen2) <- as.vector(report.phen2[,1])
  return(list(report1 = report.phen1, report2 = report.phen2))

}  # end of definition of GSEA.analysis

# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.


# G S E A -- Gene Set Enrichment Analysis

# Auxiliary functions and definitions 

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F) { 

# This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap  
# subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation 
# in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced. 
# It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
# the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for 
# all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second 
# matrix for the null distribution will still have the values for the random permutations 
# (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
# It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain 
# smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity 
# checks before trusting the code.
#
# Inputs:
#   A: Matrix of gene expression values (rows are genes, columns are samples) 
#   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
#   gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of the expression matrix 
#   nperm: Number of random permutations/bootstraps to perform 
#   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
#   sigma.correction: Correction to the signal to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package) 
#   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
#   replace: Resampling mode (replacement or not replacement). For experts only (default: F) 
#   reverse.sign: Reverse direction of gene list (default = F)
#
# Outputs:
#   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios (rows are genes, columns are permutations or bootstrap subsamplings
#   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
#   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
#   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

     A <- A + 0.00000001

     N <- length(A[,1])
     Ns <- length(A[1,])

     subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
     class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)

     order.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
     s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
     obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)

     obs.gene.labels <- vector(length = N, mode="character")
     obs.gene.descs <- vector(length = N, mode="character")
     obs.gene.symbols <- vector(length = N, mode="character")

     M1 <- matrix(0, nrow = N, ncol = nperm)
     M2 <- matrix(0, nrow = N, ncol = nperm)
     S1 <- matrix(0, nrow = N, ncol = nperm)
     S2 <- matrix(0, nrow = N, ncol = nperm)

     gc()

     C <- split(class.labels, class.labels)
     class1.size <- length(C[[1]])
     class2.size <- length(C[[2]])
     class1.index <- seq(1, class1.size, 1)
     class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)

     for (r in 1:nperm) {
        class1.subset <- sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
        class2.subset <- sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
        class1.subset.size <- length(class1.subset)
        class2.subset.size <- length(class2.subset)
        subset.class1 <- rep(0, class1.size)
        for (i in 1:class1.size) {
            if (is.element(class1.index[i], class1.subset)) {
                subset.class1[i] <- 1
            }
        }
        subset.class2 <- rep(0, class2.size)
        for (i in 1:class2.size) {
            if (is.element(class2.index[i], class2.subset)) {
                subset.class2[i] <- 1
            }
        }
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        fraction.class1 <- class1.size/Ns
        fraction.class2 <- class2.size/Ns

        if (permutation.type == 0) { # random (unbalanced) permutation
           full.subset <- c(class1.subset, class2.subset)
           label1.subset <- sample(full.subset, size = Ns * fraction.class1)
           reshuffled.class.labels1[, r] <- rep(0, Ns)
           reshuffled.class.labels2[, r] <- rep(0, Ns)
           class.labels1[, r] <- rep(0, Ns)
           class.labels2[, r] <- rep(0, Ns)
           for (i in 1:Ns) {
               m1 <- sum(!is.na(match(label1.subset, i)))
               m2 <- sum(!is.na(match(full.subset, i)))
               reshuffled.class.labels1[i, r] <- m1
               reshuffled.class.labels2[i, r] <- m2 - m1
               if (i <= class1.size) {
                 class.labels1[i, r] <- m2
                 class.labels2[i, r] <- 0
               } else {
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
               }
           }

        } else if (permutation.type == 1) { # proportional (balanced) permutation

           class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
           class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
           reshuffled.class.labels1[, r] <- rep(0, Ns)
           reshuffled.class.labels2[, r] <- rep(0, Ns)
           class.labels1[, r] <- rep(0, Ns)
           class.labels2[, r] <- rep(0, Ns)
           for (i in 1:Ns) {
               if (i <= class1.size) {
                  m1 <- sum(!is.na(match(class1.label1.subset, i)))
                  m2 <- sum(!is.na(match(class1.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
               } else {
                  m1 <- sum(!is.na(match(class2.label1.subset, i)))
                  m2 <- sum(!is.na(match(class2.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
               }
           }
        }
    }

# compute S2N for the random permutation matrix
     
     P <- reshuffled.class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- reshuffled.class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

     if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
         S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
         S2 <- ifelse(S2 == 0, 0.2, S2)
         S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
         S1 <- ifelse(S1 == 0, 0.2, S1)
         gc()
     }

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     s2n.matrix <- M1/S1

     if (reverse.sign == T) {
        s2n.matrix <- - s2n.matrix
     }
     gc()

     for (r in 1:nperm) {
        order.matrix[, r] <- order(s2n.matrix[, r], decreasing=T)            
     }

# compute S2N for the "observed" permutation matrix

     P <- class.labels1 * subset.mask
     n1 <- sum(P[,1])         
     M1 <- A %*% P
     M1 <- M1/n1      
     gc()
     A2 <- A*A        
     S1 <- A2 %*% P   
     S1 <- S1/n1 - M1*M1    
     S1 <- sqrt(abs((n1/(n1-1)) * S1))   
     gc()
     P <- class.labels2 * subset.mask
     n2 <- sum(P[,1])           
     M2 <- A %*% P           
     M2 <- M2/n2          
     gc()
     A2 <- A*A           
     S2 <- A2 %*% P      
     S2 <- S2/n2 - M2*M2 
     S2 <- sqrt(abs((n2/(n2-1)) * S2))
     rm(P)
     rm(A2)
     gc()

     if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
         S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
         S2 <- ifelse(S2 == 0, 0.2, S2)
         S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
         S1 <- ifelse(S1 == 0, 0.2, S1)
         gc()
     } 

     M1 <- M1 - M2
     rm(M2)
     gc()
     S1 <- S1 + S2
     rm(S2)
     gc()

     obs.s2n.matrix <- M1/S1
     gc()

     if (reverse.sign == T) {
        obs.s2n.matrix <- - obs.s2n.matrix
     }

     for (r in 1:nperm) {
        obs.order.matrix[,r] <- order(obs.s2n.matrix[,r], decreasing=T)            
     }

     return(list(s2n.matrix = s2n.matrix, 
                 obs.s2n.matrix = obs.s2n.matrix, 
                 order.matrix = order.matrix,
                 obs.order.matrix = obs.order.matrix))
}

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
#
# Computes the weighted GSEA score of gene.set in gene.list. 
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 
   if (weighted.score.type == 0) {
      correl.vector <- rep(1, N)
   }
   alpha <- weighted.score.type
   correl.vector <- abs(correl.vector**alpha)
   sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
   norm.tag    <- 1.0/sum.correl.tag
   norm.no.tag <- 1.0/Nm
   RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
#      ES <- max.ES
      ES <- signif(max.ES, digits = 5)
      arg.ES <- which.max(RES)
   } else {
#      ES <- min.ES
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}


OLD.GSEA.EnrichmentScore <- function(gene.list, gene.set) {  
#
# Computes the original GSEA score from Mootha et al 2003 of gene.set in gene.list 
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
   no.tag.indicator <- 1 - tag.indicator 
   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 

   norm.tag    <- sqrt((N - Nh)/Nh)
   norm.no.tag <- sqrt(Nh/(N - Nh))

   RES <- cumsum(tag.indicator * norm.tag - no.tag.indicator * norm.no.tag)      
   max.ES <- max(RES)
   min.ES <- min(RES)
   if (max.ES > - min.ES) {
      ES <- signif(max.ES, digits=5)
      arg.ES <- which.max(RES)
   } else {
      ES <- signif(min.ES, digits=5)
      arg.ES <- which.min(RES)
   }
   return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
#
# Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
# GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
# This call is intended to be used to asses the enrichment of random permutations rather than the 
# observed one.
# The weighted score type is the exponent of the correlation 
# weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
# necessary to input the correlation vector with the values in the same order as in the gene list.
#
# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   N <- length(gene.list) 
   Nh <- length(gene.set) 
   Nm <-  N - Nh 

   loc.vector <- vector(length=N, mode="numeric")
   peak.res.vector <- vector(length=Nh, mode="numeric")
   valley.res.vector <- vector(length=Nh, mode="numeric")
   tag.correl.vector <- vector(length=Nh, mode="numeric")
   tag.diff.vector <- vector(length=Nh, mode="numeric")
   tag.loc.vector <- vector(length=Nh, mode="numeric")

   loc.vector[gene.list] <- seq(1, N)
   tag.loc.vector <- loc.vector[gene.set]

   tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

   if (weighted.score.type == 0) {
      tag.correl.vector <- rep(1, Nh)
   } else if (weighted.score.type == 1) {
      tag.correl.vector <- correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else if (weighted.score.type == 2) {
      tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
      tag.correl.vector <- abs(tag.correl.vector)
   } else {
      tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
      tag.correl.vector <- abs(tag.correl.vector)
   }

   norm.tag <- 1.0/sum(tag.correl.vector)
   tag.correl.vector <- tag.correl.vector * norm.tag
   norm.no.tag <- 1.0/Nm
   tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
   tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
   tag.diff.vector <- tag.diff.vector * norm.no.tag
   peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
   valley.res.vector <- peak.res.vector - tag.correl.vector
   max.ES <- max(peak.res.vector)
   min.ES <- min(valley.res.vector)
   ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)

   return(list(ES = ES))

}


GSEA.Res2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

GSEA.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
   ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T)
   ds <- ds[-1]
   return(ds)
}

GSEA.Gct2Frame2 <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
      content <- readLines(filename)
      content <- content[-1]
      content <- content[-1]
      col.names <- noquote(unlist(strsplit(content[1], "\t")))
      col.names <- col.names[c(-1, -2)]
      num.cols <- length(col.names)
      content <- content[-1]
      num.lines <- length(content)


      row.nam <- vector(length=num.lines, mode="character")
      row.des <- vector(length=num.lines, mode="character")
      m <- matrix(0, nrow=num.lines, ncol=num.cols)

      for (i in 1:num.lines) {
         line.list <- noquote(unlist(strsplit(content[i], "\t")))
         row.nam[i] <- noquote(line.list[1])
         row.des[i] <- noquote(line.list[2])
         line.list <- line.list[c(-1, -2)]
         for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
         }
      }
      ds <- data.frame(m)
      names(ds) <- col.names
      row.names(ds) <- row.nam
      return(ds)
}

GSEA.ReadClsFile <- function(file = "NULL") { 
#
# Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

      cls.cont <- readLines(file)
      num.lines <- length(cls.cont)
      class.list <- unlist(strsplit(cls.cont[[3]], " "))
      s <- length(class.list)
      t <- table(class.list)
      l <- length(t)
      phen <- vector(length=l, mode="character")
      phen.label <- vector(length=l, mode="numeric")
      class.v <- vector(length=s, mode="numeric")
      for (i in 1:l) {
         phen[i] <- noquote(names(t)[i])
         phen.label[i] <- i - 1
      }
      for (i in 1:s) {
         for (j in 1:l) {
             if (class.list[i] == phen[j]) {
                class.v[i] <- phen.label[j]
             }
         }
      }
      return(list(phen = phen, class.v = class.v))
}

GSEA.Threshold <- function(V, thres, ceil) { 
#
# Threshold and ceiling pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        V[V < thres] <- thres
        V[V > ceil] <- ceil
        return(V)
}

GSEA.VarFilter <- function(V, fold, delta, gene.names = "NULL") { 
#
# Variation filter pre-processing for gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        cols <- length(V[1,])
        rows <- length(V[,1])
        row.max <- apply(V, MARGIN=1, FUN=max)
               row.min <- apply(V, MARGIN=1, FUN=min)
        flag <- array(dim=rows)
        flag <- (row.max /row.min > fold) & (row.max - row.min > delta)
        size <- sum(flag)
        B <- matrix(0, nrow = size, ncol = cols)
        j <- 1
        if (gene.names == "NULL") {
           for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 j <- j + 1
               }
           }
        return(B)
        } else {
            new.list <- vector(mode = "character", length = size)
            for (i in 1:rows) {
              if (flag[i]) {
                 B[j,] <- V[i,]
                 new.list[j] <- gene.names[i]
                 j <- j + 1
              }
            }
        return(list(V = B, new.list = new.list))
        }
}

GSEA.NormalizeRows <- function(V) { 
#
# Stardardize rows of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        row.mean <- apply(V, MARGIN=1, FUN=mean)
               row.sd <- apply(V, MARGIN=1, FUN=sd)
        row.n <- length(V[,1])
        for (i in 1:row.n) {
             if (row.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[i,] <- (V[i,] - row.mean[i])/row.sd[i]
           }
        }
        return(V)
}

GSEA.NormalizeCols <- function(V) { 
#
# Stardardize columns of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        col.mean <- apply(V, MARGIN=2, FUN=mean)
               col.sd <- apply(V, MARGIN=2, FUN=sd)
        col.n <- length(V[1,])
        for (i in 1:col.n) {
             if (col.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
           }
        }
        return(V)
}
