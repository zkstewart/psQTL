#!/usr/bin/env Rscript

loadRData <- function(fileName){
    load(fileName)
    mget(ls()[ls() != "fileName"])
}

if (!requireNamespace("argparser", quietly=TRUE)) {
  install.packages("argparser")
}
library(argparser)

# Establish parser
p <- arg_parser("Run sPLS-DA on features selected from 2 analyses to identify the most important QTLs")

# Add arguments
p <- add_argument(p, "call", help="Call sPLS-DA Rdata file", type = "character")
p <- add_argument(p, "depth", help="Depth sPLS-DA Rdata file", type = "character")
p <- add_argument(p, "o", help="Output file for selected variants", type = "character")
p <- add_argument(p, "--threads", help="Threads to use", type = "numeric", default = 1)
p <- add_argument(p, "--nrepeat", help="Number of repeats for stability analysis", type = "numeric", default = 10)
p <- add_argument(p, "--maxiters", help="Maximum number of iterations when tuning sPLS-DA", type = "numeric", default = 1000)

# Parse the command line arguments
args <- parse_args(p)

# Load in remaining libraries
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly=TRUE))
    BiocManager::install("mixOmics")
library(mixOmics)

if (!requireNamespace("BiocParallel", quietly=TRUE))
    BiocManager::install("BiocParallel")
library(BiocParallel)

# Set up parallel processing
if (args$threads == 1) {
  BPPARAM <- BiocParallel::SerialParam()
} else if (Sys.info()['sysname'] == "Windows") {
  BPPARAM <- BiocParallel::SnowParam(workers = args$threads)
} else {
  BPPARAM <- BiocParallel::MulticoreParam(workers = args$threads)
}

# Load each Rdata file to its own variable
callRdata <- loadRData(args$call)
depthRdata <- loadRData(args$depth)

# Extract X variables
call.X <- as.data.frame(callRdata$selected.X[,rownames(callRdata$feature.details.table)])
depth.X <- as.data.frame(depthRdata$selected.X[,rownames(depthRdata$feature.details.table)])

# Order datasets equivalently
call.X$Y <- callRdata$Y
depth.X$Y <- depthRdata$Y
depth.X <- depth.X[match(rownames(call.X), rownames(depth.X)),]

# Detect incompatible datasets
if (! all(rownames(call.X) == rownames(depth.X))){
  stop("X sample names do not match across datasets.")
}

if (! identical(call.X$Y, depth.X$Y)) {
  stop("Y variables do not match across datasets.")
}

# Extract Y variable and clean up X dataframes
Y <- call.X$Y
call.X <- call.X[, ! colnames(call.X) == "Y"]
depth.X <- depth.X[, ! colnames(depth.X) == "Y"]

# Set up list of X variables
splsda.X <- list(
  call = call.X,
  depth = depth.X
)

# Set design for multiblock sPLS-DA
design = matrix(0, ncol = 3, nrow = 3, 
                dimnames = list(c("call", "depth", "Y"), 
                                c("call", "depth", "Y")))
design[1,3] <- 1
design[3,1] <- 1
design[2,3] <- 1
design[3,2] <- 1

# Tune to obtain components
list.keepX = list(
  call = if (ncol(call.X) < 10) {
    1:ncol(call.X) } else if (ncol(call.X) < 20) {
      c(1:9, seq(10, ncol(call.X), 2)) } else {
        c(1:9,  seq(10, ncol(call.X)/2, 5))
      },
  depth = if (ncol(depth.X) < 10) {
    1:ncol(depth.X) } else if (ncol(depth.X) < 20) {
      c(1:9, seq(10, ncol(depth.X), 2)) } else {
        c(1:9,  seq(10, ncol(depth.X)/2, 5))
      }
)
tune.mbplsda <- tune.block.splsda(splsda.X, Y, test.keepX = list.keepX,
                                  ncomp = 1, folds = 2, design = design,
                                  scale = FALSE,
                                  nrepeat = args$nrepeat, max.iter = args$maxiters,
                                  BPPARAM = BPPARAM)

# Pick optimal keepX values
#select.keepX <- tune.mbplsda$choice.keepX # avoid default as it is not a sophisticated method
errors = as.data.frame(tune.mbplsda$error.rate)

balanceCall = list.keepX$call[ceiling(length(list.keepX$call) / 2)]
balanceDepth = list.keepX$depth[ceiling(length(list.keepX$depth) / 2)]

callEmphasis <- errors[paste0(max(list.keepX$call), "_1"),]
balance <- errors[paste0(balanceCall, "_", balanceDepth),]
depthEmphasis <- errors[paste0("1_", max(list.keepX$depth)),]

if (callEmphasis == min(callEmphasis, balance, depthEmphasis))
{
  callIndex <- length(list.keepX$call)
  depthIndex <- 1
  previousBest <- callEmphasis
  while ((depthIndex < length(list.keepX$depth) & depthIndex > 1) &
         (callIndex < length(list.keepX$call) & callIndex > 1))
  {
    callValue <- list.keepX$call[callIndex]
    depthValue <- list.keepX$depth[depthIndex]
    
    # Generate new values
    newDepthValue <- list.keepX$depth[depthIndex + 1]
    newDepthAlone <- errors[paste0(callValue, "_", newDepthValue),]
    
    newCallValue <- list.keepX$call[callIndex - 1]
    newCallAlone <- errors[paste0(newCallValue, "_", depthValue),]
    
    newBoth <- errors[paste0(newCallValue, "_", newDepthValue),]
    
    if (newBoth <= previousBest) {
      callIndex <- callIndex - 1
      depthIndex <- depthIndex + 1
      previousBest <- newBoth
      next
    } else if (newCallAlone <= previousBest) {
      callIndex <- callIndex - 1
      previousBest <- newCallAlone
      next
    } else if (newDepthAlone < previousBest) {
      depthIndex <- depthIndex + 
      previousBest <- newDepthAlone
      next
    }
    
    # Leapfrog down if we haven't find any incremental changes
    leapFrogged <- FALSE
    for (i in 1:(callIndex-1)) # callIndex must be >1
    {
      newCallDownValue <- list.keepX$call[i]
      newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
      if (newCallDown <= previousBest) {
        callIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in 1:(depthIndex-1))
    {
      newDepthDownValue <- list.keepX$depth[i]
      newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
      if (newDepthDown <= previousBest) {
        depthIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    
    if (leapFrogged == TRUE) {
      next
    }
    
    # Break if we didn't find any improvement
    break
  }
} else if (balance == min(callEmphasis, balance, depthEmphasis))
{
  callIndex = ceiling(length(list.keepX$call) / 2)
  depthIndex = ceiling(length(list.keepX$depth) / 2)
  previousBest <- balance
  
  while ((depthIndex < length(list.keepX$depth) & depthIndex > 1) &
         (callIndex < length(list.keepX$call) & callIndex > 1))
  {
    callValue <- list.keepX$call[callIndex]
    depthValue <- list.keepX$depth[depthIndex]
    
    # Generate down values
    newDepthDownValue <- list.keepX$depth[depthIndex - 1]
    newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
    
    newCallDownValue <- list.keepX$call[callIndex - 1]
    newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
    
    newBothDown <- errors[paste0(newCallDownValue, "_", newDepthDownValue),]
    
    # Generate up values
    newDepthUpValue <- list.keepX$depth[depthIndex + 1]
    newDepthUp <- errors[paste0(callValue, "_", newDepthUpValue),]
    
    newCallUpValue <- list.keepX$call[callIndex + 1]
    newCallUp <- errors[paste0(newCallUpValue, "_", depthValue),]
    
    newBothUp <- errors[paste0(newCallUpValue, "_", newDepthUpValue),]
    
    # Go down if we find something better OR equivalent
    if (newBothDown <= previousBest) {
      depthIndex <- depthIndex - 1
      callIndex <- callIndex - 1
      previousBest <- newBothDown
      next
    } else if (newCallDown <= previousBest) {
      callIndex <- callIndex - 1
      previousBest <- newCallDown
      next
    } else if (newDepthDown <= previousBest) {
      depthIndex <- depthIndex - 1
      previousBest <- newDepthDown
      next
    }
    
    # Go up if we find something better
    if (newBothUp < previousBest) {
      depthIndex <- depthIndex + 1
      callIndex <- callIndex + 1
      previousBest <- newBothUp
      next
    } else if (newCallUp < previousBest) {
      callIndex <- callIndex + 1
      previousBest <- newCallUp
      next
    } else if (newDepthUp < previousBest) {
      depthIndex <- depthIndex + 1
      previousBest <- newDepthUp
      next
    }
    
    # Leapfrog down if we haven't find any incremental changes
    leapFrogged <- FALSE
    for (i in 1:(callIndex-1)) # callIndex must be >1
    {
      newCallDownValue <- list.keepX$call[i]
      newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
      if (newCallDown <= previousBest) {
        callIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in 1:(depthIndex-1))
    {
      newDepthDownValue <- list.keepX$depth[i]
      newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
      if (newDepthDown <= previousBest) {
        depthIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    
    if (leapFrogged == TRUE) {
      next
    }
    
    # Break if we didn't find any improvement
    break
  }
} else {
  callIndex = 1
  depthIndex = length(list.keepX$depth)
  previousBest <- depthEmphasis
  while ((depthIndex < length(list.keepX$depth) & depthIndex > 1) &
         (callIndex < length(list.keepX$call) & callIndex > 1))
  {
    callValue <- list.keepX$call[callIndex]
    depthValue <- list.keepX$depth[depthIndex]
    
    # Generate new values
    newDepthValue <- list.keepX$depth[depthIndex - 1]
    newDepthAlone <- errors[paste0(callValue, "_", newDepthValue),]
    
    newCallValue <- list.keepX$call[callIndex + 1]
    newCallAlone <- errors[paste0(newCallValue, "_", depthValue),]
    
    newBoth <- errors[paste0(newCallValue, "_", newDepthValue),]
    
    if (newBoth < previousBest) {
      depthIndex <- depthIndex - 1
      callIndex <- callIndex + 1
      previousBest <- newBoth
      next
    } else if (newCallAlone < previousBest) {
      callIndex <- callIndex + 1
      previousBest <- newCallAlone
      next
    } else if (newDepthAlone < previousBest) {
      depthIndex <- depthIndex - 1
      previousBest <- newDepthAlone
      next
    }
    
    # Leapfrog down if we haven't find any incremental changes
    leapFrogged <- FALSE
    for (i in 1:(depthIndex-1))
    {
      newDepthDownValue <- list.keepX$depth[i]
      newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
      if (newDepthDown <= previousBest) {
        depthIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in 1:(callIndex-1))
    {
      newCallDownValue <- list.keepX$call[i]
      newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
      if (newCallDown <= previousBest) {
        callIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    
    if (leapFrogged == TRUE) {
      next
    }
    
    # Break if we didn't find any improvement
    break
  }
}
callValue <- list.keepX$call[callIndex]
depthValue <- list.keepX$depth[depthIndex]
select.keepX <- list(
    call = callValue,
    depth = depthValue
)

# Run sPLS-DA to select multi-QTLs and identify most important features
final.mbsplsda <- block.splsda(splsda.X, Y, keepX = select.keepX,
                               ncomp = 1, design = design,
                               scale = FALSE,
                               max.iter = args$maxiters)
# perf.final.mbsplsda <- perf(final.mbsplsda,
#                             folds = 2, validation = "Mfold",
#                             nrepeat = args$nrepeat,
#                             BPPARAM = BPPARAM)
if (ncol(final.mbsplsda$X) <= 5) {
    perf.final.mbsplsda <- perf(final.mbsplsda,
                                validation = "loo", ## TBD: what happens with stability?
                                BPPARAM = BPPARAM)
} else {
    perf.final.mbsplsda <- perf(final.mbsplsda,
                                folds = 2, validation = "Mfold",
                                nrepeat = args$nrepeat,
                                BPPARAM = BPPARAM)
}

# Tabulate stability values
select.name.call <- selectVar(final.mbsplsda, comp = 1)$call$name
select.name.depth <- selectVar(final.mbsplsda, comp = 1)$depth$name
select.name <- c(select.name.call, select.name.depth)

stability.table.call <- data.frame("features" = select.name.call)
stability.table.depth <- data.frame("features" = select.name.depth)
for (i in 1:args$nrepeat)
{
  stability.rep <- perf.final.mbsplsda$features$stable[[paste0("nrep", i)]]
  stability.call <- stability.rep$call$comp1
  stability.depth <- stability.rep$depth$comp1
  
  stability.table.call[[paste0("nrep", i)]] <- stability.call[match(stability.table.call$features, names(stability.call))]
  stability.table.depth[[paste0("nrep", i)]] <- stability.depth[match(stability.table.depth$features, names(stability.depth))]
}
stability.table.call$Freq <- rowMeans(stability.table.call[,seq(2, args$nrepeat+1)], na.rm=TRUE)
stability.table.call$type <- "call"
stability.table.depth$Freq <- rowMeans(stability.table.depth[,seq(2, args$nrepeat+1)], na.rm=TRUE)
stability.table.depth$type <- "depth"

stability.table <- rbind(stability.table.call, stability.table.depth)
rownames(stability.table) <- stability.table$features
stability.table <- stability.table[,c("Freq", "type"),drop=FALSE]

stability.table <- stability.table[order(stability.table$Freq, decreasing = TRUE),,drop=FALSE]

# Obtain loadings
mbsplsda.loadings.call <- as.data.frame(final.mbsplsda$loadings$call[,1,drop=FALSE])
mbsplsda.loadings.call <- mbsplsda.loadings.call[abs(rowSums(mbsplsda.loadings.call)) > 0,,drop=FALSE]
mbsplsda.loadings.call$type <- "call"

mbsplsda.loadings.depth <- as.data.frame(final.mbsplsda$loadings$depth[,1,drop=FALSE])
mbsplsda.loadings.depth <- mbsplsda.loadings.depth[abs(rowSums(mbsplsda.loadings.depth)) > 0,,drop=FALSE]
mbsplsda.loadings.depth$type <- "depth"

mbsplsda.loadings <- rbind(mbsplsda.loadings.call, mbsplsda.loadings.depth)

# Reformat loading values for ease of interpretation
mbsplsda.loadings$direction <- ifelse(mbsplsda.loadings$comp1 > 0, "right", "left")
mbsplsda.loadings$comp1 <- abs(mbsplsda.loadings$comp1)

# Match loadings order to stability values
mbsplsda.loadings <- mbsplsda.loadings[match(rownames(stability.table), rownames(mbsplsda.loadings)),,drop=FALSE]

# Join stability and loading values
feature.details.table <- cbind(stability.table, mbsplsda.loadings)
feature.details.table[c("chrom", "pos")] <- do.call(rbind, strsplit(rownames(feature.details.table), "_(?=[^_]+$)", perl=TRUE))
feature.details.table <- feature.details.table[,c("chrom", "pos", "type", "Freq", "comp1", "direction")]
colnames(feature.details.table) <- c("chrom", "pos", "type", "stability", "abs_loading", "direction")

# Write selected variants to file
write.table(feature.details.table, file=args$o, sep="\t", row.names=FALSE, quote=FALSE)
