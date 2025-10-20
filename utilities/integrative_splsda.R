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
p <- add_argument(p, "os", help="Output file for selected variants", type = "character")
p <- add_argument(p, "op", help="Output file for predicted sample groups", type = "character")
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

if (!requireNamespace("stringr", quietly=TRUE))
    install.packages("stringr")
library(stringr)

# Identify mixOmics version
mixomicsVersion = unlist(packageVersion("mixOmics"))

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

# Detect a lack of need to integrate
if (ncol(call.X) < 2) {
  stop("Integrative sPLS-DA is not possible or necessary since only one or zero 'call' variants were selected")
}
if (ncol(depth.X) < 2) {
  stop("Integrative sPLS-DA is not possible or necessary since only one or zero 'depth' CNVs were selected")
}

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
errors <- as.data.frame(tune.mbplsda$error.rate)

balanceCall <- list.keepX$call[ceiling(length(list.keepX$call) / 2)]
balanceDepth <- list.keepX$depth[ceiling(length(list.keepX$depth) / 2)]

callEmphasis <- errors[paste0(max(list.keepX$call), "_1"),]
balance <- errors[paste0(balanceCall, "_", balanceDepth),]
depthEmphasis <- errors[paste0("1_", max(list.keepX$depth)),]

if (callEmphasis == min(callEmphasis, balance, depthEmphasis))
{
  callIndex <- length(list.keepX$call)
  depthIndex <- 1
  previousBest <- callEmphasis
  while ((depthIndex < length(list.keepX$depth) & depthIndex >= 1) &
         (callIndex <= length(list.keepX$call) & callIndex > 1))
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
      depthIndex <- depthIndex + 1
      previousBest <- newDepthAlone
      next
    }
    
    # Leapfrog down if we haven't find any incremental changes
    leapFrogged <- FALSE
    for (i in (callIndex-1):1) # callIndex must be >1
    {
      if (i==0)
      {
        break
      }
      newCallDownValue <- list.keepX$call[i]
      newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
      if (newCallDown <= previousBest) {
        callIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in (depthIndex-1):1)
    {
      if (i==0)
      {
        break
      }
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
  callIndex <- ceiling(length(list.keepX$call) / 2)
  depthIndex <- ceiling(length(list.keepX$depth) / 2)
  previousBest <- balance
  
  while ((depthIndex <= length(list.keepX$depth) & depthIndex >= 1) &
         (callIndex <= length(list.keepX$call) & callIndex >= 1))
  {
    callValue <- list.keepX$call[callIndex]
    depthValue <- list.keepX$depth[depthIndex]
    
    # Generate down values
    newDepthDownValue <- list.keepX$depth[depthIndex - 1]
    newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
    if (is.na(newDepthDown)) { newDepthDown <- 0.5 }
    
    newCallDownValue <- list.keepX$call[callIndex - 1]
    newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
    if (is.na(newCallDown)) { newCallDown <- 0.5 }
    
    newBothDown <- errors[paste0(newCallDownValue, "_", newDepthDownValue),]
    if (is.na(newBothDown)) { newBothDown <- 0.5 }
    
    # Generate up values
    newDepthUpValue <- list.keepX$depth[depthIndex + 1]
    newDepthUp <- errors[paste0(callValue, "_", newDepthUpValue),]
    if (is.na(newDepthUp)) { newDepthUp <- 0.5 }
    
    newCallUpValue <- list.keepX$call[callIndex + 1]
    newCallUp <- errors[paste0(newCallUpValue, "_", depthValue),]
    if (is.na(newCallUp)) { newCallUp <- 0.5 }
    
    newBothUp <- errors[paste0(newCallUpValue, "_", newDepthUpValue),]
    if (is.na(newBothUp)) { newBothUp <- 0.5 }
    
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
    for (i in (callIndex-1):1) # callIndex must be >1
    {
      if (i==0)
      {
        break
      }
      newCallDownValue <- list.keepX$call[i]
      newCallDown <- errors[paste0(newCallDownValue, "_", depthValue),]
      if (newCallDown <= previousBest) {
        callIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in (depthIndex-1):1)
    {
      if (i==0)
      {
        break
      }
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
  callIndex <- 1
  depthIndex <- length(list.keepX$depth)
  previousBest <- depthEmphasis
  while ((depthIndex <= length(list.keepX$depth) & depthIndex > 1) &
         (callIndex < length(list.keepX$call) & callIndex >= 1))
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
    for (i in (depthIndex-1):1)
    {
      if (i==0)
      {
        break
      }
      newDepthDownValue <- list.keepX$depth[i]
      newDepthDown <- errors[paste0(callValue, "_", newDepthDownValue),]
      if (newDepthDown <= previousBest) {
        depthIndex <- i
        leapFrogged <- TRUE
        break
      }
    }
    for (i in (callIndex-1):1)
    {
      if (i==0)
      {
        break
      }
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
tryCatch(
    {
        if (mixomicsVersion[2] <= 30) {
        perf.final.mbsplsda <- perf(final.mbsplsda,
                                    validation = "loo")
        } else {
        perf.final.mbsplsda <- perf(final.mbsplsda,
                                    validation = "loo")
        }
    },
    error = function(e) {
        print(paste0("Warning: stability cannot be estimated due to error \"", conditionMessage(e),'"'))
    }
)

# Tabulate loadings values
loadings.call <- selectVar(final.mbsplsda, comp = 1)$call$value
loadings.depth <- selectVar(final.mbsplsda, comp = 1)$depth$value

# Tabulate stability values
if (nrow(loadings.call) == 1)
{
  stability.call <- data.frame(
    "Var1" = rownames(loadings.call),
    "Freq" = 1
  )
} else {
  stability.call <- tryCatch(
    {
      as.data.frame(perf.final.mbsplsda$features$stable$nrep1$call$comp1)
    },
    error = function(e) {
      data.frame(
        "Var1" = rownames(loadings.call),
        "Freq" = rep(0, nrow(loadings.call))
      )
    }
  )
}
rownames(stability.call) <- stability.call$Var1

if (nrow(loadings.depth) == 1)
{
  stability.depth <- data.frame(
    "Var1" = rownames(loadings.depth),
    "Freq" = 1
  )
} else {
  stability.depth <- tryCatch(
    {
      as.data.frame(perf.final.mbsplsda$features$stable$nrep1$depth$comp1)
    },
    error = function(e) {
      data.frame(
        "Var1" = rownames(loadings.depth),
        "Freq" = rep(0, nrow(loadings.depth))
      )
    }
  )
}
rownames(stability.depth) <- stability.depth$Var1

# Identify and fix mismatches between loadings and stability values
## This is a bug in mixOmics which appears to be corrected in later versions but can occur prior to R 4.5
if (any(is.na(match(rownames(loadings.call), rownames(stability.call)))))
{
    stability.call <- data.frame(
        "Var1" = rownames(loadings.call),
        "Freq" = rep(0, nrow(loadings.call))
    )
    rownames(stability.call) <- stability.call$Var1
    print("WARNING: a bug in older mixOmics versions does not let us obtain 'call' stability values; program will otherwise continue without issues")
}

if (any(is.na(match(rownames(loadings.depth), rownames(stability.depth)))))
{
    stability.depth <- data.frame(
        "Var1" = rownames(loadings.depth),
        "Freq" = rep(0, nrow(loadings.depth))
    )
    rownames(stability.depth) <- stability.depth$Var1
    print("WARNING: a bug in older mixOmics versions does not let us obtain 'depth' stability values; program will otherwise continue without issues")
}

# Format call and depth values
stability.call <- stability.call[match(rownames(loadings.call), rownames(stability.call)),,drop=FALSE]
stability.call <- na.omit(stability.call)
loadings.call <- loadings.call[rownames(loadings.call) %in% rownames(stability.call),,drop=FALSE]
values.call <- cbind(stability.call, loadings.call)
values.call$type <- "call"

stability.depth <- stability.depth[match(rownames(loadings.depth), rownames(stability.depth)),,drop=FALSE]
stability.depth <- na.omit(stability.depth)
loadings.depth <- loadings.depth[rownames(loadings.depth) %in% rownames(stability.depth),,drop=FALSE]
values.depth <- cbind(stability.depth, loadings.depth)
values.depth$type <- "depth"

# Join loading and stability values
feature.details.table <- rbind(values.call, values.depth)

# Reformat loading values for ease of interpretation
feature.details.table$direction <- ifelse(feature.details.table$value.var > 0, "right", "left")
feature.details.table$value.var <- abs(feature.details.table$value.var)

# Split out the locations of the feature
location.details <- str_match(rownames(feature.details.table), "^(.+)_(\\d+)_(\\d+)$")
feature.details.table[c("chrom", "start", "end")] <- location.details[,2:4]

# Sort and rename table columns and rows
feature.details.table <- feature.details.table[,c("chrom", "start", "end", "type", "Freq", "value.var", "direction")]
colnames(feature.details.table) <- c("chrom", "start", "end", "type", "stability", "abs_loading", "direction")
feature.details.table <- feature.details.table[order(feature.details.table$stability, decreasing = TRUE),,drop=FALSE]

# Write selected variants to file
write.table(feature.details.table, file=args$os, sep="\t", row.names=FALSE, quote=FALSE)

# Produce an additional output to see what class sPLS-DA predicts each sample as belonging to
mbsplsda.pred <- predict(final.mbsplsda, splsda.X)
mbsplsda.pred <- as.data.frame(mbsplsda.pred$WeightedPredict[,,])
mbsplsda.pred$predicted_group <- ifelse(mbsplsda.pred$group1 >= mbsplsda.pred$group2, "group1", "group2")
mbsplsda.pred$true_group <- Y

write.table(mbsplsda.pred, file=args$op, sep="\t", row.names=TRUE, quote=FALSE, col.names=NA)
