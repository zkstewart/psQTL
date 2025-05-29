#!/usr/bin/env Rscript

if (!requireNamespace("argparser", quietly=TRUE)) {
  install.packages("argparser")
}
library(argparser)

fix_numeric_first_char <- function(x) {
  if (grepl("^[[:digit:]]+", x)) {
    return (paste0("X", x))
  } else {
    return (x)
  }
}

# Establish parser
p <- arg_parser("Run PLS-DA in windows across a genome, using sPLS-DA to select features that contribute to the model")

# Add arguments
p <- add_argument(p, "m", help="Metadata file", type = "character")
p <- add_argument(p, "v", help="Encoded VCF to process", type = "character")
p <- add_argument(p, "ov", help="Output file for selected variants", type = "character")
p <- add_argument(p, "ob", help="Output file for Balanced Error Rate (BER) by window", type = "character")
p <- add_argument(p, "or", help="Output file for Rdata objects", type = "character")
p <- add_argument(p, "--threads", help="Threads to use", type = "numeric", default = 1)
p <- add_argument(p, "--windowSize", help="Window size (default: 100000)", type="numeric", default=100000)
p <- add_argument(p, "--berCutoff", help="BER cutoff (default: 0.4)", type="numeric", default=0.4)
p <- add_argument(p, "--MAF", help="Minor Allele Frequency (MAF) filter (default: 0.05)", type="numeric", default=0.05)
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

# Identify mixOmics version
mixomicsVersion = unlist(packageVersion("mixOmics"))
if (mixomicsVersion[2] < 30) {
    print("Older version of mixOmics detected, which means that the plotLoadings function will emit an 'Rplots.pdf' file you can ignore")
}

# Parse metadata file
metadata.table <- read.table(file=args$m, header=FALSE, sep="\t", stringsAsFactors=FALSE)
metadata.table$V1 <- unlist(lapply(metadata.table$V1, fix_numeric_first_char))

# Read in encoded values
df <- read.table(gzfile(args$v), header=TRUE, sep="\t", na.strings=".")
maf.cutoff <- ceiling((ncol(df) - 2) * args$MAF)

df <- df[rowSums(df[,! colnames(df) %in% c("chrom", "pos")], na.rm=TRUE)>0,]
rownames(df) <- make.names(paste0(df$chrom, "_", df$pos), unique=TRUE)

# Drop any df values we are not analysing [Can occur if user metadata is a subset of VCF samples]
df <- df[,c("chrom", "pos", metadata.table$V1)] # this also sorts df and metadata equivalently

# Discover incompatibilities between metadata and encoded VCF
if ((ncol(df)-2) != nrow(metadata.table)) # -2 to account for c("chrom", "pos")
{
  stop(paste0("Encoded VCF column names (", paste(colnames(df[3:ncol(df)]), collapse=","), ") do not equal metadata sample labels (", paste(metadata.table$V1, collapse=","), "); incompatibility means sPLS-DA analysis cannot continue"))
}

# Extract Y variable values
Y <- metadata.table$V2

# Iterate over chromosomes and windows to run PLS-DA
window.explanation <- data.frame("chrom" = numeric(0), "pos" = numeric(0), "BER" = numeric(0))
selected.features <- data.frame("chrom" = numeric(0), "pos" = numeric(0), "BER" = numeric(0))
for (chromosome in unique(df$chrom))
{
  chromDF <- df[df$chrom == chromosome,]
  chromLength <- max(chromDF$pos)
  numWindows <- ceiling(chromLength / args$windowSize)
  for (windowIndex in 1:numWindows)
  {
    # Extract variants in window
    windowStart <- (windowIndex-1) * args$windowSize
    windowEnd <- windowIndex * args$windowSize
    windowDF <- chromDF[chromDF$pos >= windowStart & chromDF$pos < windowEnd,]
    windowDF <- windowDF[,! colnames(windowDF) %in% c("chrom", "pos"),drop=FALSE]
    
    # Drop any variants with NA values
    windowDF <- na.omit(windowDF)
    
    # Remove duplicate/linked variants
    windowDF <- windowDF[!duplicated(windowDF),]
    
    # Remove invariant sites
    windowDF <- windowDF[! rowSums(windowDF) %in% c(0, ncol(windowDF)),,drop=FALSE]
    
    # Store BER=0.5 if insufficient variants are found in this window
    if (nrow(windowDF) < 2)
    {
      window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, windowStart, 0.5) # BER=0.5
      next
    }
    
    # Transpose to obtain X value for PLS-DA
    X <- t(windowDF)
    
    # Run PLS-DA
    window.plsda <- plsda(X, Y,
                          ncomp = 1,
                          scale = FALSE,
                          max.iter = args$maxiters)
    
    # Assess model performance
    if (ncol(window.plsda$X) <= 6) {
      window.perf <- perf(window.plsda,
                          validation = "loo")
    } else {
      window.perf <- tryCatch(
        {
           perf(window.plsda,
                folds = 2, validation = "Mfold", 
                nrepeat = NREP)
        },
        error = function(e) {
          perf(window.plsda,
               validation = "loo")
        }
      )
    }
    window.ber <- window.perf$error.rate$BER[[1]]
    
    # Locate the features contributing most to PLS-DA outcome
    if (mixomicsVersion[2] >= 30) {
        window.loadings <- plotLoadings(window.plsda, comp=1, contrib = 'max', method = 'median', plot = FALSE)
    } else {
        window.loadings <- plotLoadings(window.plsda, comp=1, contrib = 'max', method = 'median')
    }
    
    window.contribution <- window.loadings$X[window.loadings$X$GroupContrib != "tie",]
    importance.cutoff <- abs(window.contribution[1,]$importance)/2 # arbitrary, but adequate
    window.features <- window.contribution[abs(window.contribution$importance) >= importance.cutoff,]
    
    # Store window explanatory power if features are found
    if (nrow(window.features) == 0)
    {
      window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, windowStart, 0.5) # BER=0.5
      next
    }
    window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, windowStart, window.ber)
    
    # Store best selected feature
    for (row.index in 1:nrow(window.features))
    {
      feature.row <- window.features[row.index,]
      feature.pos <- unlist(strsplit(rownames(feature.row), "_"))
      feature.pos <- as.numeric(feature.pos[length(feature.pos)])
      selected.features[nrow(selected.features) + 1,] <- c(chromosome, feature.pos, window.ber)
    }
  }
}
rownames(selected.features) <- paste0(selected.features$chrom, "_", format(as.numeric(selected.features$pos), scientific = FALSE, trim = TRUE))
selected.features$pos <- as.numeric(selected.features$pos)
selected.features$BER <- as.numeric(selected.features$BER)

# Filter down feature table to those selected in windows
selected.df <- df[match(rownames(selected.features[selected.features$BER <= args$berCutoff,]), rownames(df)),]
selected.X <- t(selected.df[, ! colnames(selected.df) %in% c("chrom", "pos")])

# Raise error if we've filtered out all features with BER cutoff
if (nrow(selected.X) == 0)
{
  original.number <- ncol(selected.df) - 2 # sans the chrom and pos columns
  stop(paste0(
    "BER cutoff of ", args$berCutoff, " reduces potential features from ",
    original.number, " down to 0. Your data either has insufficient information ",
    "to enable QTL prediction or your BER cutoff may be too strict."
  ))
}

# Tune sPLS-DA to choose number of genomic features
list.keepX <- c(1:9,  seq(10, 30, 5)) # it is very improbable that more than 30 QTLs exist or can be meaningfully identified
if (mixomicsVersion[2] <= 30) {
    tune.splsda.test <- tune.splsda(selected.X, Y, test.keepX = list.keepX,
                                    ncomp = 1, folds = 2, validation = 'Mfold',
                                    scale = FALSE,
                                    nrepeat = args$nrepeat, max.iter = args$maxiters,
                                    cpus = args$threads)
} else {
    tune.splsda.test <- tune.splsda(selected.X, Y, test.keepX = list.keepX,
                                    ncomp = 1, folds = 2, validation = 'Mfold',
                                    scale = FALSE,
                                    nrepeat = args$nrepeat, max.iter = args$maxiters,
                                    BPPARAM = BPPARAM)
}
ncomp <- 1 # 'tune.splsda.test$choice.ncomp$ncomp' won't be used as a single component is sufficient when discriminating two groups
select.keepX <- tune.splsda.test$choice.keepX[1:ncomp]

# Run sPLS-DA to select multi-QTLs and identify most important features
final.splsda <- splsda(selected.X, Y, keepX = select.keepX,
                       ncomp = ncomp,
                       scale = FALSE,
                       max.iter = args$maxiters)
if (mixomicsVersion[2] <= 30) {
    perf.final.splsda <- perf(final.splsda,
                              folds = 2, validation = "Mfold",
                              nrepeat = args$nrepeat,
                              cpus = args$threads)
} else {
    perf.final.splsda <- perf(final.splsda,
                              folds = 2, validation = "Mfold",
                              nrepeat = args$nrepeat,
                              BPPARAM = BPPARAM)
}

# Tabulate stability values
select.name <- selectVar(final.splsda, comp = 1)$name
stability.table <- as.data.frame(perf.final.splsda$features$stable$comp1[select.name])
stability.table <- na.omit(stability.table) # this might not be needed anymore, but is kept for safety

rownames(stability.table) <- stability.table$Var1
stability.table <- stability.table[,c("Freq"), drop=FALSE]
stability.table <- stability.table[order(stability.table$Freq, decreasing = TRUE),,drop=FALSE]

# Obtain loadings
splsda.loadings <- as.data.frame(final.splsda$loadings$X[,1,drop=FALSE])
splsda.loadings <- splsda.loadings[abs(rowSums(splsda.loadings)) > 0,,drop=FALSE]

# Reformat loading values for ease of interpretation
splsda.loadings$direction = ifelse(splsda.loadings$comp1 > 0, "right", "left")
splsda.loadings$comp1 = abs(splsda.loadings$comp1)

# Match loadings order to stability values
splsda.loadings <- splsda.loadings[match(rownames(stability.table), rownames(splsda.loadings)),,drop=FALSE]

# Join stability and loading values
feature.details.table <- cbind(stability.table, splsda.loadings)
feature.details.table[c("chrom", "pos")] <- do.call(rbind, strsplit(rownames(feature.details.table), "_(?=[^_]+$)", perl=TRUE))
feature.details.table <- feature.details.table[,c("chrom", "pos", "Freq", "comp1", "direction")]
colnames(feature.details.table) <- c("chrom", "pos", "stability", "abs_loading", "direction")

# Adjust BER to cap at 0.5
window.explanation.adjusted <- window.explanation
window.explanation.adjusted$BER <- as.numeric(window.explanation.adjusted$BER)
window.explanation.adjusted$BER[window.explanation.adjusted$BER > 0.5] <- 0.5 - (window.explanation.adjusted$BER[window.explanation.adjusted$BER > 0.5] %% 0.5)

# Write selected variants and BER to file
write.table(feature.details.table, file=args$ov, sep="\t", row.names=FALSE, quote=FALSE)
write.table(window.explanation.adjusted, file=args$ob, sep="\t", row.names=FALSE, quote=FALSE)

# Write R objects to file for later integration
save(selected.X, Y, feature.details.table, file=args$or)
