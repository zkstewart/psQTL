#!/usr/bin/env Rscript

if (!requireNamespace("argparser", quietly=TRUE)) {
  install.packages("argparser")
}
library(argparser)
#library(tidyr)
#library(dplyr)

# Establish parser
p <- arg_parser("Run PLS-DA in windows across a genome, selecting features that contribute to the model")

# Add arguments
p <- add_argument(p, "m", help="Metadata file", type = "character")
p <- add_argument(p, "v", help="Encoded VCF to process", type = "character")
p <- add_argument(p, "ov", help="Output file for selected variants", type = "character")
p <- add_argument(p, "ob", help="Output file for Balanced Error Rate (BER) by window", type = "character")
p <- add_argument(p, "--windowSize", help="Window size (default: 100000)", type="numeric", default=100000)
p <- add_argument(p, "--berCutoff", help="BER cutoff (default: 0.4)", type="numeric", default=0.4)
p <- add_argument(p, "--MAF", help="Minor Allele Frequency (MAF) filter (default: 0.05)", type="numeric", default=0.05)

# Parse the command line arguments
args <- parse_args(p)

# Load in remaining libraries
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("mixOmics", quietly=TRUE)) {
    BiocManager::install("mixOmics")
}
library(mixOmics)

# Parse metadata file
metadata.table <- read.table(file=args$m, header=FALSE, sep="\t", stringsAsFactors=FALSE)

# Read in encoded values
df <- read.table(gzfile(args$v), header=TRUE, sep="\t")
maf.cutoff <- ceiling((ncol(df) - 2) * args$MAF)

df <- df[rowSums(df[,! colnames(df) %in% c("chrom", "pos")], na.rm=TRUE)>0,]
rownames(df) <- make.names(paste0(df$chrom, "_", df$pos), unique=TRUE)

# Order metadata to match VCF rownames
metadata.table <- metadata.table[match(colnames(df)[! colnames(df) %in% c("chrom", "pos")], metadata.table$V1),]

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
    
    # Drop any variants with NA values
    windowDF <- na.omit(windowDF)
    
    # Remove duplicate/linked variants
    windowDF <- windowDF[!duplicated(windowDF[,! colnames(windowDF) %in% c("chrom", "pos")]),]
    
    # Store BER=0.5 if insufficient variants are found in this window
    if (nrow(windowDF) < 2)
    {
      window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, (windowEnd + windowStart) / 2, 0.5) # BER=0.5
      next
    }
    
    # Transpose to obtain X value for PLS-DA
    X <- t(windowDF[, ! colnames(windowDF) %in% c("chrom", "pos")])
    
    # Run PLS-DA
    window.plsda <- plsda(X, Y, ncomp = 1)
    window.perf <- perf(window.plsda, folds = 2, validation = "Mfold", 
                        dist = "max.dist", progressBar = FALSE, nrepeat = 10)
    window.ber <- window.perf$error.rate$BER[[1]]
    
    # Locate the features contributing most to PLS-DA outcome
    window.loadings <- plotLoadings(window.plsda, comp=1, contrib = 'max', method = 'median')
    window.contribution <- window.loadings$X[window.loadings$X$GroupContrib != "tie",]
    importance.cutoff <- abs(window.contribution[1,]$importance)/2 # arbitrary, but adequate
    window.features <- window.contribution[abs(window.contribution$importance) >= importance.cutoff,]
    
    # Store window explanatory power if features are found
    if (nrow(window.features) == 0)
    {
      window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, (windowEnd + windowStart) / 2, 0.5) # BER=0.5
      next
    }
    window.explanation[nrow(window.explanation) + 1,] <- c(chromosome, (windowEnd + windowStart) / 2, window.ber)

    # Store selected features
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

# Tune sPLS-DA to choose number of genomic features
list.keepX <- c(1:9,  seq(10, 30, 5)) # it is very improbable that more than 30 QTLs exist or can be meaningfully identified
tune.splsda.test <- tune.splsda(selected.X, Y, ncomp = 2, validation = 'Mfold',
                                folds = 2, dist = 'max.dist',
                                test.keepX = list.keepX, nrepeat = 10)
ncomp <- 1 # 'tune.splsda.test$choice.ncomp$ncomp' won't be used as a single component is sufficient when discriminating two groups
select.keepX <- tune.splsda.test$choice.keepX[1:ncomp]

# Run sPLS-DA to select multi-QTLs and identify most important features
final.splsda <- splsda(selected.X, Y, ncomp = ncomp, keepX = select.keepX)
perf.final.splsda <- perf(final.splsda, folds = 2, validation = "Mfold",
                          dist = "max.dist", progressBar = FALSE, nrepeat = 10)

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
