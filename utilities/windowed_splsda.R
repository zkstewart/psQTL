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

unify_accepted_y_values <- function(Y) {
  # Unify all variants of group1
  Y[Y == "bulk1"] <- "group1"
  Y[Y == "bulk 1"] <- "group1"
  Y[Y == "b1"] <- "group1"
  Y[Y == "g1"] <- "group1"
  Y[Y == "group 1"] <- "group1"
  Y[Y == "1"] <- "group1"
  
  # Unify all variants of group2
  Y[Y == "bulk2"] <- "group2"
  Y[Y == "bulk 2"] <- "group2"
  Y[Y == "b2"] <- "group2"
  Y[Y == "g2"] <- "group2"
  Y[Y == "group 2"] <- "group2"
  Y[Y == "2"] <- "group2"
  
  return (Y)
}

balanced_error_rate <- function(Ytrue, Ypred) {
    # Prevent incompatible comparison
    if (length(Ytrue) != length(Ypred)) {
        stop("Ytrue and Ypred must have the same length")
    }

    if (length(unique(Ypred)) == 1)
    {
        BER <- 0.5
    } else {
        BER <- mean(diag(1 - (table(Ytrue, Ypred) / table(Ytrue, Ytrue))))
    }
    return (BER)
}

lda_prediction <- function(Ytrue, lm.df) {
    lda.fit <- lda(Y ~ X, data=lm.df)
    lda.pred <- predict(lda.fit, lm.df)
    posterior <- as.data.frame(lda.pred$posterior)
    predicted <- ifelse(posterior$group1 >= posterior$group2, "group1", "group2")
    return (balanced_error_rate(Ytrue, predicted))
}

auto_rerun_tuning <- function(mixomicsVersion, selected.X, Y, list.keepX, nrep, maxiters, cpus, rerun=5) {
    # Automatically rerun tuning if it fails, up to 'rerun' times
    for (i in 1:rerun)
    {
        if (mixomicsVersion[2] <= 30)
        {
            tune.splsda.test <- tryCatch({
                tune.splsda(selected.X, Y, test.keepX = list.keepX,
                            ncomp = 1, folds = 2, validation = 'Mfold',
                            scale = FALSE,
                            nrepeat = nrep, max.iter = maxiters,
                            cpus = cpus)
            }, error = function(e) {
                NA
            })

            if (! all(is.na(tune.splsda.test))) {
                break
            }
        } else {
            tune.splsda.test <- tryCatch({
                tune.splsda(selected.X, Y, test.keepX = list.keepX,
                            ncomp = 1, folds = 2, validation = 'Mfold',
                            scale = FALSE,
                            nrepeat = nrep, max.iter = maxiters,
                            BPPARAM = cpus)
            }, error = function(e) {
                NA
            })

            if (! all(is.na(tune.splsda.test))) {
                break
            }
        }
    }

    # Get tuned result if it worked, otherwise default to 1 feature
    if (! all(is.na(tune.splsda.test))) {
        select.keepX <- tune.splsda.test$choice.keepX[1:1] # 'tune.splsda.test$choice.ncomp$ncomp' won't be used as a single component is sufficient when discriminating two groups
    } else {
        select.keepX <- 1 # just pick the best feature if tuning isn't working
    }

    return (select.keepX)
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

if (!requireNamespace("dplyr", quietly=TRUE))
    install.packages("dplyr")
library(dplyr)

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

df <- df[rowSums(df[,! colnames(df) %in% c("chrom", "pos")], na.rm=TRUE)>=maf.cutoff,]
rownames(df) <- make.names(paste0(df$chrom, "_", df$pos), unique=TRUE)

# Drop any metadata samples not present in VCF
metadata.table <- metadata.table[metadata.table$V1 %in% colnames(df),,drop=FALSE]

# Extract and unify Y variable values to 'group1' and 'group2'
Y <- metadata.table$V2
Y <- unify_accepted_y_values(Y)

# Discover issues with Y variable
if (length(unique(Y)) != 2)
{
    stop(paste0("Y variable must have exactly two unique values, but found ", length(unique(Y)), " unique values: ", paste(unique(Y), collapse=", ")))
}
if (length(union(c("group1", "group2"), Y)) != 2)
{
    stop(paste0("Y variable is expected to only have 'group1' and 'group2' values, but found : ", paste(unique(Y), collapse=", ")))
}

# Drop any df values we are not analysing [Can occur if user metadata is a subset of VCF samples]
df <- df[,c("chrom", "pos", metadata.table$V1)] # this also sorts df and metadata equivalently

# Check that VCF has at least 1 variant
if (nrow(df) == 0)
{
    stop(paste0("Encoded VCF file '", args$v, "' has no variants. sPLS-DA analysis cannot continue; please check your input VCF file for issues."))
}

# Discover incompatibilities between metadata and encoded VCF
if ((ncol(df)-2) != nrow(metadata.table)) # -2 to account for c("chrom", "pos")
{
    stop(paste0("Encoded VCF column names (", paste(colnames(df[3:ncol(df)]), collapse=","), ") do not equal metadata sample labels (", paste(metadata.table$V1, collapse=","), "); incompatibility means sPLS-DA analysis cannot continue"))
}

# Iterate over chromosomes and windows to run PLS-DA
window.explanation <- list()
window.index <- 1
selected.features <- list()
selected.index <- 1
for (chromosome in unique(df$chrom))
{
  chromDF <- df[df$chrom == chromosome,]
  chromLength <- max(chromDF$pos)+1 # +1 to ensure we include the last position in the chromosome
  numWindows <- ceiling(chromLength / args$windowSize)
  for (windowIndex in 1:numWindows)
  {
    # Extract variants in window
    windowStart <- (windowIndex-1) * args$windowSize # 0-based index
    windowEnd <- windowIndex * args$windowSize
    windowDF <- chromDF[chromDF$pos >= windowStart & chromDF$pos < windowEnd,]
    windowDF <- windowDF[,! colnames(windowDF) %in% c("chrom", "pos"),drop=FALSE]
    
    # Drop any variants with NA values
    windowDF <- na.omit(windowDF)
    
    # Remove duplicate/linked variants
    windowDF <- windowDF[!duplicated(windowDF),]
    
    # Remove invariant sites
    windowDF <- windowDF[! rowSums(windowDF) %in% c(0, ncol(windowDF)),,drop=FALSE]
    
    # Skip windows with no variant presence
    if (nrow(windowDF) == 0)
    {
        window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=0.5)
        window.index <- window.index + 1
        next
    }
    
    # Perform alternate process if only 1 variant is found in this window
    if (nrow(windowDF) == 1)
    {
        # Extract details for this single-feature window
        feature.row <- windowDF[1,]
        feature.pos <- unlist(strsplit(rownames(feature.row), "_"))
        feature.pos <- as.numeric(feature.pos[length(feature.pos)])
        lm.df <- data.frame("X" = as.numeric(feature.row), "Y" = Y)
        
        # Run LDA, or skip if complete similarity or segregation exists
        group1Values <- unique(lm.df[lm.df$Y == "group1",]$X)
        group2Values <- unique(lm.df[lm.df$Y == "group2",]$X)
        if (length(group1Values) == 1 & length(group2Values) == 1)
        {
            if (group1Values == group2Values) {
                window.ber <- 0.5 # cannot distinguish in an identical region
            } else {
                window.ber <- 0 # complete segregation leads to perfect prediction
            }
        } else {
            # Use LDA to predict BER for this window
            window.ber <- tryCatch({
                    lda_prediction(Y, lm.df)
                }, error = function(e) {
                    if (e$message == "group means are numerically identical") {
                        return (0.5)
                    } else {
                        stop(paste0("Encountered unhandled error: ", e$message))
                    }
                }
            )
        }
        
        # Store window and feature
        window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=window.ber)
        window.index <- window.index + 1

        selected.features[[selected.index]] <- cbind(data.frame(chrom=chromosome, pos=feature.pos, BER=window.ber), feature.row)
        selected.index <- selected.index + 1
        next
    }

    # Detect scenario where only 1 variant site exists in a population
    if (sum(colSums(windowDF) > 0) < 2)
    {
        window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=0.5)
        window.index <- window.index + 1
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
    window.perf <- tryCatch(
        {
            if (ncol(window.plsda$X) <= 6)
            {
                perf(window.plsda, validation = "loo")
            } else {
                tryCatch(
                    {
                        if (mixomicsVersion[2] <= 30) {
                            perf(window.plsda,
                                 folds = 2, validation = "Mfold",
                                 nrepeat = args$nrepeat,
                                 cpus = args$threads)
                        } else {
                            perf(window.plsda,
                                 folds = 2, validation = "Mfold",
                                 nrepeat = args$nrepeat,
                                 BPPARAM = BPPARAM)
                        }
                    },
                    error = function(e) {
                        perf(window.plsda, validation = "loo")
                    }
                )
            }
        },
        error = function(e) {
            NA
        }
    )

    # If performance could not be calculated, skip this window
    if (all(is.na(window.perf))) {
        window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=0.5)
        window.index <- window.index + 1
        next
    }
    window.ber <- window.perf$error.rate$BER[[1]] # since we could calculate performance, we can extract BER
    
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
        window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=0.5)
        window.index <- window.index + 1
        next
    }
    window.explanation[[window.index]] <- data.frame(chrom=chromosome, pos=windowStart, BER=window.ber)
    window.index <- window.index + 1

    # Set up df for feature encodings and ensure its compatibility with the feature positions
    encodings.df <- windowDF[rownames(window.features),]
    if (nrow(encodings.df) != nrow(window.features))
    {
        stop("Error with ordering of encoded and selected features")
    }

    # Store best selected feature
    for (row.index in 1:nrow(window.features))
    {
        feature.row <- window.features[row.index,]
        feature.pos <- unlist(strsplit(rownames(feature.row), "_"))
        feature.pos <- as.numeric(feature.pos[length(feature.pos)])

        encodings.row <- encodings.df[row.index,]
        selected.features[[selected.index]] <- cbind(data.frame(chrom=chromosome, pos=feature.pos, BER=window.ber), encodings.row)
        selected.index <- selected.index + 1
    }
  }
}
window.explanation <- as.data.frame(dplyr::bind_rows(window.explanation))
selected.features <- as.data.frame(dplyr::bind_rows(selected.features))
selected.features$BER <- as.numeric(selected.features$BER)

# Filter down feature table to those selected in windows
original.number <- nrow(selected.features) # keep track of this for later diagnostic error information
selected.features <- selected.features[selected.features$BER <= args$berCutoff,,drop=FALSE]

# Raise error if we've filtered out all features with BER cutoff
if (nrow(selected.features) == 0)
{
    stop(paste0(
        "BER cutoff of ", args$berCutoff, " reduces potential features from ",
        original.number, " down to 0. Your data either has insufficient information ",
        "to enable QTL prediction with sPLS-DA or your BER cutoff may be too strict."
    ))
}

# Format features for sPLS-DA
selected.df <- selected.features[, ! colnames(selected.features) %in% c("chrom", "pos", "BER")] # drop non-feature columns
selected.df <- selected.df %>%
  mutate(across(everything(), as.numeric)) # convert feature values back into numeric, since it gets changed along the way

selected.X <- t(selected.df)
colnames(selected.X) <- paste0(selected.features$chrom, "_", selected.features$pos)

# Tune sPLS-DA to choose number of genomic features
if (ncol(selected.X) < 2)
{
    print("Only one feature selected, so no tuning of sPLS-DA is performed")
    select.keepX <- ncol(selected.X) # this will be 1
    tune.splsda.test <- NULL
} else {
    # Identify the number of features to test
    if (ncol(selected.X) < 10) {
        list.keepX <- c(1:ncol(selected.X)) # if fewer than 10 features, test all
    } else {
        max.features <- ifelse(ncol(selected.X) < 30, ncol(selected.X), 30) # it is very improbable that more than 30 QTLs exist or can be meaningfully identified
        list.keepX <- c(1:9,  seq(10, max.features, 5))
    }

    # Run tuning of sPLS-DA, with automatic rerun if it fails with 'if (diff.value < tol | iter > max.iter)'
    if (mixomicsVersion[2] <= 30) {
        select.keepX <- auto_rerun_tuning(mixomicsVersion, selected.X, Y, list.keepX, args$nrepeat, args$maxiters, args$threads)
    } else {
        select.keepX <- auto_rerun_tuning(mixomicsVersion, selected.X, Y, list.keepX, args$nrepeat, args$maxiters, BPPARAM)
    }
}

# Run sPLS-DA to select multi-QTLs and identify most important features
final.splsda <- splsda(selected.X, Y, keepX = select.keepX,
                       ncomp = 1,
                       scale = FALSE,
                       max.iter = args$maxiters)
tryCatch(
    {
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
    },
    error = function(e) {
        print(paste0("Warning: stability cannot be estimated due to error \"", conditionMessage(e),'"'))
    }
)

# Obtain loadings
splsda.loadings <- as.data.frame(final.splsda$loadings$X[,1,drop=FALSE])
splsda.loadings <- splsda.loadings[abs(rowSums(splsda.loadings)) > 0,,drop=FALSE]

# Tabulate stability values
if (nrow(splsda.loadings) == 1)
{
    stability.table <- data.frame(
        "Var1" = rownames(splsda.loadings),
        "Freq" = 1
    )
} else {
    stability.table <- tryCatch(
    {
        as.data.frame(perf.final.splsda$features$stable$comp1[
        selectVar(final.splsda, comp = 1)$name])
    },
    error = function(e) {
        data.frame(
        "Var1" = rownames(splsda.loadings),
        "Freq" = rep(0, nrow(splsda.loadings))
        )
    })
}
stability.table <- na.omit(stability.table) # this might not be needed anymore, but is kept for safety
rownames(stability.table) <- stability.table$Var1

stability.table <- stability.table[,c("Freq"), drop=FALSE]
stability.table <- stability.table[order(stability.table$Freq, decreasing = TRUE),,drop=FALSE]
stability.table <- na.omit(stability.table)

# Reformat loading values for ease of interpretation
splsda.loadings$direction = ifelse(splsda.loadings$comp1 > 0, "right", "left")
splsda.loadings$comp1 = abs(splsda.loadings$comp1)

# Match loadings order to stability values
splsda.loadings <- splsda.loadings[rownames(stability.table),,drop=FALSE]
splsda.loadings <- na.omit(splsda.loadings)

# Make sure stability and loading values match up
stability.table <- stability.table[rownames(splsda.loadings),,drop=FALSE]

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
