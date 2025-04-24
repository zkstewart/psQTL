# Table of Contents
- [Getting started](#getting-started)
- [Introduction](#introduction)
  - [Traditional bulked segregant analysis](#traditional-bulked-segregant-analysis)
  - [Per-sample segregant analysis](#per-sample-segregant-analysis)
  - [What psQTL does](#what-psqtl-does)
- [Installation](#installation)
- [How to use](#how-to-use)
  - [Formatting a metadata file](#formatting-a-metadata-file)
  - [Calling variants](#calling-variants)
  - [Predicting deletions](#predicting-deletions)
  - [Plotting results](#plotting-results)
  - [Reporting gene proximity](#reporting-gene-proximity)
- [A typical analysis pipeline](#a-typical-analysis-pipeline)
- [How to cite](#how-to-cite)

# Getting started
```
# Download this repository
git clone https://github.com/zkstewart/psQTL.git

# Make sure you have all prerequisite Python packages and external programs available
## Refer to the websites indicated in the 'Installation' section below.

# 'Prep'are for your analysis using psQTL_prep.py
python psQTL/psQTL_prep.py initialise -d /location/of/your/working/directory \
    --meta metadata.tsv \
    --bam /location/of/bam_files \
    --bamSuffix .sorted.bam
python psQTL/psQTL_prep.py call -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    --qual 30 --threads 12
python psQTL/psQTL_prep.py depth -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    --windowSize 1000 --threads 12

# 'Proc'ess your variant and deletion predictions and calculate segregation statistics
python /location/of/psQTL_proc.py call -d /location/of/your/working/directory --ignoreIdentical
python /location/of/psQTL_proc.py depth -d /location/of/your/working/directory

# 'Post'-processing steps including data plotting and report tabulation
python /mnt/c/git/psQTL/psQTL_post.py plot -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    -i call \ # or depth
    -o /location/to/write/output.pdf \
    -p line scatter histogram coverage genes \
    --annotation location/of/genome.gff3
python /mnt/c/git/psQTL/psQTL_post.py report -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    -i call \ # or depth
    -o /location/to/write/output.tsv \
    -a location/of/genome.gff3
```

# Introduction
## Traditional bulked segregant analysis
A variety of methods and pipelines exist for performing a bulked segregant analysis (BSA). In this type of experiment, organisms will be partitioned (*segregated*) into two populations (*bulks*) according to some differing phenotype, then tissue samples will be obtained from all individuals of each population. Those samples are pooled together and DNA is extracted together. Data analysis aims to identify differences in the allele frequencies of these two DNA pools to determine if there is a systematic difference that would point to the existence of one or more quantitative trait loci (QTL).

Historically, this has been necessary due to the cost of sequencing, especially when populations consist of dozens or hundreds of individuals. However, there are biases associated with this approach including but not limited to:
1. Unequal amounts of DNA obtained from each individual sample within the pooled sample may skew the allele frequency especially if the overrepresented sample is of a different genotype than most of its peers.
2. It can be hard or impossible to tell whether a population is a mix of homozygous reference (0/0) and homozygous alternate (1/1) alleles, or is heterozygous for that allele (0/1).

To address this second issue, it is common that parents would be separately sequenced to provide insight into the likely genotype of offspring. However, it is not uncommon for parent samples to be unavailable or for the true parents of some organisms to be unknown.

## Per-sample segregant analysis
Reduced sequencing costs open up the possibility for sequencing each individually obtained tissue sample. The benefits of doing so include but are not limited to:
1. The lack of sample pooling means that unequal DNA bias can be eliminated. This means that allele frequency does not need to be *estimated*, it can be *known* for the populations.
2. Samples can be individually genotyped, and we can know the exact proportion of the population that is of the varying homozygous or heterozygous genotypes.

Because of the specific knowledge of each individual's genotype, we can run an analysis without parent samples being available. And, we can use statistics that benefit from the knowledge of each sample's genotype rather than the collective allele frequency, in order to obtain greater power when predicting QTLs.

## Segregation at deletion sites
It's well established that SNPs and small indels are known to be responsible for influencing phenotype in QTLs. However, large deletions (such as those that deactivate or eliminate genes) can also be a major contributor to phenotypic difference among organisms. If such a deletion is responsible for the phenotypic segregation in two populations, you might identify it when analysing variant site segregation only through the variant's linkage to the deletion. This isn't guaranteed however, and hence it can be useful to specifically analyse deletions and how their occurrence segregates between populations.

## What psQTL does
psQTL runs a per-sample segregant analysis (PSA) to improve the statistical power of predicting QTLs relative to traditional BSA approaches. It offers three modules to streamline a PSA experiment, including:
1. It will prepare your data which includes predicting variants and deleted regions with `psQTL_prep.py`.
2. It will process your data to statistically predict variants and/or deletions that segregate between two populations with varying phenotype with `psQTL_proc.py`.
3. It will plot and report your results to allow for interrogation and understanding of the likely QTL(s) and what genes may be associated with the QTL using `psQTL_post.py`.

# Installation
psQTL is provided as a collection of Python 3.X scripts which make use of the following Python packages:
- biopython (https://biopython.org/)
- numpy (https://numpy.org/)
- pandas (https://pandas.pydata.org/)

It calls upon the HTSlib (http://www.htslib.org/) programs including:
- samtools
- bcftools
- bgzip

As well as vt (https://genome.sph.umich.edu/wiki/Vt).

psQTL should work on any system which can run Python, HTSlib, and vt; this likely precludes its use on Windows except if you use Windows Subsystem for Linux (WSL) which is compatible with these software.

You should make sure the Python packages and external program dependencies are installed and locatable within your system PATH variable.

# How to use
The psQTL pipeline is intended to proceed from **prep**aration, to **proc**essing, to **post**-processing as accomplished by the correspondingly named Python script files.

All scripts will return help information on the command line in order to assist you in running them like:

`python /location/of/psQTL_prep.py -h`

All parameters provided with a single dash e.g., `-d /location/to/run/the/analysis` need to be provided each time you call the program. Other parameters specified with double dashes e.g., `--meta metadata.tsv` are remembered or *cached* by psQTL within the analysis directory and hence only need to be provided once since the cached value is retrieved on subsequent use of psQTL ***in the analysis folder***. In other words, each analysis directory has its own cached values and parameters won't carry over between separate analyses.

## Directory setup
For any psQTL analysis, you must first *initialise* a working directory like:

`python /location/of/psQTL_prep.py initialise -d /location/to/run/the/analysis`

Using the *initialise* function, or when using the *call* or *depth* functions afterwards, you must indicate a metadata file and one of either:
1. One or more BAM files or directories containing BAM files, in order to facilitate variant calling and/or deletion prediction.
2. A VCF file you've already produced using your own methodology.

BAM files can be produced using any conventionally applied process of mapping reads to a genome, so long as they are sorted appropriately. And any valid VCF file should be acceptable.

## Formatting a metadata file
The metadata file is a simple tab-delimited text file (.TSV) with **two columns** and **no header**. The left column must contain your sample identifiers, with the right column indicating which *bulk* the sample belongs to with the possible options being *bulk1* or *bulk2*. It doesn't matter which phenotype is labelled as *bulk1* or *bulk2* as this will not influence your results - you just need to partition your samples into two populations.

For the sample identifiers, if you've provided BAM files you should ensure that the start of each file begins with that sample identifier and is immediately followed by a consistent file suffix. For a sample termed `sample1`, you should ensure that you have a BAM file that looks like `sample1{bamSuffix}` where `{bamSuffix}` is provided on the command line using the `--bamSuffix` option.

If you've provided a VCF file, then these sample identifiers should be the same as found in your VCF header line that begins with `#CHROM`.

You *can* provide a metadata file that *doesn't exactly match* your VCF or BAM files. If you do so, the VCF filtering will only consider the samples contained within your metadata file, and results produced by `psQTL_proc.py` will only use the indicated samples when calculating the segregation statistics. However, I think you should avoid doing so unless you're intentionally using this mechanism to filter your results. The program will warn you when there's a discrepancy.

## Calling variants
Variant calling is optionally assisted through the `psQTL_prep.py call` function. It automates a process involving *bcftools mpileup* followed by *bcftools call* with variant normalisation steps involving *bcftools norm* and *vt decompose_blocksub*.

The options you'll need to specify to allow this to occur are:
- BAM files or directories containing BAM files through `--bam`.
- A BAM suffix through `--bamSuffix`.
  - Directories provided to `--bam` will be checked for files ending with this suffix, and will be included in the analysis.
- A variant quality score to filter on through `--qual`.
  - The recommended value is `30` but you may increase this to make it more strict, or decrease it to reduce the amount of filtration.

Variant calling will produce a raw and filtered VCF file with standard formatting for downstream analysis.

## Predicting deletions
Deletion prediction is optionally assisted through the `psQTL_prep.py depth` function. It automates a process involving *samtools depth* from which position depth values are summed within non-overlapping window regions into a histogram. The median depth value of each chromosome within each sample is determined, and a simple heuristic employed to identify regions without deletions, as well as regions with hemizygous or homozygous deletion.

The options you'll need to specify to allow this to occur are:
- The `--bam` and `--bamSuffix` arguments; see notes in [Calling variants](#calling-variants) above.
- The bin or window size to sum read alignments within using `--windowSize`.
  - We ideally want to pick a size that isn't *too small* and hence subject to bias from random data variations, but also not *too big* such that real deletions get hidden by their adjacent non-deleted regions.
  - The recommended size is `1000` which should balance these goals.

The deletion prediction will produce a VCF-*like* file. It does not contain the normal information fields that a true VCF has, but it conforms to VCF styling from the `#CHROM` line onwards albeit with dummy values used for the `REF` and `ALT` fields. Genotypes are encoded as `0/0` for homozygous non-deletion, `0/1` for hemizygous deletion, and `1/1` for homozygous deletion.

## Calculating Euclidean Distance for QTL segregation
The variant calls or deletion predictions are **proc**essed with the `psQTL_proc.py` script. It calculates the Euclidean distance (ED) between genotypes across populations; deletion predictions are encoded as genotypes as noted above. The result provides an indication of how much the two bulks differ at each variant or deletion region, with large values indicating more segregation and values of 0 indicating an identical allele frequency between the two bulks.

When handling deletion predictions, no additional options are necessary.

When handling variant calls, you may specify the following:
- The `--ignoreIdentical` argument will specify whether you would like to ignore variant calls if all samples have an identical genotype e.g., they all have the `1/1` genotype. These variants can end up in your VCF file if all of your samples differ to the reference genome. Hence, since these *variants* are not actually *varying* across your samples, it may be best to ignore them.

The output file is in TSV format, with six columns denoted by a header row. These columns include:
- CHROM
  - The contig identifier from your genome FASTA file.
- POSI
  - The position (in bp) of your variant OR the starting position of the window that depth values were summarised into.
- variant
  - A simple description of what type of variant is being described: 'snp' or 'indel'
- bulk1_alleles
  - The number of diploid alleles which were genotyped at this position for the samples belonging to bulk 1. Since this program assumes diploidy, this is equivalent to `number of genotyped individuals * 2`
- bulk2_alleles
  - As above, but for bulk 2.
- euclideanDistance
  - The ED calculated at this position for the bulk genotype segregation.

## Plotting results
**Proc**essed results can be plotted using the `psQTL_post.py` script. It can plot variant calls using the `-i call` option, or deletion predictions using the `-i depth` option.

Importantly, you can raise ED values to a power using `--power`. By default the power value is `4` which helps to minimise noise and amplify signal, which is important when plotting smoothed ED values.

Plotting can be done on all chromosomes, on a selection of one or more specified chromosomes, or in bounded regions within specified chromosomes using the `--regions` option. For example:
- If you omit the `--regions` option, **all** chromosomes will be plotted.
- If you just provide one or more values e.g. `--regions chr1 chr2`, then you will plot only the specified chromosomes.
- You can provide ranges e.g., `--regions chr1:100000-250000 chr2:2000000-3000000` to only plot within the specified regions of those chromosomes.

Variants can be filtered at this stage to allow for exploration of important trends. Specifically, the proportion of variants which are allowed to be ungenotyped in the population bulks is controlled by `--missing`.
- This value can range from `0` (0% are allowed to be missing; **strict**) to `1` (100% are allowed to be missing; **relaxed**) with the default recommendation being `0.5`.
- If you specify a value of `0.5`, then we will filter out a variant if it's missing (i.e., hasn't been genotyped) in more than 50% of the samples in **either bulk**. Hence, if a variant has no samples genotyped (100% missing) in bulk 1, regardless of what we see in bulk 2, we will omit the variant from results plotting.

Lastly, it allows one to flexibly produce various plot types including:

- Line plot
  - The line plot shows the smoothed ED values across the genome, which can help to visualise peaks or troughs in segregation statistics.
  - The `--wma` option controls how much smoothing is applied using a Weighted Moving Average of the segregation statistic.
- Scatter plot
  - Plots the same underlying data as the line plot, but without smoothing. The line plot can be overlayed on top of the scatter plot for combined visualisation.
- Histogram plot
  - The histogram plot helps to reveal where highly segregating variants or deletion regions are most common by counting their occurrence within window/bin regions. It provides an additional check against the line plot to see if a peak is occurring because of only a small number of highly segregating variants, or if there is a large number of variants.
  - The `--bin` option controls how large the bin region is for counting the occurrence of highly segregating variants or deletion regions.
  - The `--threshold` option controls what ED value is considered to be *highly segregating* for counting in bin regions.
- Coverage plot
  - The coverage plot helps to reveal how each bulk's read depth varies across the genome. Specifically, it might show where deletions occur in one bulk but not the other.
  - The `--sampleCoverage` option lets you plot one or more samples as individual lines in the coverage plot. This might help you to see what the parent(s) of your segregating phenotypes look like when it comes to copy number variants or deletions.
- Gene plot
  - The gene plot shows the location of gene models within plotted regions. It may help to identify which genes occur in regions where segregation statistics peak or trough.
  - The `--annotation` option lets you input a GFF3 file containing annotations which will be parsed and presented.

## Reporting gene proximity
TBD

# A typical analysis pipeline
## Preparation
Prior to psQTL, you should map your reads to your genome however you normally do so, getting sorted BAM results. An example is indicated for a single sample below, but any popular program should work.

```
bwa mem -t <CPUS> -R <readgroup> <genomeFASTA> <forwardreads> <reverseread> > sample1.sam
samtools sort -@ <CPUS> -O bam -o sample1.sorted.bam -m <MEM>G sample1.sam
bamtools index sample1.sorted.bam
```

Next, you may optionally want to mark PCR duplicates. The necessity of this is debateable, but doing so is unlikely to hurt.

```
picard MarkDuplicates I=sample1.sorted.bam O=sample1.sorted.md.bam M=sample1.md_metrics.txt
bamtools index sample1.sorted.md.bam
```

Afer doing this for all of your samples, you next want to format a metadata file for the analysis. As detailed [earlier](#formatting-a-metadata-file) this file should be **tab-delimited** with **two columns** and **no header**. The left column indicates your sample IDs, and the right indicates which of two bulks it belongs to. Note that the sample ID should be the same as your BAM files, just with the suffix removed. So for example, if your file is `sample1.sorted.md.bam` then your sample ID is `sample1` and your suffix is `.sorted.md.bam`.

```
sample1 bulk1
sample2 bulk2
sample3 bulk2
```

Now you're ready to use `psQTL_prep.py` to **initialise** your working directory. Using this, we'll instruct psQTL to set up a new folder/directory for an analysis of your mapped BAM files using the metadata file you've got ready to go. Note that we're going to name our working directory `psqtl_analysis`, our metadata file is named `metadata.tsv`, and our BAM files are all located in `/location/of/bam_files` with the suffix `.sorted.md.bam`. It's important that *all* of your BAM files have that exact same suffix!

`python /location/of/psQTL_prep.py initialise -d psqtl_analysis --meta metadata.tsv --bam /location/of/bam_files --bamSuffix .sorted.md.bam`

The program will tell you that some values have been stored in various *caches* - these *caches* are what psQTL uses to remember where your file are located and what options you've specified.

With the initialisation done, we're ready to **call** variants now. The program will remember where your metadata and BAM files are located, but you will also need to let psQTL know where the genome FASTA file is located; this should be the exact same file as you used when mapping reads. Lastly, we need to provide some information on how we want to filter our variants in terms of their **qual**ity scores. Our genome file will be called `/location/of/genome.fasta`. We'll use the recommended value of `30` for quality scores and we'll also use `12` **threads** to speed things along.

`python /location/of/psQTL_prep.py call -d psqtl_analysis -f /location/of/genome.fasta --qual 30 --threads 12`

We should expect this process to take a while. If we're interested in predicting genomic deletions to see if they appear to be QTLs as well, we can use the **depth** functionality of psQTL to do that. We'll indicate the same genome FASTA file as before, and let psQTL know what window size we want to predict deletions in.

`python /location/of/psQTL_prep.py depth -d psqtl_analysis -f /location/of/genome.fasta --windowSize 1000 --threads 12`

Once the analyses are complete, we can **view** the analysis folder for some details on what we've predicted. Expect this to tell you what your files are, how many variants you called (before and after filtering), and how many windows in your genome appear to contain a deletion in one or more of your samples.

`python /location/of/psQTL_prep.py view -d psqtl_analysis`

## Processing
Either using psQTL or independently yourself, you should have one or both of 1) a variant VCF file, and 2) a deletion VCF-like file. We want to process the file(s) and calculate how much the two bulks segregate at variants or genomic deletions. The statistic used to measure that is the Euclidean distance, and we can calculate that for the variant **call** file with:

`python /location/of/psQTL_proc.py call -d psqtl_analysis`

That will produce the `psQTL_variants.ed.tsv.gz` file within your analysis directory which contains the raw Euclidean distance data for each variant in your VCF file.

Alternatively, to look at the segregation at genomic window regions using the deletion VCF-like file, we can do it with:

`python /location/of/psQTL_proc.py depth -d psqtl_analysis`

This will produce the `psQTL_depth.ed.tsv.gz` file with the raw Euclidean distance data for segregation at each genomic window region.

## Post-processing
### Plotting
After running the above processing step, you are now ready to plot the results. psQTL relies upon the visual inspection of results, with a usual process involving you plotting ED statistics across the entire genome, then manually looking for regions that have elevated segregation. You can then plot just those regions to give you a "zoomed in" look at the statistics, and iterate upon this process until you are happy with what you've found. This can involve changes in how you count variants with the histogram plot such as changes to the bin size or threshold for counting variants.

To begin, you might plot all the chromosomes for variant calls like:

```
python /location/of/psQTL_post.py plot -d psqtl_analysis -i call \
    -f /location/of/genome.fasta -o plot.pdf \
    -p line histogram
```

This will give you a plot from which you can look for peaks in the ED statistic on your line plot, or for peaks where larger numbers of variants occur in your histogram. Let's say you find that `chr2` is of interest, but that you have a lot of variants with large ED values. Let's plot just that chromosome, increase the threshold for counting a variant as being *highly segregating*, and produce a larger plot size for us to look at in more detail:

```
python /location/of/psQTL_post.py plot -d psqtl_analysis -i call \
    -f /location/of/genome.fasta -o plot_chr2.pdf \
    -p line histogram \
    --threshold 1.2 \
    --regions chr2 --width 15
```

Using this plot we can now (hypothetically) identify that the peak occurs specifically in the region from 3.8 Mbp to 4.2 Mbp. And we can see that our threshold value of `1.2` was overshooting the mark a little; let's adjust that down. Because we've zoomed in so much, we might also reduce the histogram `--bin` size down from the default of 100Kbp to 10Kbp. We will plot just the region we're interested in, and include all plot types relevant to the variant calls since this is our final product:

```
python /location/of/psQTL_post.py plot -d psqtl_analysis -i call \
    -f /location/of/genome.fasta -o plot_chr2_region.pdf \
    -p line scatter histogram gene \
    --threshold 1 --bin 10000 \
    --regions chr2:3800000-4200000
```

If we were looking at the `-i depth` data, maybe we'd want to include `-p coverage` as a visual indication of the sequencing depth within the region, too.

### Reporting
This is still to be implemented. It should let you tabulate reports to interrogate the data in the context of gene proximity to variants or deletions.

# How to cite
A publication is hopefully forthcoming which can be referred to when using this program. Until then, you can link to this repository.
