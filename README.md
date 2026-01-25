# Table of Contents
- [Getting started](#getting-started)
- [Installation](#installation)
- [How to use psQTL](#how-to-use-psQTL)
- [How to cite](#how-to-cite)

# Getting started
```
# Obtain psQTL
conda create -n psqtl bioconda::psqtl
conda activate psqtl

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
python /location/of/psQTL_proc.py ed -d /location/of/your/working/directory \
    -i call depth --threads 12
python /location/of/psQTL_proc.py splsda -d /location/of/your/working/directory \
    -i call depth --threads 12

# 'Post'-processing steps including data plotting and report tabulation
python /mnt/c/git/psQTL/psQTL_post.py plot -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    -i call depth \
    -m ed splsda \
    --ed alleles \ # or inheritance or genotypes
    -o /location/to/write/output.pdf \ # or .png or .svg file
    -p line scatter coverage genes \
    --annotation location/of/genome.gff3
python /mnt/c/git/psQTL/psQTL_post.py report -d /location/of/your/working/directory \
    -f /location/of/genome.fasta \
    -i call depth \
    -m ed splsda \
    -o /location/to/write/output.tsv \ # or .csv
    -a location/of/genome.gff3
```

# Installation
It is recommended that you set up an Anaconda or Miniconda environment to run psQTL. A conda package is provided through bioconda to simplify the installation of psQTL and its prerequisites. To ensure package version compatibility, obtain psQTL in its own environment by doing:

```
conda create -n psqtl bioconda::psqtl
```

See the [Installing psQTL wiki page](https://github.com/zkstewart/psQTL/wiki/Installing-psQTL) for more details.

# How to use psQTL
The psQTL pipeline is intended to proceed from **prep**aration, to **proc**essing, to **post**-processing as accomplished by the correspondingly named Python script files.

All scripts will return help information on the command line in order to assist you in running them like:

```
python /location/of/psQTL_prep.py -h
python /location/of/psQTL_prep.py initialise -h
python /location/of/psQTL_proc.py ed -h
python /location/of/psQTL_post.py plot -h
... etc ...
```

All parameters with a single dash e.g., `-d /location/to/run/the/analysis` need to be provided each time you call the program. Other parameters specified with double dashes are optional, with many having default values that should work in the majority of analyses.

Some double-dashed parameters which indicate file locations e.g., `--meta metadata.tsv` are remembered or *cached* by psQTL within the analysis directory and only need to be provided once since the cached value is retrieved on subsequent use of psQTL ***in the analysis folder***. In other words, each analysis directory has its own cached values, and parameters won't carry over between separate analyses.

See the [Using psQTL wiki page](https://github.com/zkstewart/psQTL/wiki/Using-psQTL) for further details.

# How to cite
A publication is hopefully forthcoming which can be referred to when using this program. Until then, you can link to this repository.
