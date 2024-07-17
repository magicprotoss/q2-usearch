---
editor_options: 
  markdown: 
    wrap: 72
---

# q2-usearch

QIIME2 plug-in for [USEARCH](https://github.com/rcedgar/usearch12/)
integration

## Introduction

For years, [USEARCH](https://drive5.com/usearch/) had been the GOAT
program for many researchers to process amplicon sequencing data
(Including us🙈). We originally wrote this plug-in for internal use, now
with USEARCH's conversion to [open-source
software](https://x.com/RobertEdgarPhD/status/1802432815566553234), we
decided to publish this plug-in for the community to use. Here are the
commands/pipelines this plug-in integrates into
[QIIME2](https://qiime2.org/):

-   Denoise valid data into ZOTU table and ZOTUs using the uniose3
    algorithm
-   Cluster valid data into OTU table and OTUs at 97% identity
    threshould using the uparse algorithm
-   Denoise then cluster valid data into OTU table and OTUs at an
    user-defined identity threshould using both the uniose3 and the
    uclust algorithm
-   Perform paired-end read merging

\*The plug-in is still in early development, thus is subject to
interface changes

## Installation

Step 1: Clone this repository to your compute node

``` bash
git clone https://github.com/magicprotoss/q2-usearch .
```

Step 2: Activate the QIIME2 conda enviroment you wish to install to

``` bash
# conda activate qiime2-amplicon-2024.5
conda activate <replace-with-your-q2-conda-env-name>
```

In case your don't know your q2 env's name, please run the following
command

``` bash
conda env list | grep qiime2
```

The env's name should appear in your terminal

``` bash
# qiime2-amplicon-2024.5     /home/navi/miniconda3/envs/qiime2-amplicon-2024.5
```

Step 3: Change directory to the project folder and execute the following
command

``` bash
cd q2-usearch && python ./setup.py install
qiime dev refresh-cache
```

Step 4: Install [seqkit2](https://bioinf.shenwei.me/seqkit/) and
[usearch12](https://github.com/rcedgar/usearch12) using mamba/conda

``` bash
mamba install -c bioconda seqkit">=2.0.0" usearch
# conda install -c bioconda seqkit">=2.0.0" usearch
```

If every thing went smoothly, you should be seeing sth. like this
printed on your terminal

``` bash
usearch
# usearch v12.0 [b1d935b], 132Gb RAM, 24 cores
# (C) Copyright 2013-24 Robert C. Edgar.
# https://drive5.com/usearch
```

Step 5: Optional Clean Up

``` bash
cd .. && rm -rf q2-usearch
```

## Methods

### 1. denoise-no-primer-pooled

-   What's under the hood?

    ![](docs/denoise_no_primer_pooled.png)

-   [API document](docs/denoise_no_primer_pooled.md)

### 2. cluster-no-primer-pooled

-   What's under the hood?

    ![](docs/cluster_no_primer_pooled.png)

-   [API document](docs/cluster_no_primer_pooled.md)

### 3. denoise-then-cluster-no-primer-pooled

-   What's under the hood?

    ![](docs/denoise_then_cluster_no_primer_pooled.png)

-   [API document](docs/denoise_then_cluster_no_primer_pooled.md)

## Tutorials on zOTU Calling

### Process 'Valid data' from sequencing centers

1.  Identify and prepare your reads files

    Nowadays the common way a sequencing center sends you data is by
    providing you a link to a shared folder. Let's take a look at the
    data structure:

    ![](docs/folder_structure.png)

    What you need is the files in the "Valid data" folder. And before
    proceeding to the next step, we recommend you to backup the all
    content in said folder and then rename all the valid reads files by
    removing excess strings other than your Sample-ID. If you're on
    Windows, one of the easiest way to do it is to use the power rename
    tool in [MS
    PowerToys](https://learn.microsoft.com/en-us/windows/powertoys/):

    ![](docs/power_rename.png)

2.  Import reads files into a QIIME2 Artifact

    Let's start by activating the QIIME2 environment first

    ``` bash
    conda activate qiime2-amplicon-2024.5
    ```

    The easiest way to import your reads files into a QIIME2 Artifact is
    by using a manifest file, which serves the same purpose as a mapping
    file in QIIME1. Since we've already renamed our reads files with
    Sample-IDs, we can generate the manifest file easily using the file
    names. There's a couple of ways to do it, there are a lot of pure
    bash solutions on the internet, like [the
    one](https://library.qiime2.org/plugins/qiime2-manifest-metadata-generator/23/)
    in QIIME2's library, or use the utility script, which is our
    preferred method. This script requires one more dependency though,
    let's install it first.

    ``` bash
    # conda install xlsxwriter
    mamba install xlsxwriter
    ```

    With that out of the way, let's run the utility script

    ``` bash
    python generate_metadata.py --input_path <path-to-your-valid-data> --from_filename
    ```

    The script generates two files in the current directory. The first
    file 'manifest.tsv' contains the Sample-IDs and absolute file-paths
    for all of our reads, this will be required when importing our
    reads.

    ![](docs/manifest_example.png)

    The second file 'metadata.xlsx' is a pre-formatted QIIME2 metadata
    file which contains the Sample-IDs and an empty column
    'default-group'.

    ![](docs/metadata_example.png)

    If we want to use QIIME2 to perform downstream analysis, it's always
    worth spending a little more time filling the blanks Once we've
    uploaded the filled metadata file, we can convert it to tsv format
    for QIIME2 to use. See
    [here](https://docs.qiime2.org/2024.5/tutorials/metadata/#metadata-formatting-requirements)
    on how to prepare QIIME2 compatible metadata.

    ``` bash
    python generate_metadata.py --to_tsv
    ```

    Now everything's prepared, we can import our reads files into a
    QIIME2 Artifact. Note since we are using 'valid data', the input
    format and schematic type here are always
    'SingleEndFastqManifestPhred33V2' and
    'SampleData[SequenceWithQuality]', respectively, regardless of our
    sequencing strategy.

    ``` bash
    qiime tools import \
        --input-path manifest.tsv \
        --input-format 'SingleEndFastqManifestPhred33V2' \
        --type 'SampleData[SequencesWithQuality]' \
        --output-path fastq-seqs.qza
    ```

3.  Denoise reads into zOTUs

    Now you've got everything prepared, let's call the plug-in and
    finish the job.

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences fastq-seqs.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

    If you haven't received your sequencing data yet, We've prepared the
    data-set used in [(Dong, Guo et al.
    2021)](https://doi.org/10.1016/j.aquaculture.2020.736199) study for
    you to [try it
    out](https://zenodo.org/records/12655543?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQyMzZlOTlkLWFkZmEtNDM3NS1hMzUzLTFlNWMxZjJlM2Q1YSIsImRhdGEiOnt9LCJyYW5kb20iOiI2NzJmNDlkYTRmNTUwZDdlYzA4ZjllMmUyYmU2MTYzNCJ9.U9S4WNLhCgTYMz37j2QUatQ2F-HYzl9suSWlFr11nMdKslJNTraD5IuxvDmBgsAzUc9s8wTsm_9UJhTs6VZ9wA).

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences ddbj_dl.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

### General Use Cases

1.  Usage on Single-End Runs

    Let's use the ["Moving Pictures"
    tutorial](https://docs.qiime2.org/2024.5/tutorials/moving-pictures/)
    as a basis, please navigate to the '**Sequence quality control and
    feature table construction**' section. Besides **Option 1: DADA2**
    and **Option 2: Deblur**, we now have a third Option.

    #### **Option 3: UNOISE3**

    The [unoise3
    command](https://drive5.com/usearch/manual/cmd_unoise3.html) uses
    the [UNOISE
    algorithm](https://drive5.com/usearch/manual/unoise_algo.html) to
    perform denoising (error-correction) of amplicon reads. After which
    chimeras are removed by performing denovo chimera identification
    using an imporved version of [uchime2
    algorithm](https://drive5.com/usearch/manual/uchime2_algo.html). The
    `q2-usearch` plug-in have wrapped the valid-data processing pipeline
    described in [(Yan, Lin et al.
    2024)](https://link.springer.com/article/10.1186/s13717-023-00480-7)
    into the `denoise-no-primer-pooled` method.

    For single-end runs, it's mandatory to perform [global
    trimming](https://drive5.com/usearch/manual/global_trimming.html) to
    all your reads. Which means to trim all your reads to a fixed
    length, so that reads from the same template should have the same
    length, or they will be splited into separate zOTUs in later stage
    of the pipeline. For this specific data-set, global trimming is
    enforced by setting the parameter `--p-trunc-len n`, which truncates
    each sequence at position `n` and discards all sequences shorter
    than `n`.

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences demux.qza \
        --p-min-size 4 \
        --p-trunc-len 120 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

    ``` bash
    qiime metadata tabulate \
        --m-input-file stats-unoise3.qza \
        --o-visualization stats-unoise3.qzv
    ```

    If we'd like to continue the tutorial using this FeatureTable
    (opposed to the feature tables generated in *Option 1 and Option
    2*), run the following commands.

    ``` bash
    mv rep-seqs-unoise3.qza rep-seqs.qza
    mv table-unoise3.qza table.qza
    ```

2.  Usage on Paired-End Runs

    For paired-end reads, we need to merge them into joined reads prior
    to denoising. For this specific reason, we'll follow the ["Atacama
    soil
    microbiome"](https://docs.qiime2.org/2024.5/tutorials/atacama-soils/)
    tutorial. Please navigate to the section before running DADA2.

    Merging Paired End Reads

    The merge-pairs method was adopted from the plug-in q2-vsearch, with
    a few tweaks to the default parameters. For details on how to set
    parameters for reads merging, please check out the [usearch
    manual](https://rcedgar.github.io/usearch12_documentation/cmd_fastq_mergepairs.html).

    ``` bash
    qiime usearch merge-pairs \
        --i-demultiplexed-seqs demux.qza \
        --o-merged-sequences merged.qza \
        --o-unmerged-sequences unmerged.qza \
        --verbose
    ```

    zOTU Calling

    ``` bash
    qiime usearch denoise-no-primer-pooled \
        --i-demultiplexed-sequences merged.qza \
        --p-min-size 4 \
        --o-representative-sequences rep-seqs-unoise3.qza \
        --o-table table-unoise3.qza \
        --o-denoising-stats stats-unoise3.qza
    ```

## What's Planned for the future?

Several Methods and pipelines were planned for future releases:

Methods:

-   Classify FeatureData[Sequences] using sintax
-   Perform OTU cluster (uclust) by using raw PacBio CCS data and DADA2
    outputs as inputs.

Pipelines:

-   Perform merging(PE) ➡️ primer-removal ➡️ denoise/OTU-cluster on
    demultiplexed raw illumina data in a single pipeline

-   Find exact matches of FeatureData[Sequences] against a given
    database using global search then classify unmatched reads using
    sintax (similar to q2-feature-classifier's
    classify-hybrid-vsearch-sklearn and dada2's assignTaxonomy()
    followed by addSpecies())

Let me know if you have any questions😉

Happy QIIMEing 🎉🎉🎉
