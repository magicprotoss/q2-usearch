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
pipelines this plug-in (plan to) integrate into
[QIIME2](https://qiime2.org/):

-   Denoise valid data into ZOTU table and ZOTUs using the uniose3
    algorithm
-   Cluster valid data into OTU table and OTUs at 97% identity
    threshould using the uparse algorithm
-   Denoise then cluster valid data into OTU table and OTUs at an
    user-defined identity threshould using both the uniose3 and the
    uclust algorithm
-   Perform paired-end read merging
-   (to do) Remove primers on SampleData[JoinedSequencesWithQuality]
-   (to do) Classify FeatureData[Sequences] using sintax
-   (to do) Classify FeatureData[Sequences] using USEAECH's
    implementation of [RDP
    classifier](https://www.doi.org/10.1128/AEM.00062-07)
-   (to do) Perform Electronic PCR on SampleData[SequencesWithQuality]
    or SampleData[JoinedSequencesWithQuality]

Pipelines to add in future release:

-   Perform merging(PE)➡️primer-removal➡️denoise/cluster on
    demultiplexed raw NGS data in a single pipeline

-   Perform primer-removal➡️denoise(dada2)➡️cluster(uclust) on raw
    PacBio CCS data in a single pipeline

-   Find exact matches of FeatureData[Sequences] against a given
    database using global search then classify unmatched reads using
    sintax (similar to q2-feature-classifier's
    classify-hybrid-vsearch-sklearn and dada2's assignTaxonomy()
    followed by addSpecies())

\*The plug-in is still in early development, thus is subject to
interface changes

## Installation

Step 1: Clone this repository to your compute node

Step 2: Activate the QIIME2 conda enviroment you wish to install to

``` bash
# conda activate qiime2-amplicon-2024.2
conda activate <replace-with-your-q2-conda-env-name>
```

In case your don't know your q2 env's name, please run the following
command

``` bash
conda env list | grep qiime2
```

The env's name should appear in your terminal

``` bash
# qiime2-amplicon-2024.2     /home/navi/miniconda3/envs/qiime2-amplicon-2024.2
```

Step 3: Change directory to the project folder and execute the following
command

``` bash
cd q2-usearch && python ./setup.py install
qiime dev refresh-cache
cd ../ && rm -rf q2-usearch
```

Step 4: Install seqkit using mamba/conda

``` bash
mamba install seqkit">=2.0.0"
# conda install seqkit">=2.0.0"
```

Step 5: [Download
USEARCH](https://drive5.com/usearch/download.html)[(version11.0.667)](https://drive5.com/downloads/usearch11.0.667_i86linux32.gz),
unzip and rename the executable to "usearch"

``` bash
wget https://github.com/rcedgar/usearch_old_binaries/blob/main/bin/usearch11.0.667_i86linux64 && mv usearch11.0.667_i86linux64 usearch
```

Step 6: Move the executable to your system's executable binary path and
add privilege to it

``` bash
# Install to the user's bin
# sudo mv usearch /usr/bin && sudo chmod +x /usr/bin/usearch
# Install system-wide
# sudo mv usearch /bin && sudo chmod +x /bin/usearch
# Install to current qiime2's conda env
mv usearch $(whereis qiime | sed 's/qiime//g')
chmod +x $(whereis qiime | sed 's/qiime//g')"usearch"
```

If every thing went smoothly, you should be seeing sth. like this
printed on your terminal

``` bash
usearch --version
# usearch v11.0.667_i86linux64
```

## Tutorials

### Before you start

-   Working environment set-up

-   To cluster or not to cluster🧐

### Step-by-step guides on common sequencing protocols

-   Valid data from illumina sequencing

-   Single-end raw data from illumina sequencing

-   Paired-end raw data from illumina sequencing

-   Raw pacbio ccs data

-   Valid pacbio ccs data

Let me know if you have any questions😉

Happy QIIMEing 🎉🎉🎉
