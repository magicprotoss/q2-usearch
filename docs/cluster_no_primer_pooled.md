## cluster-no-primer-pooled: Cluster valid data into OTU table and OTUs at 97% identity using the uparse algorithm

### Command Line Interface

``` bash
Usage: qiime usearch cluster-no-primer-pooled [OPTIONS]

  This Method Pools All Samples Together and Cluster Them into 97% OTUs using
  the Uparse Algorithm. Non-Biological Sequence (i.e. Barcodes, Primers) MUST
  be REMOVED Prior to this step. You MUST Also MERGE Your Reads If You are
  Using PAIRED-END Sequncing Protocol You Can Directly Use the 'Valid-Data'
  Provided by the Sequencing Center. Note: Nowadays 97% OTUs are Mostly Considered
  Mostly OBSELETE. 

Inputs:
  --i-demultiplexed-sequences ARTIFACT SampleData[SequencesWithQuality] |
    SampleData[JoinedSequencesWithQuality]
                         Quality screened, Adapter stripped,
                         Joined(paired-end) sequences.              [required]
Parameters:
  --p-trim-left INTEGER  Position at which sequences should be trimmed due to
    Range(0, None)       low quality. This trims the 5' end of the of the
                         input sequences.                         [default: 0]
  --p-trunc-len INTEGER  Position at which sequences should be truncated due
    Range(0, None)       to decrease in quality. This truncates the 3' end of
                         the of the input sequences. Reads that
                         are shorter than this value will be discarded. If 0
                         is provided, no truncation or length filtering will
                         be performed                             [default: 0]
  --p-min-len INTEGER    Reads with less length than this number value will
    Range(0, None)       be discarded.                           [default: 50]
  --p-max-ee NUMBER      Reads with number of expected errors higher than
    Range(0.0, None)     this value will be discarded.          [default: 1.0]
  --p-n-threads VALUE Int % Range(1, None) | Str % Choices('auto')
                         The number of threads to use for computation. If set
                         to auto, the plug-in will use (all vcores - 3)
                         present on the node.                [default: 'auto']
  --p-min-size INTEGER   The minimum abundance of a given OTU to be retained.
    Range(1, None)       Default is 2, which means unique reads are discarded.
                                                                  [default: 2]
Outputs:
  --o-table ARTIFACT FeatureTable[Frequency]
                         The resulting feature table.               [required]
  --o-representative-sequences ARTIFACT FeatureData[Sequence]
                         The resulting feature sequences. Each feature in the
                         feature table will be represented by exactly one
                         sequence.                                  [required]
  --o-denoising-stats ARTIFACT SampleData[USEARCHStats]
                         DataFrame containing statistics during each step of
                         the pipeline.                              [required]
Miscellaneous:
  --output-dir PATH      Output unspecified results to a directory
  --verbose / --quiet    Display verbose output to stdout and/or stderr
                         during execution of this action. Or silence output if
                         execution is successful (silence is golden).
  --example-data PATH    Write example data and exit.
  --citations            Show citations and exit.
  --help                 Show this message and exit.
```

### Artifact API

#### Import

``` python
from qiime2.plugins.usearch.methods import cluster_no_primer_pooled
```

#### Example Usage

``` python
otu_tab, otus, stats, = cluster_no_primer_pooled(<your_valid_data_artifact>)
```

#### **Docstring**

```         
Cluster valid data into 97% OTUs

  This Method Pools All Samples Together and Cluster Them into 97% OTUs using
  the Uparse Algorithm. Non-Biological Sequence (i.e. Barcodes, Primers) MUST
  be REMOVED Prior to this step. You MUST Also MERGE Your Reads If You are
  Using PAIRED-END Sequncing Protocol You Can Directly Use the 'Valid-Data'
  Provided by the Sequencing Center. Note: Nowadays 97% OTUs are Mostly Considered
  Mostly OBSELETE. 

Parameters
----------
demultiplexed_seqs : SampleData[SequencesWithQuality]
    The single-end demultiplexed PacBio CCS sequences to be denoised.
trim_left : Int, optional
    Position at which sequences should be trimmed due to low quality. This
    trims the 5' end of the of the input sequences.
trunc_len : Int, optional
    Position at which sequences should be truncated due to decrease in
    quality. This truncates the 3' end of the of the input sequences. Reads that
    are shorter than this value will be discarded. If 0 is provided, no
    truncation or length filtering will be performed.
min_len : Int, optional
    Minimal length for reads to be retained. Reads shorter than this value will be excluded     for subsequent denoising.
max_ee : Float, optional
    Reads with number of expected errors higher than this value will be
    discarded.
n_threads : Threads, optional
    The number of threads to use for multithreaded processing. If "auto" is
    provided, the value of ( system core count - 3 ) will be used.
min_size : Int, optional
    The minimum abundance of a given OTU to be retained.
    Default is 2, which means unique reads are discarded.

Returns
-------
table : FeatureTable[Frequency]
    The resulting feature table.
representative_sequences : FeatureData[Sequence]
    The resulting feature sequences. Each feature in the feature table will
    be represented by exactly one sequence.
stats : SampleData[USEARCHStats]
    DataFrame containing statistics during each step of the pipeline
```
