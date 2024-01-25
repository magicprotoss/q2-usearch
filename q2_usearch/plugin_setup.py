# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------
import importlib

from qiime2.plugin import Citations, Plugin, Int, Range, Str, Choices, Float, Bool
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality, Sequences, JoinedSequencesWithQuality

from q2_usearch._format import USEARCHFastaFmt, USEARCHFastQFmt, USEARCHFastaDirFmt, USEARCHFastQDirFmt
from q2_usearch._type import PooledSampleData, PooledSequencesWithQuality, PooledSequences
from qiime2.plugin import plugin


from q2_usearch import __version__
import q2_usearch

# Register Usearch stats fmt
from q2_usearch._format import USEARCHStats, USEARCHStatsFormat, USEARCHStatsDirFmt

citations = Citations.load("citations.bib", package="q2_usearch")

plugin = Plugin(
    name="usearch",
    version=__version__,
    website="https://github.com/magicprotoss/q2-usearch",
    package="q2_usearch",
    description="This QIIME 2 plugin wraps USEARCH and supports: "
                "1. Paired-end reads merging "
                "2. Denoise illumina reads with unoise3 "
                "3. De novo OTU cluster with uparse3 "
                "4. De novo ZOTU cluster with uclust after denoise with unoise3 "
                "5. Classify representative sequences with SINTAX "
                "6. Classify representative sequences with USEARCH's RDP classifier "
                "7. Perform ePCR on raw reads for meta analysis",
    short_description="Plugin for amplicon analysis with USEARCH. ",
)

plugin.methods.register_function(
    function=q2_usearch.unoise3,
    parameters={
        'minsize': Int % Range(1, None),
        'alpha': Float,
    },
    name="ZOTU picking using unoise3.",
    description="Perform ZOTU picking on quality screened, dereplicated sequences using unoise3",
    citations=[citations['edgar2016unoise2']],
    parameter_descriptions={
        'minsize': 'Minimum frequency threshold for a given ZOTU to be retained. '
        'The defalut value is 8, for further information. '
        'Visit https://drive5.com/usearch/manual/cmd_unoise3.html ',
        'alpha': 'The alpha parameter. '
        'See defination in paper https://doi.org/10.1101/081257 '
        'Do not modify this value unless you know what you are doing.',
    },
    inputs={'dereplicated_reads': PooledSampleData[PooledSequences]},
    input_descriptions={'dereplicated_reads': 'Dereplicated sequences.'},
    outputs=[('zotus', FeatureData[Sequence])],
    output_descriptions={'zotus': 'Representative sequences without chimeras.'}
)


plugin.methods.register_function(
    function=q2_usearch.uprase_otu,
    parameters={
        'minsize': Int % Range(1, None),
    },
    name="De novo 97% OTU picking using uprase3-otu.",
    description='Perform De Novo 97% OTU picking on quality screened, dereplicated sequences using uprase3-otu. ' +
    'WARNING: OTU picking on amlicons had been considered obsolete shortly after ESVs came out. ' +
    'ESV(Extracted Sequence Variant)s include zOTUs (usearch), sOTUs (Deblur), and ASVs (DADA2). ' +
    'We recommend to use this method only as a reference, not for publication. ',
    citations=[citations['edgar2013uparse']],
    parameter_descriptions={
        'minsize': 'Minimum frequency threshold for a given OTU to be retained. '
        'The defalut value is 2, for further information. '
        'Visit https://drive5.com/usearch/manual/cmd_cluster_otus.html '
    },
    inputs={'dereplicated_reads': PooledSampleData[PooledSequences]},
    input_descriptions={'dereplicated_reads': 'Dereplicated sequences.'},
    outputs=[('otus', FeatureData[Sequence])],
    output_descriptions={'otus': 'Representative sequences without chimeras.'}
)


plugin.methods.register_function(
    function=q2_usearch.pool_samples,
    parameters={
        'keep_annotations': Bool,
    },
    name="Pool all samples together and dereplicate them for otu and zotu picking with usearch.",
    description='Pool all samples together and dereplicate them into a single fasta file. ' +
    'I cannot make pooling samples and dereplication two separate steps with limited experience ' +
    'with qiime2 architecture, maybe more versitile options will be avaliable in the future :( ',
    citations=[citations['edgar2010usearch']],
    parameter_descriptions={
        'keep_annotations': 'Preseve existing annotations. '
        'This is turned off by default '
        'As you do not need seqs annotations to run usearch. '
    },
    inputs={
        'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={
        'demultiplexed_sequences': 'Quality screened, Adapter stripped, Length trimmed(single-end), Joined(paired-end) sequences.'},
    outputs=[('merged_sequences', PooledSampleData[PooledSequencesWithQuality])],
    output_descriptions={
        'merged_sequences': 'All samples merged into a single fastq file.'}
)

plugin.methods.register_function(
    function=q2_usearch.quality_control,
    parameters={
        'max_error_rate': Float,
        'truncate_single': Int,
        'min_sequence_len': Int,
    },
    name="Perform max error rate based quality control on samples.",
    description='Here is what Edgar says: ' +
    'It is important to use the USEARCH fastq_mergepairs command for read merging, ' +
    'because according to my tests, other read mergers do not correctly calculate the posterior Q scores, ' +
    'including PANDAseq, COPE, PEAR, fastq-join and FLASH. ' +
    'Some of these mergers also have high rates of false positive and incorrect merges,' +
    'including especially PANDAseq and the make.contigs() command in mothur.',
    citations=[citations['edgar2015errorfilter']],
    parameter_descriptions={
        'max_error_rate': 'Discard reads with > E total expected errors for all bases in the read '
        'after any truncation options have been applied. ',
        'truncate_single': 'Truncate sequences at the $L th base. If the sequence is shorter than L, discard. '
        'Do not apply this to merged paired-end reads. ',
        'min_sequence_len': 'Discard sequences with < $L letters. '
    },
    inputs={'merged_sequences': PooledSampleData[PooledSequencesWithQuality]},
    input_descriptions={'merged_sequences': 'Merged samples.'},
    outputs=[('clean_reads', PooledSampleData[PooledSequences])],
    output_descriptions={
        'clean_reads': 'Merged samples after quality control.'}
)

plugin.methods.register_function(
    function=q2_usearch.dereplication,
    parameters={
        'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
    },
    name="Dereplicate merged sequences for otu and zotu picking with usearch.",
    description='zaixiele zaixiele . ' +
    'zaixiele zaixiele ' +
    'zaixiele zaixiele ',
    citations=[citations['edgar2010usearch']],
    parameter_descriptions={
        'n_threads': 'zaixiele zaixiele '
        'zaixiele zaixiele '
        'zaixiele zaixiele. '
    },
    inputs={'clean_reads': PooledSampleData[PooledSequences]},
    input_descriptions={'clean_reads': 'Quality checked reads.'},
    outputs=[('dereplicated_sequences', PooledSampleData[PooledSequences])],
    output_descriptions={
        'dereplicated_sequences': 'Unique sequences in all samples.'}
)

plugin.methods.register_function(
    function=q2_usearch.build_feature_table,
    parameters={
        'rep_seqs_type': Str % Choices(['ZOTU', 'OTU']),
        'n_threads': Int % Range(1, None) | Str % Choices(['auto'])
    },
    name="Build feature table.",
    description='Pool all samples together and dereplicate them into a single fasta file. ' +
    'I cannot make pooling samples and dereplication two separate steps with limited experience ' +
    'with qiime2 architecture, maybe more versitile options will be avaliable in the future :( ',
    citations=[citations['edgar2010usearch']],
    parameter_descriptions={
        'n_threads': 'Number of threads to use when searching unique sequences. '
        'It is better to set it to auto, this step is very fast. '
        'You would probably nerver need more than 10 threads '
                     '(Which is the default cap in usearch). '
    },
    inputs={'merged_sequences': PooledSampleData[PooledSequencesWithQuality],
            'representative_sequences': FeatureData[Sequence]},
    input_descriptions={'merged_sequences': 'Quality screened, Adapter stripped, Length trimmed(single-end), Joined(paired-end) sequences.',
                        'representative_sequences': 'Unique sequences in all samples.'},
    outputs=[('feature_table', FeatureTable[Frequency])],
    output_descriptions={'feature_table': 'The resulting feature table.'}
)

plugin.methods.register_function(
    function=q2_usearch.sintax,
    parameters={
        'threads': Int % Range(1, None) | Str % Choices(['auto']),
        'confidence': Float % Range(0, 1, inclusive_end=True),
        'drop_species': Bool,
    },
    name="Taxonomy classification via sintax.",
    description='zaixiele zaixidele. ' +
    'zaixiele zaixidele ' +
    'zaixiele zaixidele ',
    citations=[citations['edgar2010usearch']],
    parameter_descriptions={
        'threads': 'Number of threads to use when searching unique sequences. '
        'If set to auto, the plug-in will use all sys vcores - 1. ',
        'confidence': 'Minimum confidence score for a taxonomy assignment to be retained. ',
        'drop_species': 'Drop species level taxonomy assignment (Recommend by Edgar). '
    },
    inputs={'query': FeatureData[Sequence],
            'reference_reads': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    input_descriptions={'query': 'zaixiele zaixiele .',
                        'reference_reads': 'zaixiele zaixiele.',
                        'taxonomy': 'zaixiele zaixiele.'},
    outputs=[('classification', FeatureData[Taxonomy])],
    output_descriptions={
        'classification': 'Taxonomy Classification via sintax.'}
)

####################
# revised methods
plugin.methods.register_function(
    function=q2_usearch.denoise_no_primer_pooled,
    parameters={
        'trim_left': Int % Range(0, None),
        'trunc_len': Int % Range(0, None),
        'min_len': Int % Range(0, None),
        'max_ee': Float % Range(0.0, None),
        'min_size': Int % Range(1, None),
        'unoise_alpha': Float % Range(0.0, None),
        'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
        'use_vsearch': Bool,
    },
    name="Pool and denoise valid-data.",
    description='This Method Pools All Samples Together and Extracts Biological Reads Using the Unoise3 Algorithm ' +
    'Non-Biological Sequence (i.e. Barcodes, Primers) MUST be REMOVED Prior to this step ' +
    'You MUST Also MERGE Your Reads If You are Using PAIRED-END Sequncing Protocol ' +
    "You Can Directly Use the 'Valid-Data' Provided by the Sequencing Center " +
    'Using Vsearch as a drop-in Replcacement is supported But with some CAVEATS, see https://github/xxx for details. ',
    citations=[citations['edgar2016unoise2']],
    parameter_descriptions={
        'trim_left': ("Position at which sequences should be trimmed due to low quality. "
                      "This trims the 5' end of the of the input sequences, "
                      "which will be the bases that were sequenced in the first cycles. "),
        'trunc_len': ("Position at which sequences should be truncated due to decrease in quality. "
                      "This truncates the 3' end of the of the input sequences, "
                      "which will be the bases that were sequenced in the last cycles. "
                      "Reads that are shorter than this value will be discarded. "
                      "If 0 is provided, no truncation or length filtering will be performed"),
        'min_len': 'Reads with less length than this number value will be discarded. ',
        'max_ee': 'Reads with number of expected errors higher than this value will be discarded. ',
        'min_size': ('The minimum abundance of input reads to be retained. '
                    'For higher sensivity, reducing minsize to 4 is reasonable. '
                    'Note: with smaller minsize, there tends to be more errors in low-abundance zotus. '),
        'unoise_alpha': 'See UNOISE2 paper for definition',
        'n_threads': ('The number of threads to use for computation. '
                      'If set to auto, the plug-in will use (all vcores - 3) present on the node.'),
        'use_vsearch': 'Use vsearch instead of usearch for computation . '
    },
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={'demultiplexed_sequences': 'Quality screened, Adapter stripped, Joined(paired-end) sequences.'},
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence]),
             ('denoising_stats', SampleData[USEARCHStats])],
    output_descriptions={
        'table': 'The resulting feature table.',
        'representative_sequences': ('The resulting feature sequences. Each '
                                     'feature in the feature table will be '
                                     'represented by exactly one sequence. '),
        'denoising_stats': 'DataFrame containing statistics during each step of the pipeline.'
        }
)


plugin.methods.register_function(
    function=q2_usearch.cluster_no_primer_pooled,
    parameters={
        'trim_left': Int % Range(0, None),
        'trunc_len': Int % Range(0, None),
        'min_len': Int % Range(0, None),
        'max_ee': Float % Range(0.0, None),
        'min_size': Int % Range(1, None),
        'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
    },
    name="Pool and cluster valid-data at 97% identity.",
    description='This Method Pools All Samples Together and Cluster Them into 97% OTUs using the Uparse Algorithm ' +
    'Non-Biological Sequence (i.e. Barcodes, Primers) MUST be REMOVED Prior to this step ' +
    'You MUST Also MERGE Your Reads If You are Using PAIRED-END Sequncing Protocol ' +
    "You Can Directly Use the 'Valid-Data' Provided by the Sequencing Center " +
    "Since Nowadays 97% OTUs are Considered Mostly OBSELETE and Uparse is Usearch Exclusive," +
    "Using Vsearch as a drop-in Replcacement is Not Supported Here. ",
    citations=[citations['edgar2013uparse']],
    parameter_descriptions={
        'trim_left': ("Position at which sequences should be trimmed due to low quality. "
                      "This trims the 5' end of the of the input sequences, "
                      "which will be the bases that were sequenced in the first cycles. "),
        'trunc_len': ("Position at which sequences should be truncated due to decrease in quality. "
                      "This truncates the 3' end of the of the input sequences, "
                      "which will be the bases that were sequenced in the last cycles. "
                      "Reads that are shorter than this value will be discarded. "
                      "If 0 is provided, no truncation or length filtering will be performed"),
        'min_len': 'Reads with less length than this number value will be discarded. ',
        'max_ee': 'Reads with number of expected errors higher than this value will be discarded. ',
        'min_size': ('The minimum abundance of a given OTU to be retained. '
                    'Default is 2, which means unique reads are discarded. '),
        'n_threads': ('The number of threads to use for computation. '
                      'If set to auto, the plug-in will use (all vcores - 3) present on the node.'),
    },
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={'demultiplexed_sequences': 'Quality screened, Adapter stripped, Joined(paired-end) sequences.'},
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence]),
             ('denoising_stats', SampleData[USEARCHStats])],
    output_descriptions={
        'table': 'The resulting feature table.',
        'representative_sequences': ('The resulting feature sequences. Each '
                                     'feature in the feature table will be '
                                     'represented by exactly one sequence. '),
        'denoising_stats': 'DataFrame containing statistics during each step of the pipeline.'
        }
)

####################

# reg pipelines
plugin.pipelines.register_function(
    function=q2_usearch.lazy_process_clean_data,
    inputs={
        'clean_data': SampleData[SequencesWithQuality],
    },
    parameters={
        'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
        'truncate_length_single': Int % Range(0, 350)
    },
    outputs=[
        ('zotu_table', FeatureTable[Frequency]),
        ('otu_table', FeatureTable[Frequency]),
        ('zotus', FeatureData[Sequence]),
        ('otus', FeatureData[Sequence]),
    ],
    input_descriptions={
        'clean_data': 'Demultiplexed, joined paired-end reads.'
    },
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use `all` to automatically use '
                     'all available cores -1. This value is used when constructing '
                     'feature tables.',
        'truncate_length_single': 'Truncate reads to a specified length. '
                                  'Only sinlge-end reads need to be truncated. '
                                  'Do NOT use this option for PAIRED-END reads. ',
    },
    output_descriptions={
        'zotu_table': 'The corrosponding ZOTU table.',
        'otu_table': 'The corrosponding OTU table.',
        'zotus': 'Extracted sequncing variant by the unoise3 algorithm.',
        'otus': 'Operational Taxonomy Units clustered (de-novo) at 97 percent indentity.',
    },
    name='Lazy Pipeline for processing Illumia Sequencing Data.',
    description=('zaixiele zaixiele '
                 'zaixiele zaixiele '
                 )
)


importlib.import_module('q2_usearch._transformer')


# Due to sample indentifier incompatibility, MUST Create usearch file formats, MUST reg them to correct shcematic types...
# Reg formats
plugin.register_formats(USEARCHFastaFmt, USEARCHFastQFmt,
                        USEARCHFastaDirFmt, USEARCHFastQDirFmt)

plugin.register_semantic_types(
    PooledSampleData, PooledSequencesWithQuality, PooledSequences)


# Reg formats to semantic types
plugin.register_artifact_class(
    PooledSampleData[PooledSequencesWithQuality],
    directory_format=USEARCHFastQDirFmt,
    description=("All samples merged into a single sequences.fastq file "
                 "Formatted with USEARCH compatibale Sample identifier.")
)

plugin.register_artifact_class(
    PooledSampleData[PooledSequences],
    directory_format=USEARCHFastaDirFmt,
    description=("All samples merged into a single sequences.fasta file "
                 "Formatted with USEARCH compatibale Sample identifier.")
)

# Register Usearch stats fmt
plugin.register_formats(USEARCHStatsFormat, USEARCHStatsDirFmt)
plugin.register_semantic_types(USEARCHStats)
plugin.register_semantic_type_to_format(
    SampleData[USEARCHStats], USEARCHStatsDirFmt)
