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
from q2_types.per_sample_sequences import SequencesWithQuality, Sequences

from q2_usearch._format import USEARCHFastaFmt, USEARCHFastQFmt, USEARCHFastaDirFmt, USEARCHFastQDirFmt
from q2_usearch._type import PooledSampleData, PooledSequencesWithQuality, PooledSequences
from qiime2.plugin import plugin


from q2_usearch import __version__
import q2_usearch

citations = Citations.load("citations.bib", package="q2_usearch")

plugin = Plugin(
    name="usearch",
    version=__version__,
    website="https://github.com/bokulich-lab/q2-plugin-name",
    package="q2_usearch",
    description="This is a template for building a new QIIME 2 plugin.",
    short_description="",
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
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality]},
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
