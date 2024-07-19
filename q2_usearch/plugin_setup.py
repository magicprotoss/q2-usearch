# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import Citations, Plugin, Int, Range, Str, Choices, Float, Bool
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality, Sequences, JoinedSequencesWithQuality, PairedEndSequencesWithQuality

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
    description="This QIIME 2 plugin wraps USEARCH and supports: \n"
                "1. Paired-end reads merging \n"
                "2. Denoise illumina reads with unoise3 \n"
                "3. De novo OTU clustering with uparse \n"
                "4. De novo OTU clustering with uclust after denoise with unoise3 \n"
                "5. Taxonomy assignment with sintax",
    short_description="Plugin for amplicon analysis with USEARCH. ",
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
    description='This Method Pools All Samples Together and Extracts Biological Reads Using the Unoise3 Algorithm. \n' +
    'Non-Biological Sequence (i.e. Barcodes, Primers) MUST be REMOVED Prior to this step. \n' +
    'You MUST Also MERGE Your Reads If You are Using PAIRED-END Sequncing Protocol. \n' +
    "You Can Directly Use the 'Valid-Data' Provided by the Sequencing Center. \n" +
    'Vsearch was supported in early development but became deprecated for shipment.',
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
    inputs={
        'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={
        'demultiplexed_sequences': 'Quality screened, Adapter stripped, Joined(paired-end) sequences.'},
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
    description='This Method Pools All Samples Together and Cluster Them into 97% OTUs using the Uparse Algorithm. \n' +
    'Non-Biological Sequence (i.e. Barcodes, Primers) MUST be REMOVED Prior to this step. \n' +
    'You MUST Also MERGE Your Reads If You are Using PAIRED-END Sequncing Protocol. \n' +
    "You Can Directly Use the 'Valid-Data' Provided by the Sequencing Center. \n" +
    "Note: Nowadays 97% OTUs are Mostly Considered Mostly OBSELETE. ",
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
    inputs={
        'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={
        'demultiplexed_sequences': 'Quality screened, Adapter stripped, Joined(paired-end) sequences.'},
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence]),
             ('stats', SampleData[USEARCHStats])],
    output_descriptions={
        'table': 'The resulting feature table.',
        'representative_sequences': ('The resulting feature sequences. Each '
                                     'feature in the feature table will be '
                                     'represented by exactly one sequence. '),
        'stats': 'DataFrame containing statistics during each step of the pipeline.'
    }
)


plugin.methods.register_function(
    function=q2_usearch.denoise_then_cluster_no_primer_pooled,
    parameters={
        'trim_left': Int % Range(0, None),
        'trunc_len': Int % Range(0, None),
        'min_len': Int % Range(0, None),
        'max_ee': Float % Range(0.0, None),
        'perc_identity': Float % Range(0.0, 1.0),
        'min_size': Int % Range(1, None),
        'unoise_alpha': Float % Range(0.0, None),
        'n_threads': Int % Range(1, None) | Str % Choices(['auto']),
        'use_vsearch': Bool,
    },
    name="Pool, denoise then cluster valid-data.",
    description='This Method Pools All Samples Together, then Extracts Biological Reads Using the Unoise3 Algorithm, ' +
    "Finally OTUs are Generated by Clustering zOTUs Using the Uclust Algorithm at the user's defined identity threshold " +
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
        'perc_identity': 'The identity threshold used for OTU clustering. ',
        'min_size': ('The minimum abundance of input reads to be retained. '
                     'For higher sensivity, reducing minsize to 4 is reasonable. '
                     'Note: with smaller minsize, there tends to be more errors in low-abundance zotus. '),
        'unoise_alpha': 'See UNOISE2 paper for definition',
        'n_threads': ('The number of threads to use for computation. '
                      'If set to auto, the plug-in will use (all vcores - 3) present on the node.'),
        'use_vsearch': 'Use vsearch instead of usearch for computation . '
    },
    inputs={
        'demultiplexed_sequences': SampleData[SequencesWithQuality] | SampleData[JoinedSequencesWithQuality]},
    input_descriptions={
        'demultiplexed_sequences': 'Quality screened, Adapter stripped, Joined(paired-end) sequences.'},
    outputs=[('table', FeatureTable[Frequency]),
             ('representative_sequences', FeatureData[Sequence]),
             ('stats', SampleData[USEARCHStats])],
    output_descriptions={
        'table': 'The resulting feature table.',
        'representative_sequences': ('The resulting feature sequences. Each '
                                     'feature in the feature table will be '
                                     'represented by exactly one sequence. '),
        'stats': 'DataFrame containing statistics during each step of the pipeline.'
    }
)

plugin.methods.register_function(
    function=q2_usearch.sintax,
    parameters={
        'strand': Str % Choices('plus', 'both'),
        'threads': Int % Range(1, None) | Str % Choices(['auto']),
        'confidence': Float % Range(0.1, 1.0)
    },
    name="Rapidly classify reads by taxon using sintax.",
    description='This Method Classifies Query Sequences By Performing BootStrapped Kmer Search \n ' +
    'Against Refrence Using the Sintax Algorithm. \n' +
    'QIIME Style Reference Reads And Taxonomy Are Used As Inputs \n.' +
    'to Facilitate Tighter Intergration with the RESCRIPt plug-in and the QIIME2 EcoSystem . ',
    citations=[citations['edgar2016sintax']],
    parameter_descriptions={
        'strand': """Align against reference sequences in forward ("plus"), or both directions ("both").""",
        'threads': ('The number of threads to use for computation. '
                    'If set to auto, the plug-in will use (all vcores - 3) present on the node.'),
        'confidence': 'Confidence threshold for limiting taxonomic depth. '
    },
    inputs={
        'query': FeatureData[Sequence],
        'reference_reads': FeatureData[Sequence],
        'reference_taxonomy': FeatureData[Taxonomy]
    },
    input_descriptions={
        'query': 'Query sequences.',
        'reference_reads': 'Reference sequences.',
        'reference_taxonomy': 'Reference taxonomy labels.'
    },
    outputs=[('classification', FeatureData[Taxonomy])],
    output_descriptions={
        'classification': 'Taxonomy classifications of query sequences.'
    }
)

####################
# modified from q2-vsearch

plugin.methods.register_function(
    function=q2_usearch.merge_pairs,
    inputs={
        'demultiplexed_seqs': SampleData[PairedEndSequencesWithQuality]
    },
    parameters={
        'truncqual': Int % Range(0, None),
        'minlen': Int % Range(0, None),
        'allowmergestagger': Bool,
        'minovlen': Int % Range(5, None),
        'maxdiffs': Int % Range(0, None),
        'percent_identity': Int % Range(0, 100),
        'minmergelen': Int % Range(0, None),
        'maxmergelen': Int % Range(0, None),
        # 'preset': Int % None | Str % Choices(['double_reigon_short_overlap', 'double_reigon_long_overlap', 'signle_reigon_long_overlap']),
        'threads': Int % Range(1, None) | Str % Choices(['auto']),
    },
    outputs=[
        ('merged_sequences', SampleData[JoinedSequencesWithQuality]),
        ('unmerged_sequences', SampleData[PairedEndSequencesWithQuality])
    ],
    input_descriptions={
        'demultiplexed_seqs': ('The demultiplexed paired-end sequences to '
                               'be merged.'),
    },
    parameter_descriptions={
        'truncqual': ('Truncate sequences post merging at the first base with the '
                      'specified quality score value or lower.'),
        'minlen': ('Sequences shorter than minlen post merging and after min quality truncation are '
                   'discarded.'),
        'allowmergestagger': ('Allow merging of staggered read pairs, default is false. '
                              'Set it to True if you are using long overlap for a single reigon (i.e. 250PE for V4). '),
        'minovlen': ('Minimum length of the area of overlap between reads '
                     'during merging.'),
        'maxdiffs': ('Maximum number of mismatches in the area of overlap '
                     'during merging.'),
        'percent_identity': ('The percent identity threshold of overlaping reigon. '),
        'minmergelen': ('Minimum length of the merged read to be retained.'),
        'maxmergelen': ('Maximum length of the merged read to be retained.'),
        # 'preset': ('Predefined parameters for commonly used 16s protocols nowadays. '
        #            'Setting this will override parameters set above. '),
        'threads': ('The number of threads to use for computation. '
                    'Default is auto, which automaticlly scales up to 10 vcores if available. ')
    },
    output_descriptions={
        'merged_sequences': 'The merged sequences.',
        'unmerged_sequences': 'The unmerged paired-end reads.'
    },
    name='Merge paired-end reads.',
    description=('Merge paired-end sequence reads using usearch\'s '
                 'merge_pairs function. See the usearch documentation for '
                 'details on how paired-end merging is performed, and for '
                 'more information on the parameters to this method.'),
    citations=[citations['edgar2010usearch']]
)

# Register Usearch stats fmt
importlib.import_module('q2_usearch._transformer')
plugin.register_formats(USEARCHStatsFormat, USEARCHStatsDirFmt)
plugin.register_semantic_types(USEARCHStats)
plugin.register_semantic_type_to_format(
    SampleData[USEARCHStats], USEARCHStatsDirFmt)
