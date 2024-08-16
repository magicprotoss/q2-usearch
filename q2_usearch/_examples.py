# ----------------------------------------------------------------------------
# Modified from q2-dada2
#
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


demux_single_url = \
    'https://data.qiime2.org/usage-examples/moving-pictures/demux.qza'

demux_paired_url = \
    'https://data.qiime2.org/usage-examples/atacama-soils/demux-full.qza'


def denoise_no_primer_pooled(use):
    demux_single = use.init_artifact_from_url('demux_single', demux_single_url)

    rep_seqs, table_unoise3, denoise_stats = use.action(
        use.UsageAction('usearch', 'denoise_no_primer_pooled'),
        use.UsageInputs(
            demultiplexed_seqs=demux_single,
            trunc_len=120
        ),
        use.UsageOutputNames(
            representative_sequences='zOTUs',
            table='zOTU_Table',
            stats='denoising_stats'
        )
    )

    rep_seqs.assert_output_type('FeatureData[Sequence]')
    table_unoise3.assert_output_type('FeatureTable[Frequency]')
    denoise_stats.assert_output_type('SampleData[USEARCHStats]')


def denoise_then_cluster_no_primer_pooled(use):
    demux_single = use.init_artifact_from_url('demux_single', demux_single_url)

    rep_seqs, table_uclust, denoise_stats = use.action(
        use.UsageAction('usearch', 'denoise_then_cluster_no_primer_pooled'),
        use.UsageInputs(
            demultiplexed_seqs=demux_paired,
            trunc_len=120,
            min_size=4
        ),
        use.UsageOutputNames(
            representative_sequences='99_OTUs',
            table='OTU_Table',
            stats='denoising_and_uclusting_stats'
        )
    )

    rep_seqs.assert_output_type('FeatureData[Sequence]')
    table_uclust.assert_output_type('FeatureTable[Frequency]')
    denoise_stats.assert_output_type('SampleData[USEARCHStats]')
    
    
def cluster_no_primer_pooled(use):
    demux_single = use.init_artifact_from_url('demux_single', demux_single_url)

    rep_seqs, table_uparse, denoise_stats = use.action(
        use.UsageAction('usearch', 'cluster_no_primer_pooled'),
        use.UsageInputs(
            demultiplexed_seqs=demux_paired,
            trunc_len=120
        ),
        use.UsageOutputNames(
            representative_sequences='97_OTUs',
            table='OTU_Table',
            stats='clustering_stats'
        )
    )

    rep_seqs.assert_output_type('FeatureData[Sequence]')
    table_uparse.assert_output_type('FeatureTable[Frequency]')
    stats.assert_output_type('SampleData[USEARCHStats]')
    
def merge_pairs(use):
    demux_paired = use.init_artifact_from_url('demux_paired', demux_paired_url)

    merged_sequences, unmerged_sequences = use.action(
        use.UsageAction('usearch', 'merge-pairs'),
        use.UsageInputs(
            demultiplexed_seqs=demux_paired
        ),
        use.UsageOutputNames(
            merged_sequences='merged_seqs',
            unmerged_sequences='unmerged_seqs'
        )
    )

    merged_sequences.assert_output_type('SampleData[JoinedSequencesWithQuality]')
    unmerged_sequences.assert_output_type('SampleData[PairedEndSequencesWithQuality]')
