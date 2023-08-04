# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------


def lazy_process_clean_data(ctx,
                            clean_data,
                            n_threads="auto",
                            truncate_length_single=0,
                            ):
    
    pool_samples = ctx.get_action('usearch', 'pool_samples')
    quality_control = ctx.get_action('usearch', 'quality_control')
    dereplication = ctx.get_action('usearch', 'dereplication')
    denoise = ctx.get_action('usearch', 'unoise3')
    cluster_otu = ctx.get_action('usearch', 'uprase_otu')
    build_feature_table = ctx.get_action('usearch', 'build_feature_table')

    pooled_samples = pool_samples(demultiplexed_sequences= clean_data)
    quality_controlled_samples = quality_control(merged_sequences = pooled_samples.merged_sequences,
                                                truncate_single = truncate_length_single)
    dereplicated_sequences = dereplication(clean_reads = quality_controlled_samples.clean_reads)

    zotus, = denoise(dereplicated_reads = dereplicated_sequences.dereplicated_sequences, minsize = 4)
    otus, = cluster_otu(dereplicated_reads = dereplicated_sequences.dereplicated_sequences)
    zotu_table, = build_feature_table(merged_sequences = pooled_samples.merged_sequences,
                                      representative_sequences = zotus,
                                      rep_seqs_type = "ZOTU",
                                      n_threads = n_threads)
    otu_table, = build_feature_table(merged_sequences = pooled_samples.merged_sequences,
                                      representative_sequences = otus,
                                      rep_seqs_type = "OTU",
                                      n_threads = n_threads)


    return zotu_table, otu_table, zotus, otus
