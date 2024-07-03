# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._illumina_pipeline import py_to_cli_interface
import re

# prep test on pacbio dataset

## options

# 1. use -otutab with -strand both (fl 1492bp is very unlikely to have accidental rev hits)
# 2. write a custom R script to call dada2 to trim primers (very slow) (dps: do not need)
# 3. build the -fastq_orient command in the pipeline (dps: do not need)

# samples and features in featrue-tab must be present in both rep-seqs and demultiplexed_sequences
# if the feature-tab is sub-sampled, issue a warning to the user

def _construct_rep_seqs_with_size(working_dir,
                                  verbose = True):
    ### inputs
    # tab: pd.DataFrame
    # rep_seqs: pd.Series
    ###
    
    # check if feature-tab is and rep-seqs index match
    samples_in_featrue_tab = feature_tab.index
    feature_in_tab_match_seqs = feature_tab.T.index.isin(rep_seqs.index).all()
    feature_in_seqs_match_tab = rep_seqs.index.isin(feature_tab.T.index).all()
    if not (feature_in_tab_match_seqs and feature_in_seqs_match_tab):
        raise ValueError("Feature table and representative sequences do not match.")
    if not feature_in_seqs_match_tab:
        print("Warning: feature in rep-seqs missing in feature-tab, did you sub-sample the feature-tab?")
        print("Warning: poceeding anyway...")
    # step 1. sort according to decreasing abundance
    size_info_df = feature_tab.T.sum(axis = 1).sort_values(ascending = False).rename("size_info").to_frame()
    # check if the feature-tab contains features with 0 observations
    check_0_observation_feature = size_info_df.loc[size_info_df['size_info'] == 0].empty
    0_observation_feature_lst = size_info_df.loc[size_info_df['size_info'] == 0].index.tolist()
    
    if not check_0_observation_feature:
        print("Warning: feature with 0 observation found in feature-tab, did you sub-sample the feature-tab?")
        size_info_df = size_info_df.loc[size_info_df['size_info'] != 0]
        for feature_id in 0_observation_feature_lst:
            print("Warning: dropping feature:", feature_id, "due to: 0 observation")
    # step 2. left join and use index to check input

    size_info_df = size_info_df.merge(rep_seqs.rename("seqs"), how='inner', left_index=True, right_index=True)

    # step 3. iterrows and use skbio to write fasta
    rep_seqs_w_size_lst = [ skbio.DNA(str(row['seqs']).upper(), metadata={'id': str(index) + ';size=' + str(round(row['size_info'])) + ';'}) for index, row in size_info_df.iterrows()]
    
    stats_df = feature_tab.T.sum().rename("reads_mapped_to_features").to_frame()
    
    mapped_features = ~feature_tab.T.eq(0)
    
    mapped_features_count = mapped_features.sum().rename("feature_count")
    
    stats_df = stats_df.merge(mapped_features_count, how='inner', left_index=True, right_index=True)
    
    return rep_seqs_w_size_lst, samples_in_featrue_tab, stats_df

def _get_samples_in_fastq_input(demultiplexed_sequences_dirpath, verbose = True):
    input_manifest_df = pd.read_csv(os.path.join(
        demultiplexed_sequences_dirpath, 'MANIFEST'), index_col=0, comment='#')
    samples_in_fastq = input_manifest_df.index
    return 

def _match_features_in_feature_tab_and_fastq(samples_in_featrue_tab, samples_in_fastq):
    samples_in_feature_tab_all_in_fastq = samples_in_featrue_tab.isin(samples_in_fastq).all()
    samples_in_fastq_all_in_feature_tab = samples_in_fastq.isin(samples_in_featrue_tab).all()
    if not samples_in_feature_tab_all_in_fastq:
        samples_not_in_fastq = samples_in_featrue_tab[~samples_in_featrue_tab.isin(samples_in_fastq)].tolist()
        raise ValueError("Samples: " + str(samples_not_in_fastq) + "in feature-tab but missing in input fastq.")
    if not samples_in_fastq_all_in_feature_tab:
        print("Warning: samples in input fastq missing in feature-tab, did you sub-sample the feature-tab?")
        samples_not_in_feature_tab = samples_in_fastq[~samples_in_fastq.isin(samples_in_featrue_tab)].tolist()
        print("Warning: proceeding anyway...")
        print("Warning: dropping the following samples from downstream analysis:", str(samples_not_in_feature_tab))

def _uclust_cli(working_dir,
                rep_seqs_w_size_lst,
                identity = 0.99,
                strand = "plus",
                verbose = True):
    # step 1. run uclust (save parent-child info into an artifact in future versions)
    asv_fasta_fp = os.path.join(working_dir, 'asvs.fasta')
    otu_fasta_fp = os.path.join(working_dir, 'otus.fasta')
    with open(asv_fasta_fp, 'w') as f:
        for seq in rep_seqs_w_size_lst:
            seq.write(f, format='fasta')
    cmd = ["usearch",
           "-cluster_smallmem", asv_fasta_fp,
           "-id", str(identity),
           "-centroids", otu_fasta_fp,
           "-sortedby", "other"
           ]
    scilence = py_to_cli_interface(cmd)
    # just leave the feature ids alone

def _add_sample_identifiers_to_input_fq(demultiplexed_sequences_dirpath,
                                        working_dir,
                                        sample_id_lst,
                                        verbose = True):
                                            
    input_manifest_df = pd.read_csv(os.path.join(
        demultiplexed_sequences_dirpath, 'MANIFEST'), index_col=0, comment='#')
        
    input_manifest_df = input_manifest_df.loc[sample_id_lst]
        
    use_temp_sample_ids = False
        
    # check if all input sample_ids meet usearch sample identifier requirements
    if len(input_manifest_df) != len([ sample_id for sample_id in input_manifest_df.index.to_list() if re.match(r'^[a-zA-Z0-9_]+$', sample_id) ]):
        use_temp_sample_ids = True
        pipeout_denoise_stats_df = pd.DataFrame(
            index=input_manifest_df.index, columns=['input_reads', 'fixed_sample_id'])
    else:
        pipeout_denoise_stats_df = pd.DataFrame(
            index=input_manifest_df.index, columns=['input_reads'])
    unzipped_seqs_dirpath = os.path.join(working_dir, "unzipped_seqs")
    os.mkdir(unzipped_seqs_dirpath)
    relabed_seqs_dirpath = os.path.join(working_dir, "relabeled_seqs")
    os.mkdir(relabed_seqs_dirpath)
    seqs_stats_dfs_dirpath = os.path.join(working_dir, "seqs_stats_dfs")
    os.mkdir(seqs_stats_dfs_dirpath)
    pooled_seqs_fp = os.path.join(working_dir, "merged.fastq")
    
    # check avalibility of seqkit
    try:
        check_seqkit_version_stdout = subprocess.run(["seqkit"], stdout=subprocess.PIPE).stdout
        if "Version: 2." in check_seqkit_version_stdout.decode("utf-8"):
            use_seqkit = True
    except FileNotFoundError:
        use_seqkit = False
    
    if not use_seqkit:
        # just keep this part for debugging purpose
        # not sure if it's necessary
        # any one to use usearch to process og data?
        with open(os.path.join(demultiplexed_sequences_dirpath, 'metadata.yml')) as f:
            phred_offset = int(str(f.readlines()[0]).split(": ")[1])
        if phred_offset == 33:
            variant = "illumina1.8"
        elif phred_offset == 64:
            variant = "illumina1.3"
    
        if verbose:
            if use_temp_sample_ids:
                print("Relabeling input seqs with unique sample identifiers, this will take a while...\n")
            else:
                print("Adding sample-id to input seqs identifiers, this will take a while...\n")
    
        with open(pooled_seqs_fp, 'wt') as pooled_seqs_fh:
            sample_num = 1
            for index, row in input_manifest_df.iterrows():
                if use_temp_sample_ids:
                    sample_id = "S" + str(sample_num)
                else:
                    sample_id = str(index)
                fn = str(row['filename'])
                gzip_reader = gzip.open(
                    os.path.join(demultiplexed_sequences_dirpath, fn), 'rt')
    
                if verbose:
                    if use_temp_sample_ids:
                        print("Now Working on sample: " + str(index))
                        print("Temporarily Relabeling this sample to: " + sample_id)
                    else:
                        print("Now Working on sample: " + sample_id)
    
                uzipped_seq_fp = os.path.join(
                    unzipped_seqs_dirpath, sample_id + ".fastq")
                relabed_seq_fp = os.path.join(
                    relabed_seqs_dirpath, sample_id + ".fastq")

                with open(uzipped_seq_fp, 'wt') as f:
                    f.write(gzip_reader.read())
                    f.flush()

                # get input seqs count
                dna_seqs_gen = skbio.io.registry.read(
                    uzipped_seq_fp, format="fastq", verify=True, variant=variant)
                i = len([seq for seq in dna_seqs_gen])

                # build relab command
                cmd = ["usearch",
                       "-fastx_relabel", uzipped_seq_fp,
                       "-prefix", sample_id + ".",
                       "-fastqout", relabed_seq_fp
                       ]

                relab_log = py_to_cli_interface(cmd, False)

                os.remove(uzipped_seq_fp)

                # merge all relabed seqs into one file

                with open(relabed_seq_fp, 'rt') as seq:
                    shutil.copyfileobj(seq, pooled_seqs_fh)
                # del relabed seq to save space
                os.remove(relabed_seq_fp)
                
                # write input seqs count to stats_df
                pipeout_denoise_stats_df.loc[index, 'input_reads'] = i
                
                if use_temp_sample_ids:
                    # keep track of samples
                    pipeout_denoise_stats_df.loc[index, 'fixed_sample_id'] = sample_id
                    sample_num = sample_num + 1
                    
            # finally swap the index to fixed ids

            pipeout_denoise_stats_df.reset_index(inplace=True)
            pipeout_denoise_stats_df.rename(
                columns={'sample-id': 'original_sample_id'}, inplace=True)
            pipeout_denoise_stats_df.set_index('fixed_sample_id', inplace=True)
            pipeout_denoise_stats_df.index.name = "sample-id"
                
    
    # seqkit makes this go burrrrr......
    else:
        
        if verbose:
            print("Adding sample-id to input seqs identifiers, seqkit makes things go burrrrrrrrrr...\n")
        
        with open(pooled_seqs_fp, 'wb') as pooled_seqs_fh:
            
            sample_num = 1
            
            for index, row in input_manifest_df.iterrows():
                
                if use_temp_sample_ids:
                    sample_id = "S" + str(sample_num)
                else:
                    sample_id = str(index)
                fn = str(row['filename'])
                input_fp = os.path.join(demultiplexed_sequences_dirpath, fn)
                relabed_seq_fp = os.path.join(relabed_seqs_dirpath, sample_id + ".fastq")
                seqs_stats_df_fp = os.path.join(seqs_stats_dfs_dirpath, sample_id + "_stats.tsv")
                
                
                if verbose:
                    if use_temp_sample_ids:
                        print("Now Working on sample: " + str(index))
                        print("Temporarily Relabeling this sample to: " + sample_id)
                    else:
                        print("Now Working on sample: " + sample_id)
                
                # build relabel command
                relab_cmd = ["seqkit",
                       "replace",
                       "-i", input_fp,
                       "-o", relabed_seq_fp,
                       "-p", ".+",
                       "-r", sample_id + ".{nr}"]
                
                silence = py_to_cli_interface(relab_cmd, verbose=False)
                
                # append to pooled_seqs and then del
                with open(relabed_seq_fp, 'rb') as seq:
                    shutil.copyfileobj(seq, pooled_seqs_fh)
                
                os.remove(relabed_seq_fp)
                
                # build count seq num command
                stats_cmd = ["seqkit",
                       "stats",
                       "-T",
                       input_fp,
                       "-o", seqs_stats_df_fp]
                       
                silence = py_to_cli_interface(stats_cmd, verbose=False)
                
                # read stats df and append to pipeout_denoise_stats_df
                stats_df = pd.read_csv(seqs_stats_df_fp, sep="\t", index_col=None, header=0)
                pipeout_denoise_stats_df.loc[index, 'input_reads'] = int(stats_df.at[0, 'num_seqs'])
                

                
                if use_temp_sample_ids:
                    # keep track of samples
                    pipeout_denoise_stats_df.loc[index, 'fixed_sample_id'] = sample_id
                    sample_num = sample_num + 1
    
    # finally swap the index to fixed ids        
    if use_temp_sample_ids:

        pipeout_denoise_stats_df.reset_index(inplace=True)
        pipeout_denoise_stats_df.rename(
            columns={'sample-id': 'original_sample_id'}, inplace=True)
        pipeout_denoise_stats_df.set_index('fixed_sample_id', inplace=True)
        pipeout_denoise_stats_df.index.name = "sample-id"
    input_stats_df = pipeout_denoise_stats_df
    
    return input_stats_df

################################################################################
def _build_otu_table(working_dir,
                     identity = 0.99,
                     n_threads = "auto",
                     verbose = True):
    
    raw_reads_fp = os.path.join(working_dir, "merged.fastq")
    otus_fp = os.path.join(working_dir, "otus.fasta")
    tsv_otutab_fp = os.path.join(working_dir, "otu_tab.tsv")
    matched_otus_fp = os.path.join(working_dir, "matched_otus.fasta")
    log_fp = os.path.join(working_dir, "otutab.log")
    node_thread_count = os.cpu_count()
    
    # build otutab cmd
    cmd = ["usearch",
           "-otutab", raw_reads_fp,
           "-otus", otus_fp
           "-otutabout", tsv_otutab_fp,
           "-dbmatched", matched_otus_fp,
           "-id", str(identity),
           "-log", log_fp
           ]
           
    if threads != "auto":
        if threads > node_thread_count:
            if verbose:
                print("Number of threads specified higher than max available on node...")
                print("Setting threads to max available on node...")
            cmd += ["-threads", str(node_thread_count)]
        else:
            cmd += ["-threads", str(threads)]
    else:
        cmd += ["-threads", str(node_thread_count - 3)]

    # run command
    # we can do stats in another function
    otutab_log = py_to_cli_interface(cmd, verbose)
    
def _prep_output_for_artifact_api(working_dir, verbose = True):
    
    matched_otus_fp = os.path.join(working_dir, "matched_otus.fasta")
    tsv_otutab_fp = os.path.join(working_dir, "otu_tab.tsv")
    
    # 1. read the output otus into a pd.series obj
    otus_lst = [ seq for seq in skbio.io.read(matched_otus_fp, format='fasta') ]
    otus = pd.Series(otus_lst, index = [ seq.metadata['id'] for seq in otus_lst ])
    
    # 2. read the otutab into a pandas dataframe obj
    with open(tsv_otutab_fp) as fh:
        tab_df = biom.Table.from_tsv(fh, None, None, None).to_dataframe()
        
    # get reads count mapped to otus
    reads_mapped_to_features_df = tab_df.sum().to_frame().astype('int')
    reads_mapped_to_features_df.columns = ["reads_mapped_to_otus"]
    reads_mapped_to_features_df.index.name = "sample-id"
    
    # get otus mapped to each sample
    mapped_otus = ~tab_df.eq(0)
    otus_mapped_to_samples_se = mapped_otus.sum().rename("otu_count")
    
    stats_df = reads_mapped_to_features_df.merge(otus_mapped_to_samples_se, left_index=True, right_index=True, how = 'inner')

    # sort featrue tab
    tab_df['sum'] = tab_df.sum(axis=1)
    tab_df.sort_values(by="sum", ascending=False, inplace=True)
    tab_df.drop(columns=["sum"], inplace=True)
    ordered_zotu_ids = tab_df.index.tolist()

    otu_table = biom.Table(tab_df.values, tab_df.index, tab_df.columns)
    
    otus.reindex(tab_df.index)
                         
    return otu_table, otus, stats_df

def cluster_features_denovo_and_reconstruct_otu_table(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt, # add merge-pairs support
                                   identity: float = 0.99,
                                   strand: str = "plus"
                                   n_threads: str = "auto",
                                   ) -> (biom.Table, pd.Series, qiime2.Metadata):
    # usearch exclusive                                   
    verbose = True
    
###############################################

    demultiplexed_sequences_dirpath = str(demultiplexed_sequences)
    
    with tempfile.TemporaryDirectory() as usearch_wd:
        # step 1. convert rep-seqs and feature-tab into rep-seqs with size info
        rep_seqs_w_size_lst, samples_in_featrue_tab, feature_stats_df = _construct_rep_seqs_with_size(usearch_wd, verbose=verbose)
        samples_in_fastq = _get_samples_in_fastq_input(demultiplexed_sequences_dirpath, verbose=verbose)
        _match_features_in_feature_tab_and_fastq(samples_in_featrue_tab, samples_in_fastq)
        # step 2. perform uclust to get otus
        uclust_stats_str = _uclust_cli(usearch_wd, rep_seqs_w_size_lst, identity=identity, strand=strand, verbose=verbose)
        # step 3. add sample identifiers to input fastq
        input_stats_df = _add_sample_identifiers_to_input_fq(
            demultiplexed_sequences_dirpath, usearch_wd, samples_in_featrue_tab, threads = n_threads, verbose=verbose)
        # step 4. build otu table
        otutab_stat_df = _build_otu_table(usearch_wd, identity=identity,
                            threads=n_threads, verbose=verbose)
        # step 5. prep output for artifact api
        otu_table, otus, map_otu_stats_df = _prep_results_for_artifact_api(
            usearch_wd, verbose=verbose)
    
        # finally prep denoise stats df
        clustering_stats = input_stats_df.merge(
            feature_stats_df, how="left", left_index=True, right_index=True)
        clustering_stats['otu_to_asv_ratio'] = (
            clustering_stats['otu_count'] / clustering_stats['feature_count']) * 100
        clustering_stats = clustering_stats.merge(
            otutab_stat_df, how="left", left_index=True, right_index=True)
        clustering_stats['sequence_recovery_rate_imporvement'] = (
            clustering_stats['reads_mapped_to_otus'] / clustering_stats['input_reads']) * 100
        
        # if sample ids were swapped during the run, we need to swap the sample ids back
        if 'original_sample_id' in clustering_stats.columns:
            
            id_map_dict = clustering_stats['original_sample_id'].to_dict()
            clustering_stats.index = clustering_stats['original_sample_id']
            clustering_stats.index.name = 'sample-id'
            clustering_stats.drop(columns=['original_sample_id'], inplace=True)
            otu_table.update_ids(id_map_dict, axis='sample', inplace=True)
            
        clustering_stats.fillna(0, inplace=True)
        
        clustering_stats = qiime2.Metadata(clustering_stats)
    
    return otu_table, otus, clustering_stats
