# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import qiime2
import pandas as pd
import numpy as np
import subprocess
import tempfile
import skbio
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt
from glob import glob
import shutil
import gzip
import os
import re
import time
import hashlib
import biom


def py_to_cli_interface(cmd, verbose=True):
    try:
        cmd_log_index = cmd.index("-log") + 1
        log_fp = cmd[cmd_log_index]
    except ValueError:
        log_fp = ""
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
        time.sleep(1)
        subprocess.run(cmd, check=True)
    else:
        time.sleep(1)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    if os.path.exists(log_fp):
        with open(log_fp, "rt") as f:
            log_lines_lst = [line.replace("\n", "").strip()
                             for line in f.readlines() if line != "\n"]
    else:
        log_lines_lst = ["no_log_file was produced during cmd run"]

    return log_lines_lst


# Pool All Samples into a single fastq

def _pool_samples(demultiplexed_sequences_dirpath, working_dir, keep_annotations: bool = False, use_vsearch: bool = False, debug = False, verbose: bool = True):
    input_manifest_df = pd.read_csv(os.path.join(
        demultiplexed_sequences_dirpath, 'MANIFEST'), index_col=0, comment='#')
        
    use_temp_sample_ids = False
        
    # check if all input sample_ids meet usearch sample identifier requirements
    if len(input_manifest_df) != len([ sample_id for sample_id in input_manifest_df.index.to_list() if re.match(r'^[a-zA-Z0-9_]+$', sample_id) ]):
        use_temp_sample_ids = True
        pipeout_denoise_stats_df = pd.DataFrame(
            index=input_manifest_df.index, columns=['prior_to_maxee_filt', 'fixed_sample_id'])
    else:
        pipeout_denoise_stats_df = pd.DataFrame(
            index=input_manifest_df.index, columns=['prior_to_maxee_filt'])
    unzipped_seqs_dirpath = os.path.join(working_dir, "unzipped_seqs")
    os.mkdir(unzipped_seqs_dirpath)
    relabed_seqs_dirpath = os.path.join(working_dir, "relabeled_seqs")
    os.mkdir(relabed_seqs_dirpath)
    seqs_stats_dfs_dirpath = os.path.join(working_dir, "seqs_stats_dfs")
    os.mkdir(seqs_stats_dfs_dirpath)
    pooled_seqs_fp = os.path.join(working_dir, "merged.fastq")
    
    if debug:
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
    
                # Could not find -fastx_relabel equivalent in vsearch
                if use_vsearch:
                    dna_seqs_gen = skbio.io.registry.read(
                        gzip_reader, format="fastq", verify=True, variant=variant)
                # skbio's fastq writer is way too slow, try to use usearch to relab the seqs if possible
                    i = 0
                    for seq in dna_seqs_gen:
                        i = i + 1
                        seq.metadata["id"] = sample_id + "." + str(i)
                        if not keep_annotations:
                            seq.metadata["description"] = ""
                        else:
                            seq.metadata["description"] = seq.metadata['description'].replace(
                                "\t", " ")
                        # pooled seqs use 1.8 encoding
                        seq.write(pooled_seqs_fh, format="fastq", variant="illumina1.8")
    
                # Usearch is much faster at relabeling seqs
                else:
    
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
                pipeout_denoise_stats_df.loc[index, 'prior_to_maxee_filt'] = i
                
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
                pipeout_denoise_stats_df.loc[index, 'prior_to_maxee_filt'] = int(stats_df.at[0, 'num_seqs'])
                

                
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

    return pipeout_denoise_stats_df


def _quality_control_cli(working_dir,
                         min_qscore=None, # Since Usearch Retains LQ reads during the final otutab stage, don't perform anthing other than maxee filtering
                         max_ee=1.0,
                         trim_left=0,
                         trunc_right=0,
                         min_len=50,  # Length filter is nessesary when dealing with valid data
                         max_ns=None,
                         threads="auto",
                         use_vsearch: bool = False,
                         verbose=True):

    pooled_seqs_fp = os.path.join(working_dir, "merged.fastq")
    filtered_reads_fp = os.path.join(working_dir, "filtered.fasta")
    log_fp = os.path.join(working_dir, "fastq_filter.log")

    # Building qc command

    if use_vsearch:
        cmd = ["vsearch"]
    else:
        cmd = ["usearch"]

    cmd += ["-fastq_filter", pooled_seqs_fp,
            "-fastaout", filtered_reads_fp,
            "-log", log_fp
            ]
    # add quality control parameters
    if min_qscore is not None:
        cmd += ["-fastq_truncqual", str(min_qscore)]
    if max_ee is not None:
        cmd += ["-fastq_maxee", str(max_ee)]
    if trim_left != 0:
        cmd += ["-fastq_stripleft", str(trim_left)]
    if trunc_right != 0:
        cmd += ["-fastq_trunclen", str(trunc_right)]
    if min_len is not None:
        cmd += ["-fastq_minlen", str(min_len)]
    if max_ns is not None:
        cmd += ["-fastq_maxns", str(max_ns)]
    # add threads settings
    if threads != "auto":
        cmd += ["-threads", str(threads)]

    # ring notification bell
    if verbose:
        print("Now performing maxEE QC on input reads...")

    # run command
    fastq_filter_log = py_to_cli_interface(cmd, verbose)

    # get stats from log file
    # fastq_filter_stats = [info for info in fastq_filter_log if "Filtered reads" in info][0]
    # fastq_filter_stats_lst = re.split("\s+", fastq_filter_stats)
    # reads_passed_filter = int(fastq_filter_stats_lst[fastq_filter_stats_lst.index("Filtered") - 1])

    # do stats on a per sample bases
    dna_seqs_gen = skbio.io.registry.read(
        filtered_reads_fp, format="fasta", verify=True)

    seq_id_lst = [seq_id.metadata['id'] for seq_id in dna_seqs_gen]

    stats_df = pd.DataFrame({'reads_passed_filter': seq_id_lst})

    stats_df['sample-id'] = stats_df["reads_passed_filter"].str.split(".", expand=True)[
        0]

    stats_df = stats_df.groupby('sample-id').count()

    return stats_df


def _dereplicate_cli(working_dir,
                     min_unique_size=None,
                     strand="plus",
                     threads="auto",
                     use_vsearch: bool = False,
                     verbose=True):

    filtered_reads_fp = os.path.join(working_dir, "filtered.fasta")
    unique_reads_fp = os.path.join(working_dir, "dereped.fasta")
    log_fp = os.path.join(working_dir, "fastx_uniques.log")

    # Building dereplication command

    if use_vsearch:
        cmd = ["vsearch"]
    else:
        cmd = ["usearch"]

    cmd += [
        "-fastx_uniques", filtered_reads_fp,
        "-fastaout", unique_reads_fp,
        "-sizeout",
        "-log", log_fp
    ]

    # Dereplication options
    if min_unique_size is not None:
        cmd += ["-minuniquesize", str(min_unique_size)]
    if strand == "both":
        cmd += ["-strand", "both"]
    if threads != "auto":
        cmd += ["-threads", str(threads)]

    # run command and get stats
    derep_log = py_to_cli_interface(cmd, verbose)

    # get stats from log file
    if use_vsearch:
        filtered_reads = "TBD"
        unique_reads = "TBD"
        singletons = "TBD"

    else:
        derep_stats = [
            info for info in derep_log if "uniques" in info and "singletons" in info][0]
        derep_stats_lst = re.split("\s+", derep_stats)
        derep_stats_lst = [text.replace(",", "") for text in derep_stats_lst]
        filtered_reads = int(derep_stats_lst[derep_stats_lst.index("seqs") - 1])
        unique_reads = int(derep_stats_lst[derep_stats_lst.index("uniques") - 1])
        singletons = int(derep_stats_lst[derep_stats_lst.index("singletons") - 1])

    return filtered_reads, unique_reads, singletons


def _unoise_cli(working_dir,
                min_size=8,
                unoise_alpha=2.0,
                use_vsearch: bool = False,
                verbose=True):

    unique_reads_fp = os.path.join(working_dir, "dereped.fasta")
    amplicon_fp = os.path.join(working_dir, "amps.fasta")
    vsearch_amplicon_fp = os.path.join(working_dir, "vsearch_amps.fasta")
    zotus_fp = os.path.join(working_dir, "zotus.fasta")
    log_fp = os.path.join(working_dir, "unoise.log")

    # Building unoise command
    # since the workflow is different for vsearch and usearch here, write two code blocks to handle them
    if not use_vsearch:
        cmd = ["usearch",
               "-unoise3", unique_reads_fp,
               "-ampout", amplicon_fp,
               "-log", log_fp
               ]

        if min_size != 8:
            cmd += ["-minsize", str(min_size)]

        if unoise_alpha != 2.0:
            cmd += ["-alpha", str(unoise_alpha)]

        unoise_log = py_to_cli_interface(cmd, verbose)

        unoise_stats = [info for info in unoise_log if "amplicons" in info][0]
        uchime_stats = [info for info in unoise_log if "good" in info][0]
        unoise_stats_lst = re.split("\s+", unoise_stats)
        uchime_stats_lst = re.split("\s+", uchime_stats)
        unoise_stats_lst = [text.replace(",", "") for text in unoise_stats_lst]
        uchime_stats_lst = [text.replace(",", "") for text in uchime_stats_lst]
        amplicons = int(unoise_stats_lst[unoise_stats_lst.index("amplicons") - 1])
        zotus = int(uchime_stats_lst[uchime_stats_lst.index("good") - 1])

    else:
        unoise_cmd = ["vsearch",
                      "--cluster_unoise", unique_reads_fp,
                      "--centroids", vsearch_amplicon_fp
                      ]

        if min_size != 8:
            unoise_cmd += ["-minsize", str(min_size)]

        if unoise_alpha != 2.0:
            unoise_cmd += ["-alpha", str(unoise_alpha)]

        py_to_cli_interface(unoise_cmd, verbose)

        uchime_cmd = ["vsearch",
                      "--uchime3_denovo", vsearch_amplicon_fp,
                      "--nonchimeras", zotus_fp,
                      "--relabel_md5"
                      ]

        py_to_cli_interface(uchime_cmd, verbose)

        amplicons = "TBD"
        zotus = "TBD"

    return amplicons, zotus

def _split_zotu_chimera(working_dir,
                        use_vsearch: bool = False,
                        verbose=True):
    amplicons_fp = os.path.join(working_dir, "amps.fasta")
    zotus_fp = os.path.join(working_dir, "zotus.fasta")
    chimeras_fp = os.path.join(working_dir, "chimeras.fasta")
    vsearch_amplicon_fp = os.path.join(working_dir, "vsearch_amps.fasta")

    if not use_vsearch:
        # used skbio and hashlib to hash the zotu ids
        dna_seqs_gen = skbio.io.registry.read(amplicons_fp, format="fasta", verify=True)
        with open(chimeras_fp, "wt") as chimeras_fh:
            with open(zotus_fp, "wt") as zotus_fh:
                # input seqs already sorted by decreasing abundance by usearch
                for seq in dna_seqs_gen:
                    if "amptype=chimera" in seq.metadata["id"]:
                        seq.metadata["id"] = hashlib.md5(str(seq).encode('utf-8')).hexdigest()
                        seq.write(chimeras_fh, format="fasta", max_width=80)
                    else:
                        seq.metadata["id"] = hashlib.md5(str(seq).encode('utf-8')).hexdigest()
                        seq.write(zotus_fh, format="fasta", max_width=80)

    else:
        # use skbio to separate chimeras from amplicons

        # Get zotu_ids
        zotu_id_lst = [zotu.metadata["id"] for zotu in skbio.io.registry.read(
            zotus_fp, format="fasta", verify=True)]
        # Get amplicon_ids
        amplicon_id_lst = [hashlib.md5(str(amplicon).upper().encode('utf-8')).hexdigest(
        ) for amplicon in skbio.io.registry.read(vsearch_amplicon_fp, format="fasta", verify=True)]
        amplicon_seqs_lst = [skbio.DNA(str(amplicon).upper(), metadata={'id': hashlib.md5(str(amplicon).encode(
            'utf-8')).hexdigest()}) for amplicon in skbio.io.registry.read(vsearch_amplicon_fp, format="fasta", verify=True)]

        # Get chimeras
        chimera_id_lst = [
            amplicon_id for amplicon_id in amplicon_id_lst if amplicon_id not in zotu_id_lst]
        chimera_seqs_lst = [
            seq for seq in amplicon_seqs_lst if seq.metadata["id"] in chimera_id_lst]

        # Write chimeras
        with open(chimeras_fp, "wt") as f:
            for seq in chimera_seqs_lst:
                seq.write(f, format="fasta", max_width=80)
                
    if verbose:
        print("Successfully split zotus and chimeras and converted to hashed ids")


def _cluster_zotus_cli(working_dir,
                       identity=0.99,
                       use_vsearch: bool = False,
                       verbose=True):
    zotus_fp = os.path.join(working_dir, "zotus.fasta")
    otus_fp = os.path.join(working_dir, "otus.fasta")
    
    # building uclust cmd
    # during the denoise step the output is sorted by size, seqs with higer abundance tend to be with lower noise
    if not use_vsearch:
        cmd = ["usearch",
               "-cluster_smallmem", zotus_fp,
               "-id", str(identity),
               "-centroids", otus_fp,
               "-sortedby", "other"
               ]
    else:
        cmd = ["vsearch",
               "--cluster_smallmem", zotus_fp,
               "--id", str(identity),
               "--centroids", otus_fp,
               "--usersort"
               ]
               
    silence = py_to_cli_interface(cmd, verbose)
    
    otus = len([seq for seq in skbio.io.registry.read(otus_fp, format="fasta", verify=True)])
    
    if verbose:
        print("Successfully clustered zotus into otus")
    
    return otus

def _build_zotu_tab_cli(working_dir,
                        threads="auto",
                        chimera_map="vsearch",
                        use_vsearch: bool = False,
                        verbose=True):

    raw_reads_fp = os.path.join(working_dir, "merged.fastq")
    filtered_reads_fp = os.path.join(working_dir, "filtered.fasta")
    zotus_fp = os.path.join(working_dir, "zotus.fasta")
    chimeras_fp = os.path.join(working_dir, "chimeras.fasta")
    tsv_otutab_fp = os.path.join(working_dir, "zotu_tab.tsv")
    tsv_chimeratab_fp = os.path.join(working_dir, "chimera_tab.tsv")
    matched_zotus_fp = os.path.join(working_dir, "matched_zotus.fasta")
    unmapped_reads_fp = os.path.join(working_dir, "unmapped.fasta")
    log_fp = os.path.join(working_dir, "otutab.log")
    node_thread_count = os.cpu_count()

    # Building otu table command
    if use_vsearch:
        cmd = ["vsearch",
               "--usearch_global", filtered_reads_fp,
               "--db", zotus_fp,
               # https://github.com/torognes/vsearch/issues/552
               "--strand", "plus"
               ]

    else:
        cmd = ["usearch",
               "-otutab", raw_reads_fp,
               "-zotus", zotus_fp
               ]

    cmd += [
        "-otutabout", tsv_otutab_fp,
        "-dbmatched", matched_zotus_fp,
        "-notmatched", unmapped_reads_fp,
        "-id", "1.0",
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

    # step2 search chimeras against unmatched fasta
    # build cmd
    if chimera_map == "usearch":
        if verbose:
            print("This step is performed just to keep a mental check...")
            print("Consider using vsearch to speed up this step...")
        chimera_cmd = ["usearch",
                       "-search_global", unmapped_reads_fp,
                       "-db", chimeras_fp,
                       "-id", "1.0",
                       "-strand", "both",
                       "-otutabout", tsv_chimeratab_fp
                       ]

        if threads != "auto":
            if threads > node_thread_count:
                if verbose:
                    print("Number of threads specified higher than max available on node...")
                    print("Setting threads to max available on node...")
                chimera_cmd += ["-threads", str(node_thread_count)]
            else:
                chimera_cmd += ["-threads", str(threads)]
        else:
            chimera_cmd += ["-threads", str(node_thread_count - 3)]
    else:
        chimera_cmd = ["vsearch",
                       "--usearch_global", unmapped_reads_fp,
                       "--db", chimeras_fp,
                       "--id", "1.0",
                       "--strand", "both",
                       "--otutabout", tsv_chimeratab_fp,
                       ]
        if threads != "auto":
            if threads > node_thread_count:
                if verbose:
                    print("Number of threads specified higher than max available on node...")
                    print("Setting threads to max available on node...")
                chimera_cmd += ["--threads", str(node_thread_count)]
            else:
                chimera_cmd += ["--threads", str(threads)]
        else:
            chimera_cmd += ["--threads", str(node_thread_count - 3)]

    # run command
    # we can do stats in another function
    if os.path.exists(chimeras_fp):
        # 2nd layer of insurance
        if len([ chimera for chimera in skbio.io.registry.read(chimeras_fp, format="fasta", verify=True)]) != 0:
            chimera_log = py_to_cli_interface(chimera_cmd, verbose)
            
            
def _uparse_cli(working_dir,
                min_size=2,
                verbose=True):
    unique_reads_fp = os.path.join(working_dir, "dereped.fasta")
    otus_fp = os.path.join(working_dir, "otus.fasta")
    chimeras_fp = os.path.join(working_dir, "chimeras.fasta")
    uparse_tab_fp = os.path.join(working_dir, "uparse_out.tsv")
    log_fp = os.path.join(working_dir, "uparse.log")
    
    # building uparse command
    cmd = ["usearch",
           "-cluster_otus", unique_reads_fp,
           "-otus", otus_fp,
           "-uparseout", uparse_tab_fp,
           "-log", log_fp
           ]
           
    if min_size != 2:
        cmd += ["-minsize", str(min_size)]
        
    silence = py_to_cli_interface(cmd, verbose)
    
    # get otu count and reformat otu_ids to qiime2 format
    dna_seqs_gen = skbio.io.registry.read(otus_fp, format="fasta", verify=True)
    otu_seqs = [ seq for seq in dna_seqs_gen ]
    otu_seq_count = len(otu_seqs)
    os.remove(otus_fp)
    with open(otus_fp, "wt") as f:
        for seq in otu_seqs:
            skbio.DNA(str(seq), metadata = {'id': hashlib.md5(str(seq).upper().encode('utf-8')).hexdigest()}).write(f, format="fasta", max_width=80)

    # indentify and retrive chimeras seqs
    # get chimera reads ids
    
    uparse_tab = pd.read_csv(uparse_tab_fp, sep="\t", header=None)
    
    try:
        chimera_seqs_id_list = uparse_tab.loc[uparse_tab[1] == "noisy_chimera", 0].str.split(";", expand = True)[0].to_list()
        chimera_seq_count = len(chimera_seqs_id_list)
        dna_seqs_gen = skbio.io.registry.read(unique_reads_fp, format="fasta", verify=True)
        chimera_seqs_lst = [ seq for seq in dna_seqs_gen if seq.metadata["id"].split(";")[0] in chimera_seqs_id_list ]
        
        # Write chimeras
        with open(chimeras_fp, "wt") as f:
            for seq in chimera_seqs_lst:
                seq.write(f, format="fasta", max_width=80)
                
    except KeyError:
        chimera_seq_count = 0
        if verbose:
            print("No chimera seqs found, skipping checking chimera in input...")

    if verbose:
        print("Denovo OTU clustering completed, chimeras splitted")
        
    return otu_seq_count , chimera_seq_count

def _build_otu_tab_cli(working_dir,
                       identity = 0.97,
                       threads="auto",
                       chimera_map="vsearch",
                       use_vsearch: bool = False,
                       verbose=True):
################################################################################
    raw_reads_fp = os.path.join(working_dir, "merged.fastq")
    filtered_reads_fp = os.path.join(working_dir, "filtered.fasta")
    otus_fp = os.path.join(working_dir, "otus.fasta")
    chimeras_fp = os.path.join(working_dir, "chimeras.fasta")
    tsv_otutab_fp = os.path.join(working_dir, "otu_tab.tsv")
    tsv_chimeratab_fp = os.path.join(working_dir, "chimera_tab.tsv")
    matched_otus_fp = os.path.join(working_dir, "matched_otus.fasta")
    unmapped_reads_fp = os.path.join(working_dir, "unmapped.fasta")
    log_fp = os.path.join(working_dir, "otutab.log")
    node_thread_count = os.cpu_count()

    # Building otu table command
    if use_vsearch:
        cmd = ["vsearch",
               "--usearch_global", filtered_reads_fp,
               "--db", otus_fp,
               # https://github.com/torognes/vsearch/issues/552
               "--strand", "plus"
               ]

    else:
        cmd = ["usearch",
               "-otutab", raw_reads_fp,
               "-otus", otus_fp
               ]

    cmd += [
        "-id", str(identity),
        "-otutabout", tsv_otutab_fp,
        "-dbmatched", matched_otus_fp,
        "-notmatched", unmapped_reads_fp,
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

    # step2 search chimeras against unmatched fasta
    # build cmd
    if chimera_map == "usearch":
        if verbose:
            print("This step is performed just to keep a mental check...")
            print("Consider using vsearch to speed up this step...")
        chimera_cmd = ["usearch",
                       "-search_global", unmapped_reads_fp,
                       "-db", chimeras_fp,
                       "-id", "1.0",
                       "-strand", "both",
                       "-otutabout", tsv_chimeratab_fp
                       ]

        if threads != "auto":
            if threads > node_thread_count:
                if verbose:
                    print("Number of threads specified higher than max available on node...")
                    print("Setting threads to max available on node...")
                chimera_cmd += ["-threads", str(node_thread_count)]
            else:
                chimera_cmd += ["-threads", str(threads)]
        else:
            chimera_cmd += ["-threads", str(node_thread_count - 3)]
    else:
        chimera_cmd = ["vsearch",
                       "--usearch_global", unmapped_reads_fp,
                       "--db", chimeras_fp,
                       "--id", "1.0",
                       "--strand", "both",
                       "--otutabout", tsv_chimeratab_fp,
                       ]
        if threads != "auto":
            if threads > node_thread_count:
                if verbose:
                    print("Number of threads specified higher than max available on node...")
                    print("Setting threads to max available on node...")
                chimera_cmd += ["--threads", str(node_thread_count)]
            else:
                chimera_cmd += ["--threads", str(threads)]
        else:
            chimera_cmd += ["--threads", str(node_thread_count - 3)]

    # run command
    # we can do stats in another function
    if os.path.exists(chimeras_fp):
        # 2nd layer of insurance
        if len([ chimera for chimera in skbio.io.registry.read(chimeras_fp, format="fasta", verify=True)]) != 0:
            chimera_log = py_to_cli_interface(chimera_cmd, verbose)
################################################################################


def _prep_results_for_artifact_api(working_dir,
                                   verbose=True):
    # get filepaths
    matched_zotus_fp = os.path.join(working_dir, "matched_zotus.fasta")
    matched_otus_fp = os.path.join(working_dir, "matched_otus.fasta")
    zotutab_fp = os.path.join(working_dir, "zotu_tab.tsv")
    otutab_fp = os.path.join(working_dir, "otu_tab.tsv")
    if os.path.exists(matched_zotus_fp) and os.path.exists(zotutab_fp):
        matched_features_fp = os.path.join(working_dir, "matched_zotus.fasta")
        tab_fp = os.path.join(working_dir, "zotu_tab.tsv")
        dt_type = "zotu"
    elif os.path.exists(matched_otus_fp) and os.path.exists(otutab_fp):
        matched_features_fp = os.path.join(working_dir, "matched_otus.fasta")
        tab_fp = os.path.join(working_dir, "otu_tab.tsv")
        dt_type = "otu"
    else:
        raise FileNotFoundError("Could not find u/vsearch output, pipeline broken...")

    chimeratab_fp = os.path.join(working_dir, "chimera_tab.tsv")

    # process feature_tab
    with open(tab_fp) as fh:
        tab_df = biom.Table.from_tsv(fh, None, None, None).to_dataframe()

    # get reads count mapped to zotus
    reads_mapped_to_features_df = tab_df.sum().to_frame().astype('int')
    if dt_type == "zotu":
        reads_mapped_to_features_df.columns = ["reads_mapped_to_zotus"]
    elif dt_type == "otu":
        reads_mapped_to_features_df.columns = ["reads_mapped_to_otus"]
    reads_mapped_to_features_df.index.name = "sample-id"

    # sort featrue tab
    tab_df['sum'] = tab_df.sum(axis=1)
    tab_df.sort_values(by="sum", ascending=False, inplace=True)
    tab_df.drop(columns=["sum"], inplace=True)
    ordered_zotu_ids = tab_df.index.tolist()

    table = biom.Table(tab_df.values, tab_df.index, tab_df.columns)

    # process chimeratab
    if os.path.exists(chimeratab_fp):
        with open(chimeratab_fp) as fh:
            chimera_df = biom.Table.from_tsv(fh, None, None, None).to_dataframe()

        # get reads count mapped to chimeras
        reads_mapped_to_chimeras_df = chimera_df.sum().to_frame().astype('int')
        reads_mapped_to_chimeras_df.columns = ["reads_mapped_to_chimeras"]
        reads_mapped_to_chimeras_df.index.name = "sample-id"
        
    else:
        reads_mapped_to_chimeras_df = pd.DataFrame({"reads_mapped_to_chimeras": 0}, index=reads_mapped_to_features_df.index)

    if verbose:
        print("Now sorting features accroding to feature tab...")

    # process zotus
    # if the query seqs is so bad, it is possible for u/vsearch to map a feature with higher frequency than in the filtered seqs
    # so we have to sort the rep-seqs here again accroding to the feature tab
    rep_seqs_lst = [ seq for seq in skbio.io.registry.read(matched_features_fp, format="fasta", verify=True) ]
    rep_seqs_id_lst = [ seq.metadata["id"] for seq in rep_seqs_lst ]
    rep_sequences = pd.Series(rep_seqs_lst, index = rep_seqs_id_lst)
    rep_sequences = rep_sequences.reindex(tab_df.index)
    if rep_sequences.isna().sum() > 0:
        raise ValueError("DEBUG: Some features in feature table is not in rep-seqs...")

    if verbose:
        if dt_type == "zotu":
            print("Successfully sorted zotutab and zotus...")
        elif dt_type == "otu":
            print("Successfully sorted otutab and otus ...")

    return table, rep_sequences, reads_mapped_to_features_df, reads_mapped_to_chimeras_df

def denoise_no_primer_pooled(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                                   trim_left: int = 0,
                                   trunc_len: int = 0,
                                   min_len: int = 50,
                                   max_ee: float = 1.0,
                                   n_threads: str = "auto",
                                   min_size: int = 8,
                                   unoise_alpha: float = 2.0,
                                   use_vsearch: bool = False,
                                   ) -> (biom.Table, pd.Series, qiime2.Metadata):
                                       
    verbose = True

    demultiplexed_sequences_dirpath = str(demultiplexed_sequences)
    with tempfile.TemporaryDirectory() as usearch_wd:
        input_stats_df = _pool_samples(
            demultiplexed_sequences_dirpath, usearch_wd, use_vsearch=use_vsearch, verbose=verbose)
        # need to sep for each sample as well
        filter_stats_df = _quality_control_cli(usearch_wd, trim_left=trim_left, trunc_right=trunc_len,
                                               min_len=min_len, max_ee=max_ee, use_vsearch=use_vsearch, threads=n_threads, verbose=verbose)
    
        filtered_reads_count, unique_reads_count, singletons_count, = _dereplicate_cli(
            usearch_wd, use_vsearch=use_vsearch, threads=n_threads, verbose=verbose)
        amplicons_count, zotus_count, = _unoise_cli(
            usearch_wd, min_size=min_size, unoise_alpha=unoise_alpha, use_vsearch=use_vsearch, verbose=verbose)
    
        denoise_stats_str = "Total Reads: " + str(filtered_reads_count) + " ;Unique Reads :" + str(unique_reads_count) + \
            " ;Singletons: " + str(singletons_count) + " ;Amplicons: " + \
            str(amplicons_count) + " ;ZOTUs: " + str(zotus_count)
    
        _split_zotu_chimera(usearch_wd, use_vsearch=use_vsearch, verbose=verbose)
        _build_zotu_tab_cli(usearch_wd, use_vsearch=use_vsearch,
                            threads=n_threads, verbose=verbose)
        table, representative_sequences, reads_mapped_to_zotus_df, reads_mapped_to_chimeras_df = _prep_results_for_artifact_api(
            usearch_wd, verbose=verbose)
    
        # finally prep denoise stats df
        denoise_stats_df = input_stats_df.merge(
            filter_stats_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_passed_filter'] = (
            denoise_stats_df['reads_passed_filter'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_zotus_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_zotus'] = (
            denoise_stats_df['reads_mapped_to_zotus'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_chimeras_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_chimeras'] = (
            denoise_stats_df['reads_mapped_to_chimeras'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df["denoise_stats_pooled_mode"] = denoise_stats_str
        
        # if sample ids were swapped during the run, we need to swap the sample ids back
        if 'original_sample_id' in denoise_stats_df.columns:
            
            id_map_dict = denoise_stats_df['original_sample_id'].to_dict()
            denoise_stats_df.index = denoise_stats_df['original_sample_id']
            denoise_stats_df.index.name = 'sample-id'
            denoise_stats_df.drop(columns=['original_sample_id'], inplace=True)
            table.update_ids(id_map_dict, axis='sample', inplace=True)
            
        denoise_stats_df.fillna(0, inplace=True)
        
        denoising_stats = qiime2.Metadata(denoise_stats_df)
        
    return table, representative_sequences, denoising_stats

# do we need to expose additional uparse parameters here?
def cluster_no_primer_pooled(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                                   trim_left: int = 0,
                                   trunc_len: int = 0,
                                   min_len: int = 50,
                                   max_ee: float = 1.0,
                                   n_threads: str = "auto",
                                   min_size: int = 2,
                                   ) -> (biom.Table, pd.Series, qiime2.Metadata):
                                       
    verbose = True
    
    if verbose:
        print("Since Usearch version 9.0.2132, the ability to directly cluster OTUs to custom identity threshould had been removed. ")
        print("The reason for this is that indentity threshold other than 0.97 mess up the chimera detection step. ")
        print("Further expalnation can be found here: https://drive5.com/usearch/manual/uparse_otu_radius.html")
        print("BTW uparse is also usearch exclusive, no vsearch support here.")

    demultiplexed_sequences_dirpath = str(demultiplexed_sequences)
    with tempfile.TemporaryDirectory() as usearch_wd:
        input_stats_df = _pool_samples(
            demultiplexed_sequences_dirpath, usearch_wd, verbose=verbose)
        # need to sep for each sample as well
        filter_stats_df = _quality_control_cli(usearch_wd, trim_left=trim_left, trunc_right=trunc_len,
                                               min_len=min_len, max_ee=max_ee, threads=n_threads, verbose=verbose)
        # discarding singletons seemed not nessaary since the clust_otu command does it anyway...
        # need to check your self
        filtered_reads_count, unique_reads_count, singletons_count, = _dereplicate_cli(
            usearch_wd, threads=n_threads, verbose=verbose)
        otu_seq_count, chimera_seq_count, = _uparse_cli(
            usearch_wd, min_size=min_size, verbose=verbose)
    
        denoise_stats_str = "Total Reads: " + str(filtered_reads_count) + " ;Unique Reads :" + str(unique_reads_count) + \
            " ;Singletons: " + str(singletons_count) + " ;OTUs: " + str(otu_seq_count) + \
            " :Chimeras: " + str(chimera_seq_count)

        _build_otu_tab_cli(usearch_wd,
                            threads=n_threads, verbose=verbose)
        table, representative_sequences, reads_mapped_to_otus_df, reads_mapped_to_chimeras_df = _prep_results_for_artifact_api(
            usearch_wd, verbose=verbose)
    
        # finally prep denoise stats df
        denoise_stats_df = input_stats_df.merge(
            filter_stats_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_passed_filter'] = (
            denoise_stats_df['reads_passed_filter'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_otus_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_otus'] = (
            denoise_stats_df['reads_mapped_to_otus'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_chimeras_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_chimeras'] = (
            denoise_stats_df['reads_mapped_to_chimeras'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df["denoise_stats_pooled_mode"] = denoise_stats_str
        
        # if sample ids were swapped during the run, we need to swap the sample ids back
        if 'original_sample_id' in denoise_stats_df.columns:
            
            id_map_dict = denoise_stats_df['original_sample_id'].to_dict()
            denoise_stats_df.index = denoise_stats_df['original_sample_id']
            denoise_stats_df.index.name = 'sample-id'
            denoise_stats_df.drop(columns=['original_sample_id'], inplace=True)
            table.update_ids(id_map_dict, axis='sample', inplace=True)
            
        denoise_stats_df.fillna(0, inplace=True)
        
        denoising_stats = qiime2.Metadata(denoise_stats_df)
        
    return table, representative_sequences, denoising_stats

def denoise_then_cluster_no_primer_pooled(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                                   trim_left: int = 0,
                                   trunc_len: int = 0,
                                   min_len: int = 50,
                                   max_ee: float = 1.0,
                                   perc_identity: float = 0.99,
                                   n_threads: str = "auto",
                                   min_size: int = 8,
                                   unoise_alpha: float = 2.0,
                                   use_vsearch: bool = False,
                                   ) -> (biom.Table, pd.Series, qiime2.Metadata):
                                       
    verbose = True

    demultiplexed_sequences_dirpath = str(demultiplexed_sequences)
    with tempfile.TemporaryDirectory() as usearch_wd:
        input_stats_df = _pool_samples(
            demultiplexed_sequences_dirpath, usearch_wd, use_vsearch=use_vsearch, verbose=verbose)
        # need to sep for each sample as well
        filter_stats_df = _quality_control_cli(usearch_wd, trim_left=trim_left, trunc_right=trunc_len,
                                               min_len=min_len, max_ee=max_ee, use_vsearch=use_vsearch, threads=n_threads, verbose=verbose)
    
        filtered_reads_count, unique_reads_count, singletons_count, = _dereplicate_cli(
            usearch_wd, use_vsearch=use_vsearch, threads=n_threads, verbose=verbose)
        amplicons_count, zotus_count, = _unoise_cli(
            usearch_wd, min_size=min_size, unoise_alpha=unoise_alpha, use_vsearch=use_vsearch, verbose=verbose)
    
        _split_zotu_chimera(usearch_wd, use_vsearch=use_vsearch, verbose=verbose)
        
        otus_count = _cluster_zotus_cli(usearch_wd, identity=perc_identity, use_vsearch=use_vsearch, verbose=verbose)
        
        denoise_stats_str = "Total Reads: " + str(filtered_reads_count) + " ;Unique Reads :" + str(unique_reads_count) + \
            " ;Singletons: " + str(singletons_count) + " ;Amplicons: " + \
            str(amplicons_count) + " ;ZOTUs: " + str(zotus_count) + " ;OTUs: " + str(otus_count)
        _build_otu_tab_cli(usearch_wd, identity=perc_identity, use_vsearch=use_vsearch,
                            threads=n_threads, verbose=verbose)
        table, representative_sequences, reads_mapped_to_otus_df, reads_mapped_to_chimeras_df = _prep_results_for_artifact_api(
            usearch_wd, verbose=verbose)
    
        # finally prep denoise stats df
        denoise_stats_df = input_stats_df.merge(
            filter_stats_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_passed_filter'] = (
            denoise_stats_df['reads_passed_filter'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_otus_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_otus'] = (
            denoise_stats_df['reads_mapped_to_otus'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df = denoise_stats_df.merge(
            reads_mapped_to_chimeras_df, how="left", left_index=True, right_index=True)
        denoise_stats_df['percent_of_input_mapped_to_chimeras'] = (
            denoise_stats_df['reads_mapped_to_chimeras'] / denoise_stats_df['prior_to_maxee_filt']) * 100
        denoise_stats_df["denoise_stats_pooled_mode"] = denoise_stats_str
        
        # if sample ids were swapped during the run, we need to swap the sample ids back
        if 'original_sample_id' in denoise_stats_df.columns:
            
            id_map_dict = denoise_stats_df['original_sample_id'].to_dict()
            denoise_stats_df.index = denoise_stats_df['original_sample_id']
            denoise_stats_df.index.name = 'sample-id'
            denoise_stats_df.drop(columns=['original_sample_id'], inplace=True)
            table.update_ids(id_map_dict, axis='sample', inplace=True)
            
        denoise_stats_df.fillna(0, inplace=True)
        
        denoising_stats = qiime2.Metadata(denoise_stats_df)
        
    return table, representative_sequences, denoising_stats

