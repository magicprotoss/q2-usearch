# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import subprocess
import tempfile
import skbio
import os
from q2_types.feature_data import DNAFASTAFormat
from q2_usearch._format import USEARCHFastQFmt
from glob import glob
import biom
import pandas as pd
import hashlib


def run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def _build_feature_table(merged_sequences,
                         representative_sequences,
                         rep_seqs_type,
                         n_threads):

    merged_sequences_fp = str(merged_sequences)
    representative_sequences_fp = str(representative_sequences)

    # generate id list file using ids from rep-seqs
    rep_seqs_gen = skbio.read(
        representative_sequences_fp, format='fasta', constructor=skbio.DNA)

    rep_seqs_lst = []
    for seq_obj in rep_seqs_gen:
        rep_seqs_lst.append(seq_obj)

    id_lst = []
    for idx in range(len(rep_seqs_lst)):
        id = rep_seqs_lst[idx].metadata['id']
        id_lst.append(id)

    # generate hashed ids (need tov rethink method registration)
    # id_map['hashed_id'] = ''

    # for idx in range(len(id_map)):
    #     id_map.loc[id_map.index[idx], 'hashed_id'] = hashlib.md5(
    #         str(id_map.loc[id_map.index[idx], 'seq']).encode('utf-8')).hexdigest()

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        tsv_feature_tab_fp = temp_dir + "/feature_tab.tsv"
        # Building otu-tab cmd
        if n_threads == "auto":
            n_threads = os.cpu_count() - 1
        if rep_seqs_type == "ZOTU":
            db_type = "-zotus"
            identity = 1.0
        elif rep_seqs_type == "OTU":
            db_type = "-otus"
            identity = 0.97
        otutab_cmd = ['usearch',
                      '-otutab', merged_sequences_fp,
                      db_type, representative_sequences_fp,
                      '-strand', 'plus',
                      '-id', str(identity),
                      '-otutabout', tsv_feature_tab_fp,
                      '-threads', str(n_threads)]
        run_command(otutab_cmd)

        # Convert feature tab from classic QIIME format to biom 2.1 format
        with open(tsv_feature_tab_fp) as tsv_tab_file:
            feature_table = biom.Table.from_tsv(tsv_tab_file, None, None, None)

            # Fix feature table ids
            feature_table = feature_table.sort_order(
                id_lst, axis="observation")

    return feature_table


def build_feature_table(merged_sequences: USEARCHFastQFmt,
                        representative_sequences: DNAFASTAFormat,
                        rep_seqs_type: str,
                        n_threads: str = "auto") -> biom.Table:
    return _build_feature_table(merged_sequences=merged_sequences,
                                representative_sequences=representative_sequences,
                                rep_seqs_type=rep_seqs_type,
                                n_threads=n_threads)
