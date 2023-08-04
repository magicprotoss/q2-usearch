# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import subprocess
import tempfile
import os
from q2_types.feature_data import DNAFASTAFormat
from q2_usearch._format import USEARCHFastQFmt
from glob import glob
import biom

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

    return feature_table


def build_feature_table(merged_sequences: USEARCHFastQFmt,
                        representative_sequences: DNAFASTAFormat,
                        rep_seqs_type: str,
                 n_threads : str = "auto") -> biom.Table:
    return _build_feature_table(merged_sequences = merged_sequences,
                                representative_sequences = representative_sequences,
                                rep_seqs_type = rep_seqs_type,
                                n_threads = n_threads)