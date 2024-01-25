# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

# F... QIIME2 Currently does not have a underline for format for a single fastq file...

import subprocess
import tempfile
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt
from ._format import USEARCHFastQFmt
from glob import glob
import shutil
import gzip
import pandas as pd
import csv


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


def _pool_samples(demultiplexed_sequences, keep_annotations):
    merged_sequences = USEARCHFastQFmt()
    merged_sequences_fp = str(merged_sequences)
    demultiplexed_sequences_dir_path = str(demultiplexed_sequences)
    # Create temp dir for usearch to work around
    with tempfile.TemporaryDirectory() as working_dir:

        # Need to further imporve by writing usearch idnetifier exception
        manifest_fp = demultiplexed_sequences_dir_path + "/MANIFEST"
        # Fix from stack overflow, it seemed early versions of qiime2 use csv for MANIFEST file
        manifest_df = pd.read_csv(manifest_fp, comment='#', sep=None)
        mapping = manifest_df.iloc[:, 0:2]
        id_map = pd.DataFrame()
        id_map[['id', 'fn']] = mapping

        # unzip and relabel all samples
        print("unzipping all samples...")

        for index, row in id_map.iterrows():
            fp = demultiplexed_sequences_dir_path + "/" + str(row['fn'])
            prefix = str(row['id']) + "."
            with gzip.open(fp, "rb") as unzipped_seqs:
                with tempfile.NamedTemporaryFile() as rf:
                    rf.write(unzipped_seqs.read())
                    rf.flush()  # flush the data to disk
                    with tempfile.NamedTemporaryFile() as wf:
                        relabel_cmd = ['usearch', '-fastx_relabel',
                                       rf.name, '-prefix', prefix, '-fastqout', wf.name]
                        if keep_annotations == True:
                            relabel_cmd += ["-keep_annots"]
                        run_command(relabel_cmd)
                        wf.flush()  # flush the data to disk
                        with open(merged_sequences_fp, "ab") as merged_seqs:
                            merged_seqs.write(wf.read())
                            merged_seqs.flush()

    return merged_sequences


def pool_samples(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                 keep_annotations: bool = False) -> USEARCHFastQFmt:
    return _pool_samples(demultiplexed_sequences=demultiplexed_sequences, keep_annotations=keep_annotations)
