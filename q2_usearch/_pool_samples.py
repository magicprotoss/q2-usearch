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
import gzip
import pandas as pd


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
    glob_gz_path = demultiplexed_sequences_dir_path + "/*.gz"
    # Create temp dir for usearch to work around
    with tempfile.TemporaryDirectory() as working_dir:

        glob_fq_path = working_dir + "/*.fastq"

        # Need to further imporve by writing usearch idnetifier exception
        manifest_fp = demultiplexed_sequences_dir_path + "/MANIFEST"
        manifest_df = pd.read_csv(manifest_fp)
        mapping = manifest_df.iloc[:, 0:2]
        id_map = pd.DataFrame()
        id_map[['id', 'fn']] = mapping
        id_map['old_id'] = id_map['fn'].str.split(".", expand=True)[0]
        id_map.drop(['fn'], axis=1, inplace=True)

        # unzip all samples
        print("unzipping all samples...")

        for zipped_seqs in glob(glob_gz_path):
            zipped_fn = zipped_seqs.split("/")[-1].split(".")[0]
            fn = id_map.loc[id_map['old_id'] ==
                            zipped_fn].reset_index().at[0, 'id']
            file_name = fn + ".fastq"
            with gzip.open(zipped_seqs, "rb") as unzipped_seqs:
                with open(working_dir + "/" + file_name, "wb") as wf:
                    wf.write(unzipped_seqs.read())

        # Fix sample labels
        print("Fixing sample labels...")

        for seq_fp in glob(glob_fq_path):
            seq_file_name = seq_fp.split("/")[-1]
            relabeled_seq_file_name = seq_file_name.split(".")[
                0] + ".relabeled"
            seq_labels = seq_file_name + "."

            relabeled_seq_fp = working_dir + "/" + relabeled_seq_file_name

            relabel_cmd = ["usearch", "-fastx_relabel", seq_fp,
                           "-prefix", seq_labels, "-fastqout", relabeled_seq_fp]

            if keep_annotations == True:
                relabel_cmd += ["-keep_annots"]

            del_original_cmd = ["rm", "-f", seq_fp]

            run_command(relabel_cmd)
            run_command(del_original_cmd)

        glob_relabeled_path = working_dir + "/*.relabeled"

        # merge all samples
        print("Now merging all samples...")

        with open(merged_sequences_fp, "wb") as merged_seqs:
            for seqs in glob(glob_relabeled_path):
                with open(seqs, "rb") as seqs_in:
                    merged_seqs.write(seqs_in.read())

    return merged_sequences


def pool_samples(demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
                 keep_annotations: bool = False) -> USEARCHFastQFmt:
    return _pool_samples(demultiplexed_sequences=demultiplexed_sequences, keep_annotations=keep_annotations)
