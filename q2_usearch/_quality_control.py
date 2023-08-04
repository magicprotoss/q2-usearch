# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

# F... QIIME2 Currently does not have a underline for format for a single fastq file...

import subprocess
from ._format import USEARCHFastQFmt, USEARCHFastaFmt

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

def _quality_control(merged_sequences,
                     max_error_rate,
                     truncate_single,
                     min_sequence_len):
    clean_reads = USEARCHFastaFmt()
    clean_reads_fp = str(clean_reads)

    merged_sequences_fp = str(merged_sequences)

    # Building qc command
    cmd = ["usearch",
           "-fastq_filter",merged_sequences_fp,
           "-fastaout", clean_reads_fp,
           "-fastq_maxee", str(max_error_rate),
           "-relabel", "Clean"]
    if truncate_single != 0:
        cmd += ["-fastq_trunclen", str(truncate_single)]
    if truncate_single == 0 and min_sequence_len != 0:
        cmd += ["-fastq_minlen", str(min_sequence_len)]

    run_command(cmd)

    return clean_reads

def quality_control(merged_sequences: USEARCHFastQFmt,
                 max_error_rate : float = 1.0,
                 truncate_single: int = 0,
                 min_sequence_len: int = 0) -> USEARCHFastaFmt:
    return _quality_control(merged_sequences = merged_sequences,
                            max_error_rate = max_error_rate,
                            truncate_single = truncate_single,
                            min_sequence_len = min_sequence_len)