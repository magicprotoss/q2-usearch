# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import subprocess
import tempfile
from ._format import USEARCHFastaFmt

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

def _dereplication(clean_reads, n_threads):
    clean_reads_fp = str(clean_reads)
    dereplicated_sequences = USEARCHFastaFmt()
    dereplicated_sequences_fp = str(dereplicated_sequences)

    # F input fmt = output fmt means bug out...
    with tempfile.TemporaryDirectory() as working_dir:

        usearch_output = working_dir + "/dereplicated_seqs.fasta"

        # Building dereplication command
        dereplication_cmd = ['usearch',
                            '-fastx_uniques', clean_reads_fp,
                            '-fastaout', usearch_output,
                            '-relabel', 'Uni',
                            '-sizeout']
        if n_threads != "auto":
            dereplication_cmd += ['-threads', str(n_threads)]
        run_command(dereplication_cmd)

        with open(dereplicated_sequences_fp, 'wb') as des:
            with open(usearch_output, 'rb') as uo:
                des.write(uo.read())

    return dereplicated_sequences

def dereplication(clean_reads: USEARCHFastaFmt,
                 n_threads: str = "auto") -> USEARCHFastaFmt:
    return _dereplication(clean_reads = clean_reads,
                            n_threads = n_threads)