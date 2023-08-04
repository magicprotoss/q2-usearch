# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import subprocess
from q2_types.feature_data import DNAFASTAFormat
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

def _uprase_otu(sequences_fp, minsize = 2, alpha = 2.0):
    # creating temporary files
    otus = DNAFASTAFormat()
    otus_fp = str(otus)

    # building command
    cmd = ["usearch", "-cluster_otus", sequences_fp, "-otus", otus_fp, "-relabel", "Otu"]
    if minsize != 2:
        cmd += ["-minsize", str(minsize)]

    run_command(cmd)
    return otus

def uprase_otu(dereplicated_reads: USEARCHFastaFmt,
            minsize: int = 2) -> DNAFASTAFormat:
    dereplicated_reads_fp = str(dereplicated_reads)
    return _uprase_otu(dereplicated_reads_fp, minsize)