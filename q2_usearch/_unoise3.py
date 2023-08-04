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

def _unoise3(dereplicated_reads, minsize, alpha):
    dereplicated_reads_fp = str(dereplicated_reads)
    # creating temporary files
    zotus = DNAFASTAFormat()
    zotus_fp = str(zotus)

    # building unoise command
    cmd = ["usearch", "-unoise3", dereplicated_reads_fp, "-zotus", zotus_fp]
    if minsize != 8:
        cmd += ["-minsize", str(minsize)]
    if alpha != 2.0:
        cmd += ["-unoise_alpha", str(alpha)]

    # Usearch output fasta files are folded, try making unfolding them
    run_command(cmd)
    # if otu_tab command accepts non Zotu labels, add a relabeling step here
    return zotus

def unoise3(dereplicated_reads: USEARCHFastaFmt,
            minsize: int = 8,
            alpha: float = 2.0) -> DNAFASTAFormat:
    return _unoise3(dereplicated_reads, minsize, alpha)

