# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
#
# I hope Dr.Edgar won't get mad though...
# ----------------------------------------------------------------------------

import qiime2

from ._format import USEARCHFastaFmt, USEARCHFastQFmt, USEARCHFastaDirFmt, USEARCHFastQDirFmt
from .plugin_setup import plugin
from q2_types.feature_data import PairedDNASequencesDirectoryFormat


@plugin.register_transformer

def _1(data: PairedDNASequencesDirectoryFormat) -> USEARCHFastaFmt:
    ff = USEARCHFastaFmt()
    with ff.open() as fh:
        data.left_dna_sequences.write_data(fh)
        data.right_dna_sequences.write_data(fh)
    return ff

# # None of these work, qiime will refuse to work without transformer
# # however, reg a dummy working tranformer _1 will solve the problem for now...
def _2(data: USEARCHFastaDirFmt) -> USEARCHFastaFmt:
    ff = USEARCHFastaFmt()
    with ff.open() as fh:
        data.write(fh)
    return ff

def _3(data: USEARCHFastQDirFmt) -> USEARCHFastQFmt:
    ff = USEARCHFastQFmt()
    with ff.open() as fh:
        data.write(fh)
    return ff


