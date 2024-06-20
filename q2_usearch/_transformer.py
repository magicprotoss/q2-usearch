# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
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

# Register Usearch stats fmt
from q2_usearch._format import USEARCHStatsFormat

@plugin.register_transformer
def _4(ff: USEARCHStatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _5(obj: qiime2.Metadata) -> USEARCHStatsFormat:
    ff = USEARCHStatsFormat()
    obj.save(str(ff))
    return ff

