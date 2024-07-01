# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
# ----------------------------------------------------------------------------

import qiime2
from .plugin_setup import plugin
from q2_usearch._format import USEARCHStatsFormat

@plugin.register_transformer
def _1(ff: USEARCHStatsFormat) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))


@plugin.register_transformer
def _2(obj: qiime2.Metadata) -> USEARCHStatsFormat:
    ff = USEARCHStatsFormat()
    obj.save(str(ff))
    return ff

