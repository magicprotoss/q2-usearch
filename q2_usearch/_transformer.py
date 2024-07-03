# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
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

