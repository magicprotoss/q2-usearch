# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# register usearch_stats_format

from qiime2.plugin import SemanticType, model
from q2_types.sample_data import SampleData

USEARCHStats = SemanticType('USEARCHStats', variant_of=SampleData.field['type'])

class USEARCHStatsFormat(model.TextFileFormat):
    def validate(*args):
        pass

USEARCHStatsDirFmt = model.SingleFileDirectoryFormat(
    'USEARCHStatsDirFmt', 'stats.tsv', USEARCHStatsFormat)
