# ----------------------------------------------------------------------------
# This is a simple pulgin for usearch intergration in qiime2
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
