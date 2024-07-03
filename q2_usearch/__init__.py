# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._illumina_pipeline import denoise_no_primer_pooled, cluster_no_primer_pooled, denoise_then_cluster_no_primer_pooled

# modified from q2-vsearch
from ._merge_pairs import merge_pairs

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

__all__ = ['denoise_no_primer_pooled', 'cluster_no_primer_pooled', 'denoise_then_cluster_no_primer_pooled', 'merge_pairs']
