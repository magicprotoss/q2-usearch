# ----------------------------------------------------------------------------
# Copyright (c) 2022, <developer name>.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._unoise3 import unoise3
from ._uprase_otu import uprase_otu
from ._pool_samples import pool_samples
from ._quality_control import quality_control
from ._dereplication import dereplication
from ._pipeline import lazy_process_clean_data
from ._sintax import sintax
# from ._dereplication import dereplication
from ._build_feature_table import build_feature_table

# from ._format import (USEARCHFastaFmt, USEARCHFastQFmt, USEARCHFastaDirFmt, USEARCHFastQDirFmt)

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

__all__ = ['unoise3', 'uprase_otu', 'pool_samples', 'build_feature_table',
           'quality_control', 'dereplication', 'lazy_process_clean_data', 'sintax']
