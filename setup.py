# ----------------------------------------------------------------------------
# Copyright (c) 2024, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name="q2-usearch",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="magicprotoss; biodps",
    author_email="magicprotoss@hotmail.com",
    description="This is a plug-in for USEARCH integration into QIIME2.",
    url="https://github.com/magicprotoss/q2-usearch",
    entry_points={
        "qiime2.plugins": ["q2-usearch=q2_usearch.plugin_setup:plugin"]
    },
    package_data={
        "q2_usearch": ["citations.bib"],
    },
    zip_safe=False,
)
