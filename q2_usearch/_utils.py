# ----------------------------------------------------------------------------
# Copyright (c) 2026, magicprotoss;biodps.
#
# Distributed under the terms of the Modified BSD License.
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import hashlib
import re

import biom
import pandas as pd
import skbio


def fix_labels(feature_table: biom.Table, rep_seqs: pd.Series,
               to_md5: bool = False, to_legacy: bool = False,
               tidy_only: bool = False, label: str = None
               ) -> tuple[biom.Table, pd.Series]:
    """Tidy a table and its sequences, optionally assigning matching new IDs."""
    if sum((to_md5, to_legacy, tidy_only)) != 1:
        raise ValueError(
            "Please select one method only(--p-to-md5 / --p-to-legacy / "
            "--p-tidy-only).")
    if to_legacy:
        if label is None or re.fullmatch(r'[A-Za-z0-9_]+', label) is None:
            raise ValueError(
                '"--p-to-legacy" requires setting a nonempty label via '
                '"--p-label" containing only '
                "ASCII letters (A-Z, a-z), digits (0-9), and underscores (_).")
    elif label is not None:
        raise ValueError(
            'Set "--p-label" only when "--p-to-legacy" was selected.')

    if rep_seqs.index.has_duplicates:
        raise ValueError(
            "Representative sequence IDs must be unique. "
            "check if you imported the wrong file into QIIME2.")
    feature_ids = feature_table.ids(axis='observation')
    sequence_ids = set(rep_seqs.index)
    missing = [feature_id for feature_id in feature_ids
               if feature_id not in sequence_ids]
    if missing:
        raise ValueError(
            "All feature table IDs must be present in rep_seqs, including "
            "zero-abundance features. Missing IDs: " + ", ".join(missing))

    empty_samples = [sample_id for sample_id, total in zip(
        feature_table.ids(axis='sample'), feature_table.sum(axis='sample'))
        if total == 0]
    if empty_samples:
        raise ValueError(
            "Some input samples have a total abundance of 0 in the table, "
            "please check the data.\n"
            "Empty sample IDs: " + ", ".join(empty_samples))

    totals = feature_table.sum(axis='observation')
    # Python's stable sort preserves original table order when totals tie.
    retained_ids = [feature_id for feature_id, total in sorted(
        zip(feature_ids, totals), key=lambda item: item[1], reverse=True)
        if total > 0]
    if not retained_ids:
        raise ValueError(
            "No features with nonzero abundance remain after tidying. "
            "Please check the input table.")

    sequences = [str(rep_seqs.loc[feature_id]).upper()
                 for feature_id in retained_ids]
    if to_md5:
        new_ids = [hashlib.md5(sequence.encode('utf-8')).hexdigest()
                   for sequence in sequences]
        if len(set(new_ids)) != len(new_ids):
            raise ValueError(
                "Duplicate MD5 feature IDs: retained sequences are identical "
                "after uppercasing, or their MD5 hashes collide. "
                "Resolve the duplicates before using to_md5; features are "
                "not merged automatically.")
    elif to_legacy:
        new_ids = [f"{label}{i}" for i in range(1, len(retained_ids) + 1)]
    else:
        new_ids = retained_ids

    fixed_table = feature_table.filter(
        retained_ids, axis='observation', inplace=False)
    fixed_table = fixed_table.sort_order(retained_ids, axis='observation')
    fixed_table.update_ids(dict(zip(retained_ids, new_ids)), axis='observation')
    fixed_seqs = pd.Series(
        [skbio.DNA(sequence, metadata={'id': feature_id})
         for feature_id, sequence in zip(new_ids, sequences)],
        index=new_ids, name=rep_seqs.name, dtype=object)
    return fixed_table, fixed_seqs
