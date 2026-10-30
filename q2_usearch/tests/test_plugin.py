"""Exercise fix-labels through the registered QIIME 2 action interface."""

import hashlib
import itertools

import biom
import numpy as np
import pandas as pd
import pytest
import skbio
from qiime2 import Artifact
from qiime2.sdk import PluginManager


@pytest.fixture(scope='module')
def fix_labels_action():
    return PluginManager().plugins['usearch'].methods['fix_labels']


@pytest.fixture(scope='module')
def label_artifacts():
    table = biom.Table(
        np.array([[1, 0], [0, 0], [2, 3]]),
        ['low', 'zero', 'high'], ['sample-a', 'sample-b'])
    sequence_ids = ['low', 'extra', 'zero', 'high']
    sequences = pd.Series(
        [skbio.DNA(sequence, metadata={'id': identifier})
         for identifier, sequence in zip(
             sequence_ids, ['AAAA', 'GGGG', 'CCCC', 'ACGT'])],
        index=sequence_ids, dtype=object)
    return {
        'feature_table': Artifact.import_data('FeatureTable[Frequency]', table),
        'rep_seqs': Artifact.import_data('FeatureData[Sequence]', sequences),
    }


@pytest.mark.parametrize('to_md5,to_legacy,tidy_only',
                         itertools.product([False, True], repeat=3))
def test_registered_mode_constraints(
        fix_labels_action, label_artifacts, to_md5, to_legacy, tidy_only):
    parameters = {
        'to_md5': to_md5, 'to_legacy': to_legacy, 'tidy_only': tidy_only,
    }
    if to_legacy:
        parameters['label'] = 'ASV'

    if sum((to_md5, to_legacy, tidy_only)) != 1:
        # This wording establishes that TypeMap rejects the combination,
        # rather than relying solely on the function's defensive check.
        with pytest.raises(ValueError, match='No solution'):
            fix_labels_action(**label_artifacts, **parameters)
        return

    results = fix_labels_action(**label_artifacts, **parameters)
    results.fixed_feature_table.validate(level='max')
    results.fixed_rep_seqs.validate(level='max')
    assert str(results.fixed_feature_table.type) == 'FeatureTable[Frequency]'
    assert str(results.fixed_rep_seqs.type) == 'FeatureData[Sequence]'
    table = results.fixed_feature_table.view(biom.Table)
    sequences = results.fixed_rep_seqs.view(pd.Series)

    if to_md5:
        expected_ids = [hashlib.md5(sequence).hexdigest()
                        for sequence in [b'ACGT', b'AAAA']]
    elif to_legacy:
        expected_ids = ['ASV1', 'ASV2']
    else:
        expected_ids = ['high', 'low']

    assert list(table.ids(axis='observation')) == expected_ids
    assert list(sequences.index) == expected_ids
    assert list(table.ids(axis='sample')) == ['sample-a', 'sample-b']
    np.testing.assert_array_equal(table.matrix_data.toarray(), [[2, 3], [1, 0]])
    assert [str(sequence) for sequence in sequences] == ['ACGT', 'AAAA']
    assert [sequence.metadata['id'] for sequence in sequences] == expected_ids


def test_registered_defaults_require_a_mode(fix_labels_action, label_artifacts):
    with pytest.raises(ValueError, match='No solution'):
        fix_labels_action(**label_artifacts)


@pytest.mark.parametrize('label', [None, ''])
def test_registered_legacy_requires_label(
        fix_labels_action, label_artifacts, label):
    with pytest.raises(ValueError, match='requires setting a nonempty label'):
        fix_labels_action(**label_artifacts, to_legacy=True, label=label)


@pytest.mark.parametrize('mode', ['to_md5', 'tidy_only'])
def test_registered_label_rejected_outside_legacy(
        fix_labels_action, label_artifacts, mode):
    with pytest.raises(ValueError,
                       match='Set "--p-label" only when "--p-to-legacy"'):
        fix_labels_action(**label_artifacts, **{mode: True}, label='ASV')


@pytest.mark.parametrize('mode', ['to_md5', 'to_legacy', 'tidy_only'])
@pytest.mark.parametrize('zero_sample_ids', [
    ['empty-one'], ['empty-one', 'empty-two'],
])
def test_registered_zero_count_samples_rejected(
        fix_labels_action, label_artifacts, mode, zero_sample_ids):
    table = label_artifacts['feature_table'].view(biom.Table)
    data = np.column_stack([
        table.matrix_data.toarray(),
        np.zeros((table.shape[0], len(zero_sample_ids))),
    ])
    zero_sample_table = biom.Table(
        data, table.ids(axis='observation'),
        [*table.ids(axis='sample'), *zero_sample_ids])
    feature_table = Artifact.import_data('FeatureTable[Frequency]',
                                         zero_sample_table)
    rep_seqs = label_artifacts['rep_seqs']
    # Zero-count samples are format-valid, so rejection must occur in the action.
    feature_table.validate(level='max')
    rep_seqs.validate(level='max')
    original_sequences = rep_seqs.view(pd.Series)

    with pytest.raises(ValueError) as error:
        fix_labels_action(feature_table=feature_table, rep_seqs=rep_seqs,
                          **{mode: True},
                          label='ASV' if mode == 'to_legacy' else None)

    message = str(error.value)
    for sample_id in zero_sample_ids:
        assert sample_id in message
    assert 'sample-a' not in message
    assert 'sample-b' not in message
    assert feature_table.view(biom.Table) == zero_sample_table
    sequences = rep_seqs.view(pd.Series)
    assert list(sequences.index) == list(original_sequences.index)
    assert [str(sequence) for sequence in sequences] == [
        str(sequence) for sequence in original_sequences]
