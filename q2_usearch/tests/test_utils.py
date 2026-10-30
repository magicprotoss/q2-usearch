import hashlib
import itertools

import biom
import numpy as np
import pandas as pd
import pytest
import skbio
from scipy import sparse

from q2_usearch._utils import fix_labels


@pytest.fixture
def labels_data():
    # Ties deliberately straddle the most abundant feature in input order.
    # The sample named rep_seqs reproduced a column collision in the old method.
    table = biom.Table(
        sparse.csr_matrix([[2, 3, 1], [0, 0, 0], [9, 1, 1], [5, 0, 1]]),
        ['tie-first', 'zero', 'high', 'tie-second'],
        ['rep_seqs', 'sample.two', 'sample-three'])
    sequences = pd.Series(
        [skbio.DNA(seq, metadata={'id': identifier})
         for identifier, seq in [('extra', 'GGGG'), ('tie-second', 'TTTT'),
                                 ('high', 'ACGT'), ('zero', 'CCCC'),
                                 ('tie-first', 'AAAA')]],
        index=['extra', 'tie-second', 'high', 'zero', 'tie-first'])
    return table, sequences


@pytest.mark.parametrize('mode,label', [
    ('tidy_only', None), ('to_legacy', 'ASV'), ('to_legacy', 'OTU_'),
    ('to_legacy', 'asv'), ('to_legacy', 'AsV_09'),
    ('to_legacy', '123'), ('to_legacy', '_'),
    # This prefix is literal in legacy mode, not a special mode selector.
    ('to_legacy', 'md5'), ('to_md5', None),
])
def test_tidy_order_counts_and_synchronized_ids(labels_data, mode, label):
    table, sequences = labels_data
    original_table = table.copy()
    original_sequence_ids = list(sequences.index)
    original_sequence_text = [str(seq) for seq in sequences]
    original_metadata = [dict(seq.metadata) for seq in sequences]

    result_table, result_sequences = fix_labels(
        table, sequences, **{mode: True}, label=label)

    ordered_sequence_text = ['ACGT', 'AAAA', 'TTTT']
    if mode == 'to_legacy':
        expected_ids = [f'{label}{i}' for i in range(1, 4)]
    elif mode == 'to_md5':
        expected_ids = [hashlib.md5(seq.encode('utf-8')).hexdigest()
                        for seq in ordered_sequence_text]
    else:
        expected_ids = ['high', 'tie-first', 'tie-second']

    assert isinstance(result_table, biom.Table)
    assert sparse.issparse(result_table.matrix_data)
    assert isinstance(result_sequences, pd.Series)
    assert list(result_table.ids(axis='observation')) == expected_ids
    assert list(result_sequences.index) == expected_ids
    assert list(result_table.ids(axis='sample')) == [
        'rep_seqs', 'sample.two', 'sample-three']
    np.testing.assert_array_equal(result_table.matrix_data.toarray(),
                                  [[9, 1, 1], [2, 3, 1], [5, 0, 1]])
    np.testing.assert_array_equal(result_table.sum(axis='sample'),
                                  original_table.sum(axis='sample'))
    assert [str(seq) for seq in result_sequences] == ordered_sequence_text
    assert all(isinstance(seq, skbio.DNA) for seq in result_sequences)
    assert [seq.metadata['id'] for seq in result_sequences] == expected_ids

    # Relabeling must not mutate callers' table, index, or sequence metadata.
    assert table == original_table
    assert list(sequences.index) == original_sequence_ids
    assert [str(seq) for seq in sequences] == original_sequence_text
    assert [seq.metadata for seq in sequences] == original_metadata


@pytest.mark.parametrize('modes', [
    modes for modes in itertools.product([False, True], repeat=3)
    if sum(modes) != 1
])
def test_exactly_one_mode_required(labels_data, modes):
    table, sequences = labels_data
    with pytest.raises(ValueError):
        fix_labels(table, sequences, to_md5=modes[0], to_legacy=modes[1],
                   tidy_only=modes[2], label='ASV' if modes[1] else None)


@pytest.mark.parametrize('label', [
    None, '', ' ', '\t\n', 'AS V', 'ASV\tfoo', 'ASV\nfoo',
    'ASV-', 'ASV.', 'ASV/', 'ASV#', 'ASV;',
    'é', '特征', '１２３', 'ASV🧬',
])
def test_legacy_requires_nonempty_ascii_alphanumeric_or_underscore_label(
        labels_data, label):
    with pytest.raises(ValueError):
        fix_labels(*labels_data, to_legacy=True, label=label)


@pytest.mark.parametrize('mode', ['to_md5', 'tidy_only'])
@pytest.mark.parametrize('label', ['ASV', ''])
def test_label_rejected_outside_legacy(labels_data, mode, label):
    with pytest.raises(ValueError):
        fix_labels(*labels_data, **{mode: True}, label=label)


@pytest.mark.parametrize('mode', ['to_md5', 'to_legacy', 'tidy_only'])
@pytest.mark.parametrize('missing_id', ['high', 'zero'])
def test_missing_table_sequence_rejected_before_cleanup(
        labels_data, mode, missing_id):
    table, sequences = labels_data
    sequences = sequences.drop(missing_id)
    with pytest.raises(ValueError, match=missing_id):
        fix_labels(table, sequences, **{mode: True},
                   label='ASV' if mode == 'to_legacy' else None)


def test_missing_sequence_error_lists_all_missing_ids(labels_data):
    table, sequences = labels_data
    with pytest.raises(ValueError) as error:
        fix_labels(table, sequences.drop(['high', 'zero']), tidy_only=True)
    assert 'high' in str(error.value)
    assert 'zero' in str(error.value)


def test_md5_uses_uppercase_sequence():
    table = biom.Table(np.array([[1]]), ['old'], ['sample'])
    sequences = pd.Series(['acgt'], index=['old'])
    result_table, result_sequences = fix_labels(table, sequences, to_md5=True)
    expected = hashlib.md5(b'ACGT').hexdigest()
    assert list(result_table.ids(axis='observation')) == [expected]
    assert list(result_sequences.index) == [expected]
    assert str(result_sequences.iloc[0]) == 'ACGT'


def test_md5_rejects_duplicate_retained_sequences(labels_data):
    table, sequences = labels_data
    sequences.loc['tie-first'] = skbio.DNA('ACGT')
    with pytest.raises(ValueError):
        fix_labels(table, sequences, to_md5=True)


def test_md5_ignores_duplicates_in_discarded_sequences(labels_data):
    table, sequences = labels_data
    sequences.loc['zero'] = skbio.DNA('ACGT')
    sequences.loc['extra'] = skbio.DNA('ACGT')
    result_table, result_sequences = fix_labels(table, sequences, to_md5=True)
    assert result_table.shape == (3, 3)
    assert result_sequences.index.is_unique


@pytest.mark.parametrize('mode', ['to_legacy', 'tidy_only'])
def test_duplicate_dna_allowed_when_ids_are_distinct(labels_data, mode):
    table, sequences = labels_data
    sequences.loc['tie-first'] = skbio.DNA('ACGT')
    result_table, result_sequences = fix_labels(
        table, sequences, **{mode: True},
        label='ASV' if mode == 'to_legacy' else None)
    assert result_table.shape == (3, 3)
    assert result_sequences.index.is_unique
    assert [str(seq) for seq in result_sequences].count('ACGT') == 2


def test_duplicate_sequence_ids_rejected(labels_data):
    table, sequences = labels_data
    sequences = pd.concat([sequences, sequences.loc[['high']]])
    with pytest.raises(ValueError):
        fix_labels(table, sequences, tidy_only=True)


def test_all_zero_table_rejected():
    table = biom.Table(np.zeros((1, 2)), ['zero'], ['a', 'b'])
    sequences = pd.Series([skbio.DNA('ACGT')], index=['zero'])
    with pytest.raises(ValueError):
        fix_labels(table, sequences, tidy_only=True)


@pytest.mark.parametrize('mode', ['to_md5', 'to_legacy', 'tidy_only'])
@pytest.mark.parametrize('zero_sample_ids', [
    ['empty-one'], ['empty-one', 'empty-two'],
])
def test_zero_count_samples_rejected_without_mutating_inputs(
        mode, zero_sample_ids):
    sample_ids = ['positive-one', *zero_sample_ids, 'positive-two']
    data = np.array([[2, 3], *[[0, 0] for _ in zero_sample_ids], [1, 4]]).T
    table = biom.Table(sparse.csr_matrix(data), ['f1', 'f2'], sample_ids)
    sequences = pd.Series(
        [skbio.DNA('ACGT', metadata={'id': 'f1'}),
         skbio.DNA('AAAA', metadata={'id': 'f2'})], index=['f1', 'f2'])
    original_table = table.copy()
    original_sequence_ids = list(sequences.index)
    original_sequence_text = [str(sequence) for sequence in sequences]
    original_metadata = [dict(sequence.metadata) for sequence in sequences]

    with pytest.raises(ValueError) as error:
        fix_labels(table, sequences, **{mode: True},
                   label='ASV' if mode == 'to_legacy' else None)

    message = str(error.value)
    for sample_id in zero_sample_ids:
        assert sample_id in message
    assert 'positive-one' not in message
    assert 'positive-two' not in message
    assert table == original_table
    assert list(sequences.index) == original_sequence_ids
    assert [str(sequence) for sequence in sequences] == original_sequence_text
    assert [sequence.metadata for sequence in sequences] == original_metadata


def test_missing_sequence_checked_before_zero_count_samples():
    table = biom.Table(np.array([[1, 0], [0, 0]]),
                       ['present', 'missing'], ['positive-one', 'empty-one'])
    sequences = pd.Series([skbio.DNA('ACGT')], index=['present'])
    with pytest.raises(ValueError) as error:
        fix_labels(table, sequences, tidy_only=True)
    assert 'missing' in str(error.value)
    assert 'empty-one' not in str(error.value)
