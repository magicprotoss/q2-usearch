"""Regression coverage for Illumina filtering and per-sample statistics."""

import numpy as np
import pandas as pd
import pytest

from q2_usearch._illumina_pipeline import _percent_of_input
from q2_usearch import _illumina_pipeline as pipeline


def test_blank_sample_preserved_without_changing_other_percentages():
    # Pooling stores input counts as objects; left merges omit blank counts.
    totals = pd.Series([100, 0, 50],
                       index=['sample-a', 'BLANK', 'sample-b'], dtype=object)
    counts = pd.Series([75, np.nan, 10], index=totals.index)
    original_counts = counts.copy()
    original_totals = totals.copy()

    percentages = _percent_of_input(counts, totals)

    pd.testing.assert_series_equal(
        percentages, pd.Series([75.0, 0.0, 20.0], index=totals.index))
    without_blank = _percent_of_input(counts.drop('BLANK'),
                                      totals.drop('BLANK'))
    pd.testing.assert_series_equal(percentages.drop('BLANK'), without_blank)
    assert np.isfinite(percentages).all()
    pd.testing.assert_series_equal(counts, original_counts)
    pd.testing.assert_series_equal(totals, original_totals)


def test_zero_retained_reads_and_missing_mapping_rows_report_zero():
    totals = pd.Series([100, 50], index=['filtered-out', 'unmapped'],
                       dtype=object)
    counts = pd.Series([0, np.nan], index=totals.index)

    pd.testing.assert_series_equal(
        _percent_of_input(counts, totals),
        pd.Series([0.0, 0.0], index=totals.index))


def test_unknown_input_count_is_distinct_from_blank_sample():
    totals = pd.Series([0, np.nan], index=['BLANK', 'unknown'], dtype=object)
    counts = pd.Series([np.nan, 0], index=totals.index)

    percentages = _percent_of_input(counts, totals)

    assert percentages.loc['BLANK'] == 0.0
    assert pd.isna(percentages.loc['unknown'])


@pytest.mark.parametrize('invalid_column', ['counts', 'totals'])
def test_malformed_counts_are_not_silently_treated_as_zero(invalid_column):
    values = {'counts': pd.Series([1]), 'totals': pd.Series([10])}
    values[invalid_column] = pd.Series(['invalid'])

    with pytest.raises(ValueError):
        _percent_of_input(**values)


@pytest.mark.parametrize('use_vsearch', [False, True])
@pytest.mark.parametrize('settings', [{'trunc_right': 10000}, {'max_ee': 0}])
def test_no_reads_passed_filter_has_actionable_error(
        tmp_path, monkeypatch, use_vsearch, settings):
    def empty_filter_output(cmd, verbose):
        (tmp_path / 'filtered.fasta').write_text('')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', empty_filter_output)

    with pytest.raises(ValueError, match='No reads passed the quality filter') as error:
        pipeline._quality_control_cli(
            str(tmp_path), use_vsearch=use_vsearch, verbose=False, **settings)

    assert 'across all samples' in str(error.value)
    assert 'trunc_len=' in str(error.value)
    assert 'max_ee=' in str(error.value)


@pytest.mark.parametrize('use_vsearch', [False, True])
def test_surviving_reads_keep_per_sample_filter_counts(
        tmp_path, monkeypatch, use_vsearch):
    (tmp_path / 'merged.fastq').write_text('')

    def surviving_filter_output(cmd, verbose):
        (tmp_path / 'filtered.fasta').write_text(
            '>sample_a.1\nACGT\n>sample_a.2\nACGT\n>sample_b.1\nTGCA\n')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', surviving_filter_output)

    stats = pipeline._quality_control_cli(
        str(tmp_path), use_vsearch=use_vsearch, verbose=False)

    expected = pd.DataFrame(
        {'reads_passed_filter': [2, 1]},
        index=pd.Index(['sample_a', 'sample_b'], name='sample-id'))
    pd.testing.assert_frame_equal(stats, expected)
