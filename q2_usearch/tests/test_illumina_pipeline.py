"""Regression coverage for per-sample Illumina pipeline statistics."""

import numpy as np
import pandas as pd
import pytest

from q2_usearch._illumina_pipeline import _percent_of_input


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
