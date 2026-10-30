"""Regression coverage for Illumina filtering, features and statistics."""

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


@pytest.mark.parametrize('use_vsearch', [False, True])
def test_empty_unoise_output_stops_before_chimera_detection(
        tmp_path, monkeypatch, use_vsearch):
    commands = []

    def empty_output(cmd, verbose):
        commands.append(cmd)
        output = 'vsearch_amps.fasta' if use_vsearch else 'amps.fasta'
        (tmp_path / output).write_text('')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', empty_output)

    with pytest.raises(ValueError, match='No features survived') as error:
        pipeline._unoise_cli(str(tmp_path), min_size=100000,
                             use_vsearch=use_vsearch, verbose=False)

    assert 'min_size=100000' in str(error.value)
    assert len(commands) == 1


def test_empty_vsearch_nonchimera_output_stops_pipeline(tmp_path, monkeypatch):
    commands = []

    def output(cmd, verbose):
        commands.append(cmd)
        if '--cluster_unoise' in cmd:
            (tmp_path / 'vsearch_amps.fasta').write_text('>amp1\nACGT\n')
        else:
            assert '--uchime3_denovo' in cmd
            (tmp_path / 'zotus.fasta').write_text('')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', output)

    with pytest.raises(ValueError, match='No features survived'):
        pipeline._unoise_cli(str(tmp_path), use_vsearch=True, verbose=False)

    assert len(commands) == 2


def test_all_chimera_amplicons_stop_before_empty_reference_is_used(tmp_path):
    (tmp_path / 'amps.fasta').write_text(
        '>amp1;amptype=chimera;size=8\nACGT\n')

    with pytest.raises(ValueError, match='No features survived'):
        pipeline._split_zotu_chimera(str(tmp_path), verbose=False)

    assert (tmp_path / 'chimeras.fasta').stat().st_size > 0
    assert (tmp_path / 'zotus.fasta').stat().st_size == 0


def test_empty_uparse_output_stops_before_parsing_empty_report(
        tmp_path, monkeypatch):
    def empty_output(cmd, verbose):
        (tmp_path / 'otus.fasta').write_text('')
        (tmp_path / 'uparse_out.tsv').write_text('')
        return []

    def unexpected_read(*args, **kwargs):
        pytest.fail('An empty UPARSE report must not be parsed')

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', empty_output)
    monkeypatch.setattr(pipeline.pd, 'read_csv', unexpected_read)

    with pytest.raises(ValueError, match='No features survived') as error:
        pipeline._uparse_cli(str(tmp_path), min_size=100000, verbose=False)

    assert 'min_size=100000' in str(error.value)


@pytest.mark.parametrize('use_vsearch', [False, True])
@pytest.mark.parametrize('stage', ['cluster_input', 'cluster_output',
                                   'zotu_mapping', 'otu_mapping'])
def test_empty_feature_files_never_reach_downstream_commands(
        tmp_path, monkeypatch, use_vsearch, stage):
    commands = []
    (tmp_path / 'zotus.fasta').write_text(
        '>zotu1\nACGT\n' if stage == 'cluster_output' else '')
    (tmp_path / 'otus.fasta').write_text('')

    def output(cmd, verbose):
        commands.append(cmd)
        assert stage == 'cluster_output', 'Empty reference reached backend'
        (tmp_path / 'otus.fasta').write_text('')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', output)

    with pytest.raises(ValueError, match='No features survived'):
        if stage.startswith('cluster_'):
            pipeline._cluster_zotus_cli(str(tmp_path), use_vsearch=use_vsearch,
                                       verbose=False)
        else:
            feature_type = 'zOTU' if stage == 'zotu_mapping' else 'OTU'
            pipeline._build_feature_tab_cli(
                str(tmp_path), feature_type=feature_type,
                use_vsearch=use_vsearch, verbose=False)

    assert len(commands) == (1 if stage == 'cluster_output' else 0)


@pytest.mark.parametrize('use_vsearch', [False, True])
@pytest.mark.parametrize('alpha', [2.0, 1.5])
def test_nonempty_unoise_output_uses_backend_specific_alpha_option(
        tmp_path, monkeypatch, use_vsearch, alpha):
    (tmp_path / 'dereped.fasta').write_text('>seq1;size=20\nACGT\n')
    commands = []

    def output(cmd, verbose):
        commands.append(cmd)
        if '-ampout' in cmd:
            (tmp_path / 'amps.fasta').write_text('>amp1\nACGT\n')
        elif '--centroids' in cmd:
            (tmp_path / 'vsearch_amps.fasta').write_text('>amp1\nACGT\n')
        else:
            assert '--nonchimeras' in cmd
            (tmp_path / 'zotus.fasta').write_text('>zotu1\nACGT\n')
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', output)

    pipeline._unoise_cli(str(tmp_path), unoise_alpha=alpha,
                         use_vsearch=use_vsearch, verbose=False)

    denoising_cmd = commands[0]
    expected_option = '--unoise_alpha' if use_vsearch else '-alpha'
    if alpha == 2.0:
        assert expected_option not in denoising_cmd
    else:
        assert denoising_cmd[denoising_cmd.index(expected_option) + 1] == str(alpha)
    assert '--alpha' not in denoising_cmd
    assert len(commands) == (2 if use_vsearch else 1)
    assert not (tmp_path / 'dereped.fasta').exists()


def test_surviving_nonchimera_feature_keeps_split_counts(tmp_path):
    (tmp_path / 'amps.fasta').write_text(
        '>amp1;size=20\nACGT\n>amp2;amptype=chimera;size=8\nTGCA\n')

    assert pipeline._split_zotu_chimera(str(tmp_path), verbose=False) == (2, 1)
    assert (tmp_path / 'zotus.fasta').stat().st_size > 0
    assert (tmp_path / 'chimeras.fasta').stat().st_size > 0


def test_missing_feature_file_retains_file_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        pipeline._check_nonempty_features(
            str(tmp_path / 'missing.fasta'), 'denoising', min_size=8)


def test_nonempty_malformed_fasta_retains_parser_error(tmp_path):
    (tmp_path / 'amps.fasta').write_text('>feature_without_sequence\n')

    with pytest.raises(pipeline.skbio.io.FASTAFormatError):
        pipeline._split_zotu_chimera(str(tmp_path), verbose=False)
