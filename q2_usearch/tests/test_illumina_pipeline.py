"""Regression coverage for Illumina filtering, features and statistics."""

import csv
import gzip
from pathlib import Path

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


@pytest.mark.parametrize('debug', [False, True], ids=['seqkit', 'skbio'])
@pytest.mark.parametrize('sample_ids', [
    ['sample_a', 'sample_b', 'safe123'],
    ['0001', '1', '01'],
    ['NA', 'NULL', 'nan'],
    ['sample.a#one', 'sample-b', '样本'],
    ['S2', 'S1', 'S3'],
], ids=['safe', 'numeric', 'na-like', 'special', 'temporary-like'])
def test_pooling_always_uses_temporary_ids_and_preserves_originals(
        tmp_path, monkeypatch, debug, sample_ids):
    input_dir = tmp_path / 'input'
    work_dir = tmp_path / 'work'
    input_dir.mkdir()
    work_dir.mkdir()
    # The third sample is blank; its mapping and zero count must survive.
    records = [
        [('old-a', 'ACGT', 'I?5!'), ('old-b', 'TGCA', '!5?I')],
        [('old-c', 'AGCT', '5678')],
        [],
    ]
    originals = {}
    with (input_dir / 'MANIFEST').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample-id', 'filename', 'direction'])
        for number, (sample_id, reads) in enumerate(zip(sample_ids, records)):
            filename = f'input-{number}.fastq.gz'
            fastq = ''.join(f'@{rid} annotation\n{seq}\n+\n{qual}\n'
                            for rid, seq, qual in reads)
            with gzip.open(input_dir / filename, 'wt') as fastq_handle:
                fastq_handle.write(fastq)
            originals[filename] = (input_dir / filename).read_bytes()
            writer.writerow([sample_id, filename, 'forward'])
    (input_dir / 'metadata.yml').write_text('phred-offset: 33\n')
    commands = []

    def fake_seqkit(cmd, verbose=False):
        commands.append(cmd)
        assert cmd[0] == 'seqkit'
        output = Path(cmd[cmd.index('-o') + 1])
        if cmd[1] == 'replace':
            input_path = Path(cmd[cmd.index('-i') + 1])
            replacement = cmd[cmd.index('-r') + 1]
            with gzip.open(input_path, 'rt') as handle:
                lines = handle.read().splitlines()
            for number, start in enumerate(range(0, len(lines), 4), 1):
                lines[start] = '@' + replacement.replace('{nr}', str(number))
            output.write_text('\n'.join(lines) + ('\n' if lines else ''))
        else:
            assert cmd[1] == 'stats'
            with gzip.open(cmd[3], 'rt') as handle:
                count = len(handle.read().splitlines()) // 4
            pd.DataFrame({'num_seqs': [count]}).to_csv(
                output, sep='\t', index=False)
        return []

    monkeypatch.setattr(pipeline, 'py_to_cli_interface', fake_seqkit)

    stats = pipeline._pool_samples(str(input_dir), str(work_dir),
                                   debug=debug, verbose=False)

    assert stats.index.name == 'sample-id'
    assert stats.index.tolist() == ['S1', 'S2', 'S3']
    assert stats['original_sample_id'].tolist() == sample_ids
    assert stats['prior_to_maxee_filt'].tolist() == [2, 1, 0]
    assert 'fixed_sample_id' not in stats.columns
    pooled_reads = list(pipeline.skbio.io.read(
        str(work_dir / 'merged.fastq'), format='fastq',
        variant='illumina1.8'))
    assert [read.metadata['id'] for read in pooled_reads] == [
        'S1.1', 'S1.2', 'S2.1']
    expected_reads = [read for sample_reads in records for read in sample_reads]
    assert [str(read) for read in pooled_reads] == [
        sequence for _, sequence, _ in expected_reads]
    for actual, (_, _, quality) in zip(pooled_reads, expected_reads):
        np.testing.assert_array_equal(
            actual.positional_metadata['quality'].to_numpy(),
            [ord(character) - 33 for character in quality])
    for filename, contents in originals.items():
        assert (input_dir / filename).read_bytes() == contents
    if debug:
        assert commands == []
    else:
        assert len(commands) == 6
        assert [cmd[cmd.index('-r') + 1] for cmd in commands
                if cmd[1] == 'replace'] == [
                    'S1.{nr}', 'S2.{nr}', 'S3.{nr}']


@pytest.mark.parametrize('action_name,use_vsearch', [
    ('denoise_no_primer_pooled', False),
    ('cluster_no_primer_pooled', False),
    ('denoise_then_cluster_no_primer_pooled', False),
    ('denoise_no_primer_pooled', True),
    ('denoise_then_cluster_no_primer_pooled', True),
])
def test_all_actions_restore_original_sample_ids_in_table_and_statistics(
        monkeypatch, action_name, use_vsearch):
    # Swapping S1/S2 must change labels without swapping the associated data.
    original_ids = ['S2', 'S1', '0001', 'NA', 'BLANK']
    temporary_ids = ['S1', 'S2', 'S3', 'S4', 'S5']
    index = pd.Index(temporary_ids, name='sample-id')
    input_stats = pd.DataFrame({
        'prior_to_maxee_filt': [10, 20, 30, 40, 0],
        'original_sample_id': original_ids,
    }, index=index)
    filtered = pd.DataFrame({'reads_passed_filter': [9, 18, 27, 36]},
                            index=index[:4])
    table_data = np.array([[5, 6, 7, 8], [0, 1, 2, 3]])
    table = pipeline.biom.Table(table_data, ['f1', 'f2'], temporary_ids[:4])
    representatives = pd.Series(['ACGT', 'TGCA'], index=['f1', 'f2'])
    mapped_column = ('reads_mapped_to_zotus'
                     if action_name == 'denoise_no_primer_pooled'
                     else 'reads_mapped_to_otus')
    mapped = pd.DataFrame({mapped_column: table_data.sum(axis=0)},
                          index=index[:4])
    chimeras = pd.DataFrame({'reads_mapped_to_chimeras': [1]}, index=index[1:2])
    monkeypatch.setattr(pipeline, '_pool_samples',
                        lambda *args, **kwargs: input_stats.copy())
    monkeypatch.setattr(pipeline, '_quality_control_cli',
                        lambda *args, **kwargs: filtered.copy())
    monkeypatch.setattr(pipeline, '_dereplicate_cli',
                        lambda *args, **kwargs: (90, 4, 1))
    monkeypatch.setattr(pipeline, '_unoise_cli',
                        lambda *args, **kwargs: None)
    monkeypatch.setattr(pipeline, '_split_zotu_chimera',
                        lambda *args, **kwargs: (3, 2))
    monkeypatch.setattr(pipeline, '_uparse_cli',
                        lambda *args, **kwargs: (2, 0))
    monkeypatch.setattr(pipeline, '_cluster_zotus_cli',
                        lambda *args, **kwargs: 2)
    monkeypatch.setattr(pipeline, '_build_feature_tab_cli',
                        lambda *args, **kwargs: None)
    monkeypatch.setattr(
        pipeline, '_prep_results_for_artifact_api',
        lambda *args, **kwargs: (table.copy(), representatives.copy(),
                                mapped.copy(), chimeras.copy()))
    kwargs = {} if action_name == 'cluster_no_primer_pooled' else {
        'use_vsearch': use_vsearch}

    actual_table, actual_representatives, actual_metadata = getattr(
        pipeline, action_name)('unused-input', **kwargs)

    assert actual_table.ids(axis='sample').tolist() == original_ids[:4]
    np.testing.assert_array_equal(actual_table.matrix_data.toarray(), table_data)
    pd.testing.assert_series_equal(actual_representatives, representatives)
    stats = actual_metadata.to_dataframe()
    assert stats.index.tolist() == original_ids
    assert stats['prior_to_maxee_filt'].tolist() == [10, 20, 30, 40, 0]
    assert stats.loc[original_ids[:4], mapped_column].tolist() == [5, 7, 9, 11]
    assert stats.loc['S1', 'reads_mapped_to_chimeras'] == 1
    assert 'original_sample_id' not in stats.columns
    assert 'fixed_sample_id' not in stats.columns
    assert stats.loc['BLANK'].drop('denoise_stats_pooled_mode').eq(0).all()
    assert 'BLANK' not in actual_table.ids(axis='sample')
