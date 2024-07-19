import pandas as pd
import numpy as np
import re
import os
import skbio
import subprocess
import tempfile
import hashlib


def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


def _get_input_seqs_ids_and_dump_to_fasta(working_dir, query_se):
    with open(os.path.join(working_dir, 'query.fasta'), 'wt') as fh:
        for index, item in query_se.items():
            item.write(fh, format="fasta", max_width=80)
    empty_df_w_input_seqs_labs = pd.DataFrame(index=query_se.index)
    return empty_df_w_input_seqs_labs


def _split_tax_into_ranks(tax, sep):
    tax_lst = str(tax).split(sep)
    return tax_lst


def _split_tax_into_ranks_and_get_max_levels(tax_df_in, sep):
    # need to confirm if sintax accepts non-7-rank systems
    tax_rank_split_se = tax_df_in['Taxon'].apply(_split_tax_into_ranks, sep=sep)
    max_levels = tax_rank_split_se.apply(len).max()
    if max_levels > 7:
        raise KeyError('according to the doc, sintax only supports up to 7 levels')
    return max_levels, tax_rank_split_se


def _replace_q2_split_w_usearch_split_and_remove_leading_trailing_blanks(rank_in):
    rank_out = re.sub(r"(?<=\b[dpcofgs])\w*__", ':', str(rank_in).strip())
    return rank_out


def _replace_non_7bit_ascii_chars(rank_in):
    rank_out = rank_in

    # let's KISS here...
    if re.search(r'[^A-Za-z0-9_.]', rank_in):
        if re.match(r'^[kdpcofgs]', str(rank_in).strip()):
            rank_out = rank_in[0] + ':' + \
                hashlib.md5(str(rank_in).encode('utf-8')).hexdigest()
        else:
            rank_out = np.nan
    return rank_out


def _detect_empty_ph_ranks(rank_in):
    result = False
    if rank_in == np.nan:
        result = True
    else:
        pattern = r"^[dpcofgs]:$"
        result = bool(re.match(pattern, rank_in))
    return result


def _join_levels_for_usearch(tax_series_in):
    # This pile of 💩 is here for one reason:
    # if a non authoritive database contains a place holder in a parent rank and also have a non-empty child rank
    usearch_tax_anno_str = ';tax='
    for level, (index, rank) in enumerate(tax_series_in.items(), 1):
        detection_res = _detect_empty_ph_ranks(rank)
        if detection_res:
            break
        else:
            if level == 1:
                usearch_tax_anno_str = usearch_tax_anno_str + rank
            else:
                usearch_tax_anno_str = usearch_tax_anno_str + ',' + rank
    usearch_tax_anno_str = usearch_tax_anno_str + ';'
    return usearch_tax_anno_str


def _make_tmp_tax_mapping_df(tax_df_in):

    tax_df_out = tax_df_in.copy()

    max_level, tax_rank_split_se = _split_tax_into_ranks_and_get_max_levels(
        tax_df_in, ';')

    ori_tax_cols_lst = ['q2_' + 'level' + '_' + str(i) for i in range(1, max_level + 1)]
    u_tax_cols_lst = ['usearch_' + 'level' + '_' +
                      str(i) for i in range(1, max_level + 1)]

    # 💩 super slow but works
    tax_df_out[ori_tax_cols_lst] = tax_rank_split_se.apply(pd.Series)

    tax_df_out[u_tax_cols_lst] = tax_df_out[ori_tax_cols_lst]

    for u_level in u_tax_cols_lst:
        tax_df_out[u_level] = tax_df_out[u_level].apply(
            _replace_q2_split_w_usearch_split_and_remove_leading_trailing_blanks)
        tax_df_out[u_level] = tax_df_out[u_level].apply(_replace_non_7bit_ascii_chars)
    tax_df_out['usearch_tax'] = tax_df_out[u_tax_cols_lst].apply(
        _join_levels_for_usearch, axis=1)

    return tax_df_out


def _convert_q2_seqs_and_taxa_to_utax(working_dir, reference_reads, reference_taxonomy, verbose):
    ref_reads_se = reference_reads
    ref_reads_se.name = 'Seqs'
    ref_reads_se.index.name = 'Feature ID'
    ref_taxa_df = reference_taxonomy
    # check if dumping tax_df to pickle is nessesary with low spec pcs
    # silva 138.1 only took 72m mem, no need here
    if verbose:
        print("Building usearch compatible fasta db file, this could take a while...")

    tmp_taxa_map_df = _make_tmp_tax_mapping_df(ref_taxa_df)

    op_fa = os.path.join(working_dir, 'ref_seqs_tax.fa')
    with open(op_fa, 'wt') as fh:
        for index, item in ref_reads_se.items():
            tax_info = tmp_taxa_map_df.at[index, 'usearch_tax']
            seq_to_dump = skbio.DNA(str(item).upper(), metadata={
                                    'id': index + tax_info})
            seq_to_dump.write(fh, format="fasta", max_width=80)

    return tmp_taxa_map_df


def _build_udb():
    # seemed unnessasary, sintax builds one on the fly, and tmp dirs don't presist in a q2 pipeline
    pass


def _run_sintax(working_dir, query_seqs_fp, strand, threads, verbose):
    # build sintax command
    cmd = ['usearch', '-sintax', query_seqs_fp, '-db', os.path.join(
        working_dir, 'ref_seqs_tax.fa'), '-tabbedout', os.path.join(working_dir, 'sintax.tsv')]

    if strand == 'plus':
        cmd += ['-strand', 'plus']
    else:
        cmd += ['-strand', 'both']

    if threads != 1:
        cmd += ['-threads', str(threads)]

    run_commands([cmd])


def _rm_conf_value_and_trim_fp_ranks(x, cut_off=float):
    if x == None:
        x_str = np.nan
    else:
        sintan_conf_col_loci = x.rfind('(')
        x_str = x[:sintan_conf_col_loci]
        x_conf = float(x[sintan_conf_col_loci + 1:].replace(')', ''))
        if x_conf < cut_off:
            x_str = np.nan
    return x_str


def _get_conf_value_deepest_rank(se_in, cut_off=float):

    conf = np.nan
    for index, item in se_in.items():
        sintan_conf_col_loci = item.rfind('(')
        conf_item = float(item[sintan_conf_col_loci + 1:].replace(')', ''))
        if conf_item < cut_off:
            break
        conf = conf_item

    return conf


def _split_utax_and_get_conf_lr(usearch_tax_df_in, confidence):

    tax_rank_split_df = pd.DataFrame(index=usearch_tax_df_in.index)

    max_level, usearch_tax_rank_split_df = _split_tax_into_ranks_and_get_max_levels(
        usearch_tax_df_in, ',')

    u_tax_cols_lst = ['usearch_' + 'level' + '_' +
                      str(i) for i in range(1, max_level + 1)]

    # 💩 super slow but works
    tax_rank_split_df[u_tax_cols_lst] = usearch_tax_rank_split_df.apply(pd.Series)

    id_conf_df = pd.DataFrame(index=usearch_tax_df_in.index)

    for index, row in tax_rank_split_df.iterrows():
        id_conf_df.at[index, 'Confidence'] = _get_conf_value_deepest_rank(
            row, confidence)

    tax_rank_split_df = tax_rank_split_df.applymap(
        _rm_conf_value_and_trim_fp_ranks, cut_off=confidence)

    return tax_rank_split_df, id_conf_df


def _map_utax_to_q2_tax(tax_rank_split_in, taxa_map_df):
    # split input and map into sep dfs, left join and concat back
    levels = len(tax_rank_split_in.columns)
    tax_rank_split_mapped_to_q2_tax = tax_rank_split_in.copy()
    for i in range(1, levels + 1):
        key = "level_" + str(i)
        taxa_map_df_sub_level = taxa_map_df[['usearch_' + key, 'q2_' + key]]
        taxa_map_df_sub_level_uni = taxa_map_df_sub_level.drop_duplicates().dropna()
        taxa_map_df_sub_level_uni_dict = taxa_map_df_sub_level_uni.set_index(
            'usearch_' + key).iloc[:, 0].to_dict()
        tax_rank_split_mapped_to_q2_tax['usearch_' + key] = tax_rank_split_in['usearch_' + key].replace(
            taxa_map_df_sub_level_uni_dict)
    return tax_rank_split_mapped_to_q2_tax


def _join_q2_tax(q2_tax_rank_split_in):
    q2_tax = pd.DataFrame()
    for index, row in q2_tax_rank_split_in.iterrows():
        tax_str = '; '.join(row.dropna().values)
        if len(tax_str) == 0:
            tax_str = 'Unclassified'
        tmp_df = pd.DataFrame({'Taxon': {index: tax_str}})
        q2_tax = pd.concat([q2_tax, tmp_df])
    q2_tax.index.name = 'Feature ID'

    return q2_tax


def _comp_plus_minus_res_and_opt_final_res(empty_df_w_input_seqs_labs, q2_taxs, q2_tax_rank_splits, id_conf_dfs):
    q2_tax = pd.DataFrame(index=empty_df_w_input_seqs_labs.index,
                          columns=['Taxon', 'Confidence'])
    # purge 💩 here later...
    for index in empty_df_w_input_seqs_labs.index.to_list():
        plus_hits_index = q2_taxs['plus'].loc[q2_taxs['plus']
                                              ['Taxon'] != 'Unclassified'].index
        minus_hits_index = q2_taxs['minus'].loc[q2_taxs['minus']
                                                ['Taxon'] != 'Unclassified'].index
        plus_hit = index in plus_hits_index
        minus_hit = index in minus_hits_index
        if plus_hit and minus_hit:
            plus_depth = len(q2_tax_rank_splits['plus'].loc[index, :].dropna())
            minus_depth = len(q2_tax_rank_splits['minus'].loc[index, :].dropna())
            plus_conf = id_conf_dfs['plus'].at[index, 'Confidence']
            minus_conf = id_conf_dfs['minus'].at[index, 'Confidence']
            if plus_depth > minus_depth:
                q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                conf_assign = plus_conf
            elif plus_depth < minus_depth:
                q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
                conf_assign = minus_conf
            else:
                if plus_conf > minus_conf:
                    q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                    conf_assign = plus_conf
                elif plus_conf < minus_conf:
                    q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
                    conf_assign = minus_conf
                else:  # 不会真出回文吧
                    q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                    conf_assign = plus_conf
        elif plus_hit:
            q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
            conf_assign = id_conf_dfs['plus'].at[index, 'Confidence']
        elif minus_hit:
            q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
            conf_assign = id_conf_dfs['minus'].at[index, 'Confidence']
        else:
            q2_tax.at[index, 'Taxon'] = 'Unclassified'
            conf_assign = np.nan

        q2_tax.at[index, 'Confidence'] = conf_assign

    return q2_tax


def _collect_sintax_anno_to_q2_anno(working_dir, taxa_map_df, empty_df_w_input_seqs_labs, strand, confidence, verbose):

    # read sintax res into pd.DataFrame
    usearch_tax = pd.read_csv(os.path.join(
        working_dir, 'sintax.tsv'), sep='\t', header=None, index_col=0)

    usearch_tax.index.name = 'Feature ID'

    usearch_tax = usearch_tax.iloc[:, 0:2]

    usearch_tax.columns = ['Taxon', 'Strand']

    # rewrite in a less 💩 way here later...
    if strand != 'plus':
        usearch_tax_both = {}
        usearch_tax_both['plus'] = usearch_tax.loc[usearch_tax['Strand']
                                                   == '+', 'Taxon'].to_frame()
        usearch_tax_both['minus'] = usearch_tax.loc[usearch_tax['Strand']
                                                    == '-', 'Taxon'].to_frame()
        usearch_tax_rank_split_dfs = {}
        id_conf_dfs = {}
        q2_tax_rank_splits = {}
        q2_taxs = {}
        for key, value in usearch_tax_both.items():

            usearch_tax_rank_split_dfs[key], id_conf_dfs[key] = _split_utax_and_get_conf_lr(
                value, confidence)

            q2_tax_rank_splits[key] = _map_utax_to_q2_tax(
                usearch_tax_rank_split_dfs[key], taxa_map_df)

            q2_taxs[key] = _join_q2_tax(q2_tax_rank_splits[key])

        q2_tax = _comp_plus_minus_res_and_opt_final_res(
            empty_df_w_input_seqs_labs, q2_taxs, q2_tax_rank_splits, id_conf_dfs)

        # purge 💩 here later...
        for index in empty_df_w_input_seqs_labs.index:
            plus_hit = index in q2_taxs['plus'].index
            minus_hit = index in q2_taxs['minus'].index
            if plus_hit and minus_hit:
                plus_depth = len(q2_tax_rank_splits['plus'].loc[index, :].dropna())
                minus_depth = len(q2_tax_rank_splits['minus'].loc[index, :].dropna())
                plus_conf = id_conf_dfs['plus'].at[index, 'Confidence']
                minus_conf = id_conf_dfs['minus'].at[index, 'Confidence']
                if plus_depth > minus_depth:
                    q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                    conf_assign = plus_conf
                elif plus_depth < minus_depth:
                    q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
                    conf_assign = minus_conf
                else:
                    if plus_conf > minus_conf:
                        q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                        conf_assign = plus_conf
                    elif plus_conf < minus_conf:
                        q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
                        conf_assign = minus_conf
                    else:  # 不会真出回文吧
                        q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                        conf_assign = plus_conf
            elif plus_hit:
                q2_tax.at[index, 'Taxon'] = q2_taxs['plus'].at[index, 'Taxon']
                conf_assign = id_conf_dfs['plus'].at[index, 'Confidence']
            elif minus_hit:
                q2_tax.at[index, 'Taxon'] = q2_taxs['minus'].at[index, 'Taxon']
                conf_assign = id_conf_dfs['minus'].at[index, 'Confidence']
            else:
                q2_tax.at[index, 'Taxon'] = 'Unclassified'
                conf_assign = np.nan

            q2_tax.at[index, 'Confidence'] = conf_assign

    else:
        usearch_tax = usearch_tax['Taxon'].to_frame()

        usearch_tax_rank_split_df, id_conf_df = _split_utax_and_get_conf_lr(
            usearch_tax, confidence)

        q2_tax_rank_split = _map_utax_to_q2_tax(usearch_tax_rank_split_df, taxa_map_df)

        q2_tax = _join_q2_tax(q2_tax_rank_split)

        q2_tax = pd.merge(q2_tax, id_conf_df, left_index=True,
                          right_index=True, how='inner')

        q2_tax = pd.merge(empty_df_w_input_seqs_labs, q2_tax,
                          left_index=True, right_index=True, how='left')

        for index, row in q2_tax.iterrows():
            if row['Taxon'] is None:
                q2_tax.at[index, 'Taxon'] == 'Unclassified'

    return q2_tax


def sintax(query: pd.Series,
           reference_reads: pd.Series,
           reference_taxonomy: pd.DataFrame,
           # limited test suggest it's common for sintax to report a better match in rev-comp using plus only 16s as input
           # maybe throw in orinet as a precaution? warn user?
           strand: str = 'plus',
           threads: str = "auto",
           confidence: float = 0.8
           ) -> (pd.DataFrame):

    verbose = True

    if threads == "auto":
        threads = os.cpu_count() - 3

    with tempfile.TemporaryDirectory() as usearch_wd:

        empty_df_w_input_seqs_labs = _get_input_seqs_ids_and_dump_to_fasta(
            usearch_wd, query)

        taxa_map_df = _convert_q2_seqs_and_taxa_to_utax(
            usearch_wd, reference_reads, reference_taxonomy, verbose)

        query_fp = os.path.join(usearch_wd, 'query.fasta')

        _run_sintax(usearch_wd, query_fp, strand, threads, verbose)

        classification = _collect_sintax_anno_to_q2_anno(
            usearch_wd, taxa_map_df, empty_df_w_input_seqs_labs, strand, confidence, verbose)

    return classification
