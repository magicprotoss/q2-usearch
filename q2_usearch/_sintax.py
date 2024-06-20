import os
import time
import pandas as pd
import numpy as np
import re
import tempfile
from q2_types.feature_data import DNAFASTAFormat
import skbio
import subprocess


def run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)

# a udb file as an optional input

def sintax(query: DNAFASTAFormat,
           reference_reads: DNAFASTAFormat,
           taxonomy: pd.DataFrame,
           threads: str = "auto",
           confidence: float = 0.8,
           drop_species: bool = True) -> pd.DataFrame:
    if threads == "auto":
        threads = os.cpu_count() - 1
    return _sintax(query, reference_reads, taxonomy, threads, confidence, drop_species)


def _sintax(query, reference_reads, taxonomy, threads, confidence, drop_species):
    query_fp = str(query)
    reference_reads_fp = str(reference_reads)
    # Construct pandas df from taxonomy file
#    temp_taxonomy_df = taxonomy.reset_index()
    temp_taxonomy_df = taxonomy
    # Fix taxonomy formatting for usearch

    temp_taxonomy_df.loc[:, 'Taxon'] = temp_taxonomy_df.loc[:, 'Taxon'].replace(
        "__", ":", regex=True).replace("; ", ",", regex=True).replace(";", ",", regex=True)
    # usearch really hates these characters... substitute them with something else, then fix them back later
    temp_taxonomy_df.loc[:, 'Taxon'] = temp_taxonomy_df.loc[:, 'Taxon'].replace(
        " ", "+", regex=True).replace("_", "^", regex=True)

    # debug
    temp_taxonomy_df.to_csv(
        '/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintaxtemp_taxonomy_df.csv')
    ########################### bug here ############################

    # silva db ...... why ?
    temp_taxonomy_df.loc[:, 'Taxon'] = temp_taxonomy_df.loc[:, 'Taxon'].replace(
        "\(", "\[", regex=True).replace("\)", "\]", regex=True)
    # For silva db only, sintax assume one child belong to one parent, we need to strip all unclutered taxa
    if temp_taxonomy_df.astype(str).apply(lambda x: x.str.contains('uncultured')).any().any() == True:
        print('\n Looks like you are using silva db, this (spicy) is (events) really (did) not (happen) recommended (https://drive5.com/usearch/manual/silva_preprint_response.html), stripping unclutured taxa... \n')
        # why this won't work ??
        temp_taxonomy_df.loc[:, 'Taxon'].replace(
            'uncultured', '', inplace=True)
    if drop_species == True:
        temp_taxonomy_df[['Taxon_tmp', 'species_to_drop']
                         ] = temp_taxonomy_df['Taxon'].str.split("s:", expand=True)
        temp_taxonomy_df = temp_taxonomy_df.drop(
            columns=['Taxon', 'species_to_drop'])
        temp_taxonomy_df.rename(columns={'Taxon_tmp': 'Taxon'}, inplace=True)
    # temp_taxonomy_df['label'] = temp_taxonomy_df.loc[:, 'Feature ID':'Taxon'].apply(
    #     lambda x: ';tax='.join(x.dropna().astype(str)), axis=1)
    # Iterate through rows in temp_taxonomy_df once and create a dictionary mapping id to label
    id_label_mapping = {index: row['Taxon']
                        for index, row in temp_taxonomy_df.iterrows()}
    # Grab Original seqs into a list for later use
    seqs_gen = skbio.io.read(
        reference_reads_fp, format='fasta', constructor=skbio.DNA)

    seqs_list = []
    for seq in seqs_gen:
        seqs_list.append(seq)
    # Check if input reads and taxanomy annotation len match
    if len(temp_taxonomy_df) != len(seqs_list):
        raise ValueError(
            "Mismatch in Ref-Seqs and Taxonomy ID, please use RESCRIPT in LCA or Majority mode to prep reference reads and taxonomy file!!!")
    with tempfile.TemporaryDirectory() as temp_dir:
        # Prep fasta with taxonomy annotation in labels for usearch
        with open(temp_dir + "/ref_seqs.fasta", 'a') as id_fixed:
            print("")
            print(
                "Combining Seqs ID and Taxonomy Annotation into a single fasta file for usearch...")
            print("")
            # Iterate through items in seqs_dict and write to the file
            for seq in seqs_list:
                label = id_label_mapping.get(seq.metadata['id'])
                if label is not None:
                    # Fix empty annotations in gtdb, silva and gg
                    # Split the string by colon and filter out empty elements
                    split_list = label.split(',')
                    result_list = [
                        item for item in split_list if len(item) > 2]
                    # Reconstruct the string by joining the filtered list with colon
                    label = ','.join(result_list)
                    # Since usearch really don't like '_' and '.' in seqs id, replace them with ''
                    fixed_id = str(seq.metadata['id']).replace(
                        '_', '').replace('.', '') + ";tax=" + label
                    # Fix taxa annotation to end with ';'
                    if fixed_id[-1] == ",":
                        fixed_id = fixed_id[:-1] + ";"
                    if fixed_id[-1] != ";":
                        fixed_id = fixed_id + ";"
                    seq.metadata['id'] = fixed_id
                # if label is not None:
                #     header = ">" + id + ";tax=" + label + ";"
                #     id_fixed.write(header + "\n")
                #     id_fixed.write(seqs + "\n")
                else:
                    raise ValueError(
                        "Mismatch in Ref-Seqs and Taxonomy ID, please use RESCRIPT in LCA or Majority mode to prep reference reads and taxonomy file!!!")
                # F*** your home-made solution bugged, but skbio API works, why TT
                seq.write(id_fixed)

            # # glob over the fasta file to find matching IDs
            # for i in range(len(seqs_list)):
            #     if seqs_list[i].metadata['id'] != old_id:
            #         continue
            #     else:
            #         seqs_binary = seqs_list[i].values
            #         seqs = seqs_binary.tobytes().decode()

            # # # debug
            # subprocess.run(
            #     ['cp', temp_dir + "/ref_seqs.fasta", '/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax'])

            # Build udb for usearch to work with
            build_udb_cmd = ['usearch', '-makeudb_usearch', temp_dir +
                             "/ref_seqs.fasta", '-output', temp_dir + "/ref_seqs.udb"]
            run_command(build_udb_cmd)

            # Run sintax
            sintax_cmd = ['usearch', '-sintax', query_fp, '-db', temp_dir + "/ref_seqs.udb",
                          '-strand', 'plus', '-tabbedout', temp_dir + "/reads.sintax", '-threads', str(threads)]
            run_command(sintax_cmd)
            results_fp = temp_dir + "/reads.sintax"
            # # debug
            # subprocess.run(
            #     ['cp', results_fp, '/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax'])

            classification_to_fix_order = _sintax_to_q2_classification(
                results_fp, confidence, drop_species)

    # We have a F***ing problems here
    # The sintax output is not in the same order as the input quer
    # Fix taxonomy dataframe order
    # Get the original order from the query seqs
    with open(query_fp, "rb") as seqs_file:
        seqs_gen = skbio.io.read(
            seqs_file, format='fasta', constructor=skbio.DNA)
        query_seqs_list = []
        for seq in seqs_gen:
            query_seqs_list.append(seq)
        query_seqs_id_list = []
        for i in range(len(query_seqs_list)):
            query_seqs_id_list.append(query_seqs_list[i].metadata['id'])
    # Sort the taxonomy dataframe using original id list
    # # Debug
    # print(len(query_seqs_id_list))
    # print(classification_to_fix_order)

    classification = classification_to_fix_order.loc[query_seqs_id_list]
    # # Debug
    # classification.to_csv(
    #     "/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax/debug.csv")
    return classification


def _sintax_to_q2_classification(results_fp, confidence, drop_species):
    # Read the data from the file
    sintax = pd.read_table(results_fp, sep="\t", header=None)
    sintax.rename(columns={0: 'Feature ID', 1: 'Taxon'}, inplace=True)

    # Extract the columns and split Taxon into separate columns
    taxa_working = sintax.iloc[:, 0:2].copy()

    # Drop species annotation (recommended by edgar)
    if drop_species == True:
        taxa_working[['d', 'p', 'c', 'o', 'f', 'g']
                     ] = taxa_working['Taxon'].str.split(",", expand=True)
        levels = range(2, 8)
    else:
        taxa_working[['d', 'p', 'c', 'o', 'f', 'g', 's']
                     ] = taxa_working['Taxon'].str.split(",", expand=True)
        levels = range(2, 9)

    # Fix usearch style annotation back to qiime style
    # .iloc[:, 2:8] works the same as range(2, 8), you've got the idea...
    if drop_species == True:
        taxa_working.iloc[:, 2:8] = taxa_working.iloc[:,
                                                      2:8].replace(":", "__", regex=True).replace("\+", " ", regex=True).replace("\^", "_", regex=True)
    else:
        taxa_working.iloc[:, 2:9] = taxa_working.iloc[:,
                                                      2:9].replace(":", "__", regex=True).replace("\+", " ", regex=True).replace("\^", "_", regex=True)

    # Add a new column 'Confidence' initialized with NaN
    taxa_working['Confidence'] = pd.NA

    # # debug
    # taxa_working.to_csv(
    #     "/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax/taxa_working.csv")

    # Still need to review this part, chat gpt sometimes mess up
    # Iterate through each row and process columns 2 to 8
    for i in range(len(taxa_working)):
        for j in levels:
            #            score = (re.sub(r'\D', '', str(taxa_working.iloc[i, j])))
            tmp_lst_1 = str(taxa_working.iloc[i, j]).split('(')
            try:
                tmp_lst_2 = tmp_lst_1[1].split(')')
                score = float(tmp_lst_2[0])
            except IndexError:
                score = float(0)
            if score != 0:
                #                score = float(int(score)) / 10000
                if score >= confidence:
                    taxa_working.at[i, 'Confidence'] = score
                    taxa_working.at[i, taxa_working.columns[j]] = re.sub(
                        r'\(.*\)', '', taxa_working.at[i, taxa_working.columns[j]])
                else:
                    taxa_working.at[i, taxa_working.columns[j]] = pd.NA
            else:
                taxa_working.at[i, taxa_working.columns[j]] = pd.NA

    # Special treatment to silva db
    taxa_working.iloc[:, 2:] = taxa_working.iloc[:,
                                                 2:].replace("\[", "\(", regex=True).replace("\]", "\)", regex=True)

    if drop_species == True:
        # Concatenate columns d:s to create the 'Taxon' column
        taxa_working['Taxon'] = taxa_working.loc[:, 'd':'g'].apply(
            lambda x: ';'.join(x.dropna().astype(str)), axis=1)
        # Drop the intermediate columns
        taxa_working.drop(columns=['d', 'p', 'c', 'o', 'f', 'g'], inplace=True)
    else:
        # Concatenate columns d:s to create the 'Taxon' column
        taxa_working['Taxon'] = taxa_working.loc[:, 'd':'s'].apply(
            lambda x: ';'.join(x.dropna().astype(str)), axis=1)
        # Drop the intermediate columns
        taxa_working.drop(
            columns=['d', 'p', 'c', 'o', 'f', 'g', 's'], inplace=True)

    # Set the index to Feature ID
    taxa_working.set_index('Feature ID', inplace=True)

    # Sintax won't report unclassified sequences, we will report them here
    # Replace empty strings with pd.NA in the 'Taxon' and 'Confidence' columns

    taxa_working['Taxon'] = taxa_working['Taxon'].replace('', pd.NA)

    # Set 'Unclassified' for the 'Taxon' column and 1.00 for the 'Confidence' column where both are empty
    taxa_working.loc[taxa_working['Taxon'].isna() & taxa_working['Confidence'].isna(), [
        'Taxon', 'Confidence']] = ['Unclassified', 1.00]

    # # debug
    # taxa_working.to_csv(
    #     "/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax/will_unclassified_work/sintax.csv")

    return taxa_working

    # Chat GPT META !!! This works !!!
    # print(taxa_working)

# os.chdir("/home/navi/synonas/My_Testing_Ground/q2-usearch-test/sintax/will_unclassified_work/644da2b9-af37-4ebf-8790-707cf7664db1/data")
# with open("dna-sequences.fasta", "rb") as seqs_file:
#         seqs_gen = skbio.io.read(
#             seqs_file, format='fasta', constructor=skbio.DNA)
#         query_seqs_list = []
#         for seq in seqs_gen:
#             query_seqs_list.append(seq)
#         query_seqs_id_list = []
#         for i in range(len(query_seqs_list)):
#             query_seqs_id_list.append(query_seqs_list[i].metadata['id'])

#  query_seqs_id_list


# # Testing
# from qiime2 import Artifact
# import skbio
# os.chdir("/home/navi/synonas/My_Testing_Ground/sintax")

# seqs_dict = {}
# seqs_gen = skbio.io.read("test/dna-sequences.fasta",
#                          format='fasta', constructor=skbio.DNA)
# for seqs in seqs_gen:
#     seqs_dict[seqs.metadata['id']] = seqs

# len(seqs_dict)

# print(str((seqs_dict['NR_025773.1'])))

# for key, value in seqs_dict.items():
#     with open("test.fna", "a") as seqs_file:
#         seqs_file.write(">" + key + "\n")
#         seqs_file.write(str(value) + "\n")

# seqs = Artifact.load("/mnt/amplicon-db/ncbi-targeted-loci/seqs.qza")
# seqs.export_data("test")

# # testing
# import pandas as pd
# from qiime2 import Artifact
# temp_taxonomy_df = Artifact.load("/mnt/amplicon-db/ncbi-targeted-loci/taxa.qza").view(pd.DataFrame)
# temp_taxonomy_df.loc[:, 'Taxon'] = temp_taxonomy_df.loc[:, 'Taxon'].replace(
#     "__", ":", regex=True).replace("; ", ",", regex=True)
# temp_taxonomy_df[['Taxon_tmp', 'species_to_drop']
#                 ] = temp_taxonomy_df['Taxon'].str.split("s:", expand=True)

# temp_taxonomy_df = temp_taxonomy_df.drop(
#     columns=['Taxon','species_to_drop'])

# temp_taxonomy_df.rename(columns={'Taxon_tmp': 'Taxon'}, inplace=True)
