#!/usr/bin/env python
# This script will work on a .out directory, and read the .out files
# to produce an output fasta file of all the types of WP ids matching a given
# query name, or something (ie. pngc)

# The intent of this is to produce an output fasta file, a list of these related sequences,
# for further analysis by the perl script mkBlastClusters.pl

from __future__ import print_function

import argparse
import os
import sys
import pandas as pd
import yaml

from collections import defaultdict

from util import parse_faa, write_fasta_sequences
# use the below line to add/change constants to be imported
from util import PPM_MATCH_LIST, PNGC_MATCH_LIST

# for reading outfiles
COLS = ['target_name', 't_accession', 'tlen', 'query_name' ]


def get_outfile_data(outfile_path, valid_target_names):
    # Extract just GCF_xyz from the given file path
    gcf_id = outfile_path.split('/')[-1].split('.')[0]

    with open(outfile_path, 'r') as out_file:
        # use first 4 columns, fourth column is query name, first column is target name

        input_df = pd.read_csv(outfile_path,
            comment='#',
            header=None,
            delimiter='\s+',
            usecols=range(4)
            )
        input_df.columns = COLS  # add the column names to the dataframe
        input_df = input_df[['query_name', 'target_name']]

    # Get the unique query names from input_df where target_name is in valid_target_names
    valid_rows = input_df[input_df['target_name'].isin(valid_target_names)]

    del input_df
    valid_query_names = valid_rows['query_name'].unique()

    # Now I want to build a dictionary of GCF_ to list of WP_ ids.
    wp_dict = {gcf_id: set(valid_query_names)}
    return wp_dict


def replace_extra_with_taxonomy(wp_data, gcf_id):
    # well, so, i guess i'm hoping i made them be faa_data objects
    # and i can just replace the extra_data with gcf_id mapped to taxonomy id
    # do they all share a gcf_id... no! no they don't!
    with open(TAXONOMY_MAP_PATH, 'r') as f:
        taxid_map = yaml.load(f)
    # so am i crazy, or does it just have 1 gcf for this bad boy.
    return wp_data
 

def write_fasta_output(outfile_path, fasta_path, taxonomy=False):
    wp_data = {}

    # So, first, we iterate over each .out file in the given path. (each!)
    out_files = [filename for filename in os.listdir(outfile_path) if filename.endswith('.out')]
    for filename in out_files:
        file_path = '/'.join([outfile_path, filename])

        # NOTE: if you want to CHANGE that it's looking for PNGC or some other target name
        # then change PNGC_MATCH_LIST here to some other MATCH_LIST -
        # you'll need to import it from util, above, as well
        wp_data.update(get_outfile_data(file_path, PPM_MATCH_LIST))

    for gcf_id, wp_list in wp_data.items():
        print('#',gcf_id)
        write_fasta_sequences(gcf_id, wp_list, fasta_path, taxonomy=taxonomy)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
            description=('Creates a new fasta file containing all sequences '
            'referenced in a given protein_distances.py output file.'))
    argparser.add_argument('outfile_path',
            help='Path to the directory of .out files to read as input.')
    argparser.add_argument('fasta_path',
            help='Path to directory where relevant fasta files are located.')
    argparser.add_argument('--taxonomy', action='store_true',
            help=''''Supply this to ensure that we replace fasta extra data with
            a Taxonomy number for its respective GCF''')
    args = argparser.parse_args()

    write_fasta_output(args.outfile_path, args.fasta_path, taxonomy=args.taxonomy)
