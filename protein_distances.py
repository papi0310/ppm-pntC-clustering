#!/usr/bin/env python
# Some kind of program to measure distances between proteins that have hits
# based on a certain .out file (70-cutga-out output.)

from __future__ import print_function

import sys
import os
import re
import gzip
import argparse
import traceback
import itertools
import subprocess
import pandas as pd
from StringIO import StringIO
from collections import defaultdict

# Lists of the protein names that are PNGC or PEP_MUTASE
from util import PPM_MATCH_LIST, PNGC_MATCH_LIST, parse_lastcol
pieces = "Complete"

INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-' + pieces
GFF_PATH   = '/research/gmh/GENOME_DB/gff-' + pieces
LIST_PATH = 'WHEREVER THE OUTPUT TO PPM_MOTIF_FILTER IS'
OUTFILE_COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full' ]

### NOTE: CHANGE THIS IF YOU WANT THE OUTPUT COLUMNS TO BE DIFFERENTLY NAMED
### IE. PERHAPS CHANGE PNGC to LICC and PPM to "OTHER MATCH" or something -
### - MAKE SURE IT CONTAINS NO SPACES!!!
DISTANCES_COLS = ['PPM', 'PNGC', 'Distance', 'PPM_Target_Name']

# This format string a - calls str() on each arg, and then
# formats it left-aligned with 20 spaces
DATA_FMT = '{!s:<20}'


def parse_motif_matches(list_path):
    match_data = {}
    with open(list_path, 'r') as f:
        raw_lines = f.readlines()
    for line in raw_lines:
        gcf_id = line.split(':')[0]
        wp_ids = line.split(':')[1].split(',')
        match_data[gcf_id] = wp_ids
    return match_data


def main(no_motif_filter=False):
    # Print the column names
    print(''.join([DATA_FMT.format(colname) for colname in DISTANCES_COLS]))
    use_motif_filter = not no_motif_filter
    if use_motif_filter:
        motif_match_data = parse_motif_matches(LIST_PATH)

    # Reading all the .out files in gff-Complete
    out_files = os.listdir(INPUT_PATH)
    for fname in out_files:
        gcf_id = fname.split('.')[0]
        motif_matches = lambda s: s in motif_match_data.get(gcf_id, [])

        try:
            input_df = pd.read_csv(os.path.join(INPUT_PATH, fname),
                comment='#',
                header=None,
                delimiter='\s+',
                usecols=range(7),
                names=OUTFILE_COLS,
                )
            query_names = input_df[['query_name', 'target_name']]
            # We are matching the target name (i.e. PEP_mutase) with the query name (i.e. WP_05123543.01) 

            ### NOTE: TO CHANGE THE MATCH LISTS, CHANGE THEIR USES HERE, 
            ### ON LINE 83 HERE and LINE 94 BELOW

            # Now we're only creating pairs between PPM and PNGC rows.
            pep_mutase_rows = query_names.loc[
                    query_names['target_name'].isin(PPM_MATCH_LIST)]
            pep_mutase_rows = pep_mutase_rows.drop_duplicates(subset='query_name')

            # We'll use the motif_filter ONLY if the arguments specify that we should-
            # the no_motif_filter argument says NOT to use it - so we'll do it only
            # if that is false.
            if use_motif_filter:
                pep_mutase_rows = pep_mutase_rows.loc[
                        pep_mutase_rows.apply(lambda row: motif_matches(row.query_name), axis=1)]

            pngc_rows = query_names.loc[
                    query_names['target_name'].isin(PNGC_MATCH_LIST)]
            # Drop duplicates based on query_name/protein id
            pngc_rows = pngc_rows.drop_duplicates(subset='query_name')

            # Dataframe.values turns it into a list of lists, where each inner list is the row of data
            # itertools.product does the cartesian product between the two sets.
            pairs = itertools.product(pep_mutase_rows.values, pngc_rows.values)

            # Load the GFF file
            gff_fname = fname.split('.')[0] + '.gff.gz'
            with gzip.GzipFile(os.path.join(GFF_PATH, gff_fname), 'r') as gff_file:
                gff_text = gff_file.read()
            gff_lines = [line for line in gff_text.split('\n') if 'CDS' in line]
            #gff_lines = [line for line in gff_text.split('\n') if 'Protein' in line]
            gff_lines = [line for line in gff_lines if 'protein_id' in line]
            gff_text = '\n'.join(gff_lines)
            # Only interested in the lines containing gene# and the protein_id
            gff_df = pd.read_csv(StringIO(gff_text),
                comment='#',
                header=None,
                delimiter='\t',
                )

        except Exception as e:
            traceback_str = traceback.format_exc(sys.exc_info())
            commented_tb_str = '\n'.join(['# {}'.format(l) for l in traceback_str.splitlines()])
            print(commented_tb_str, file=sys.stderr)
            continue

        last_col = gff_df[gff_df.columns[-1]]
        # this function turns the last column into a dictionary from its tag=value; format
        last_col_data = last_col.apply(parse_lastcol)

        # Use apply to ensure that the rows correspond properly, and if any lambda fails, it gets 'None'
        # This is a bit more verbose, but it won't throw exceptions.
        gff_df['parent_gene'] = last_col_data.apply(lambda s: int(s['Parent'].replace('gene', '')) if 'Parent' in s else None)
        gff_df['protein_id'] = last_col_data.apply(lambda s: s.get('Name'))
        # 'Parent' gives us the gene number we care about and 'Name' is the identifier
        # Then, we find the ones where the 'name' from last_col data == 'WP_xyz'
        # And then for thsoe columns, we subtract for the corresponding pairs
        # to get the gene distance using pairs

        results = []
        for pair in pairs:
            # We select 'left' and 'right' pair members based on whether they match
            # the target and query names of the pair that we have from input_df
            try:
                ppm_gff = gff_df[gff_df.protein_id == pair[0][0]].iloc[0]
            except IndexError as e:
                print('# Error - no matching ppm protein found in gff file: {}, {}'.format(pair[0], e), file=sys.stderr)
                continue
            try:
                pngc_gff = gff_df[gff_df.protein_id == pair[1][0]].iloc[0]
            except IndexError as e:
                print('# Error - no matching pngc protein found in gff file: {}, {}'.format(pair[1], e), file=sys.stderr)
                continue

            # '0' here referring to the first column in the GFF file - the DNA/gene id
            # eg: NZ_xyzxyzxyz - if these are not from the same gene/DNA, skip it
            if ppm_gff[0] != pngc_gff[0]:
                continue

            distance = abs(ppm_gff.parent_gene - pngc_gff.parent_gene)
            results.append((pair, distance))

        print('# {}'.format(fname))
        for r in results:
            # r looks like ( (ppm_data, pngc_data), distance )
            ppm = r[0][0]
            pngc = r[0][1]
            distance = r[1]
            print(''.join(
                [DATA_FMT.format(x) for x in [ppm[0], pngc[0], distance, ppm[1]]]
            ))


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--no_motif_filter', action='store_true', help='')
    args = arg_parser.parse_args()

    main(no_motif_filter=args.no_motif_filter)
