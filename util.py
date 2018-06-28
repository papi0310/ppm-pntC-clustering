from __future__ import print_function

from subprocess import call
from collections import defaultdict
import os
import sys
import pandas as pd
import gzip
import yaml

# Constants for now - From process_hmmscan
SETNAME = 'Contig'
DB_NAME = 'faa-{}'.format(SETNAME)
TEMP_PATH ='/Vagabundo/monica/temp/'
#OUTPUT_DIR ='CUTGA-OUT-{}'.format(DB_NAME)
OUTPUT_DIR = 'pngc+ppm/CUTGA-OUT-{}'.format(DB_NAME)
ZIP_PATH = '/research/gmh/GENOME_DB/{}'.format(DB_NAME)

TAXONOMY_MAP_PATH = 'CompleteTaxid.yml'

# From/for sort_matches.py and protein_distances.py
PPM_MATCH_LIST =['PEP_mutase']
PNGC_MATCH_LIST =['NTP_transferase', 'NTP_transf_3', 'IspD']
LICA_MATCH_LIST = ['Choline_kinase']
GLUTAMINE_LIST = ['PEP-utilizers', 'PPDK_N']
LICB_MATCH_LIST = ['EamA']
SWAG_MATCH_LIST = ['PEP_mutase', 'Choline_kinase', 'PEP-utilizers', 'PPDK_N', 'EamA']

HMM_FILE = '/Vagabundo/monica/Proteins/models.hmm'
#HMM_FILE = '/Vagabundo/monica/Proteins/cluster_models'
PFAM_DB = '/Vagabundo/monica/Proteins/Pfam-A.hmm'

MULTIPROCESSING_FACTOR = 100 # We'll run 100 per batch


class ParseError(Exception):
    pass


class FastaData(object):
    def __init__(self, wp_id, sequence_data, extra_data=''):
        self.wp = wp_id
        self.sequence = sequence_data
        self.extra_data = extra_data


class MatchData(object):
    def __init__(self, left_wp_id, right_wp_id, distance, extra_data=''):
        self.left_wp = left_wp_id
        self.right_wp = right_wp_id
        self.distance = distance
        self.extra_data = extra_data


# You could pass in a lambda/function to cond and have that also filter out
# false results from idir (ie. is file in idir also in sort_matches both.txt)
def files_remaining(idir, odir, cond=lambda x: true):
    """Pass None if you want no filtering to happen...
    Although the default value is also a pass-thru
    """
    files = [f for f in os.listdir(idir) if f not in os.listdir(odir)]
    if cond is not None:
        files = [f for f in files if cond(f)]
    return files


def parse_lastcol(col_text):
    """ Breaking the last column of GFF file into segments eg.
    Parent=gene14
    data = {key:value for a=b;y=z} in col_text
    """
    data = col_text.split(';')
    for seg in data:
        if seg.count('=') > 1:
            raise ParseError('Malformed lastcol string: {}'.format(col_text))
    data = [segment.split('=') for segment in data if segment]
    data = {segment[0]:segment[1] for segment in data}
    return data


BLASTDB = '/research/gmh/GENOME_DB/blastpDB-Complete/GCF_{}'
def blast_motif_match(qname):
    command = 'blastdbcmd -db ' + BLASTDB + ' -target_only -entry ' + qname + ' -outfmt %s'
    output = subprocess.check_output(command, shell = True)
    has_motif = re.search('EDK\w{3,7}NS',output)
    return bool(has_motif)


def parse_faa(faa_path):
    if faa_path.endswith('.gz'):
        with gzip.open(faa_path, 'r') as f:
            file_contents = f.read()
    else:
        with open(faa_path, 'r') as f:
            file_contents = f.read()
    raw_sequence_list = file_contents.split('>')

    del file_contents # Save memory! Helps the computer run better!

    # Now, I need to turn each item into the form WP_id : sequence
    faa_data = {}

    # We know that the WP_id is the first space delimited part of the first line.
    for raw_seq in raw_sequence_list:
        if not raw_seq: # we know for sure the first result will be empty
            continue
        seq_lines = raw_seq.split('\n')
        first_line_data = seq_lines[0].split()
        wp_id = first_line_data[0]
        extra_data = ' '.join(seq_lines[0].split()[1:]) if len(first_line_data) > 1 else ''
        sequence = ''.join(seq_lines[1:])
        fd = FastaData(wp_id, sequence, extra_data=extra_data)
        faa_data[wp_id] = fd

    return faa_data


def write_fasta_sequence(fasta_data, line_length=80, output_file=sys.stdout):
    print('>{} {}'.format(fasta_data.wp, fasta_data.extra_data), file=output_file)
    for i in range(0, len(fasta_data.sequence), line_length):
        print(fasta_data.sequence[i:i+line_length], file=output_file)


def write_fasta_sequences(gcf_id, wp_ids, fasta_path, taxonomy=False):
    try:
        faa_data = parse_faa('/'.join([fasta_path, gcf_id+'.faa.gz']))
    except IOError as e:
        try:
            faa_data = parse_faa('/'.join([fasta_path, gcf_id+'.faa']))
        except IOError as e:
            print('# IOError: (write_fasta_sequences) {}: gcf_id {}'.format(e, gcf_id))
            return

    if taxonomy:
        with open(TAXONOMY_MAP_PATH, 'r') as tm_path:
            taxid_map = yaml.load(tm_path)
        for k, v in faa_data.items():
            v.extra_data = taxid_map.get(gcf_id, 'NOT_FOUND')

    for wp_id in wp_ids:
        try:
            fasta_data = faa_data[wp_id]
        except KeyError as e:
            print('# KeyError: (write_fasta_sequences) {}: gcf_id {}'.format(e, gcf_id))
            return

        write_fasta_sequence(fasta_data)


def run_hmmscan(input_path, output_path, model_db_file=PFAM_DB):
    """
    Run hmmscan on a given input (.faa) file, get output (.out) at output_path
    @input_path: must point to a valid fasta file (UNZIPPED/RAW)
    @output_path: path where we will put --domtblout output
    """

    DEVNULL = open('/dev/null', 'w')
    command = "/usr/local/biotools/bin/hmmscan"
    options = ['--cut_ga', '--domtblout', output_path]
    command_array = [command] + options + [model_db_file, input_path]
    code = call(command_array, stdout=DEVNULL)
    DEVNULL.close()
    # TODO/NOTE: should probably raise exception if code != 0
    return code


def get_family(out_path):
    """
    Look up the out_path file as a .out file, and retrieve the 'family' name of
    the first row (highest e value) (NOTE: ?)
    """
    # use first 2 columns, (for no apparent reason)
    # fourth column is query name, first column is target name
    try:
        input_df = pd.read_csv(out_path,
            comment='#',
            header=None,
            delimiter='\s+',
            usecols=range(2)
            )
        first_row = input_df.loc[0].values
    except pd.io.common.EmptyDataError:
        return 'HYPOTHETICAL'

    target_name = first_row[0]  # Column 0 is the Target Name
    return target_name


def parse_distances_file(distances_path):
    with open(distances_path, 'r') as distance_file:
        distances_raw = distance_file.readlines()

    if not distances_raw[0].startswith('# '): # skip column headers, if they exist
        distances_raw = distances_raw[1:]

    # Now I want to build a dictionary of GCF_ to list of matches
    match_data = defaultdict(set)
    gcf_id = None

    for line in distances_raw:
        if line.startswith('#'): # comment
            try:
                line_data = line.split()[1:]
            except Exception as e:
                print("Failed to parse comment line, empty comment. Error: {}".format(e),
                        file=sys.stderr)
                continue
            if line_data[0].startswith('GCF_'):
                gcf_id = line_data[0].split('.')[0] # remove the .out from the end

        else: # it's a data line
            if gcf_id is None:
                raise Exception("Error: Can't start adding data before finding a valid GCF id")
            line_data = line.split()
            match = MatchData(line_data[0], line_data[1], line_data[2], line_data[3:])
            match_data[gcf_id].add(match)
    return match_data
