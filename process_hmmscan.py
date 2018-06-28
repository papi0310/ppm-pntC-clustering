#!/sw/bin/python
# DO NOT RUN - WORK IN PROGRESS

# process_hmmscan - Simultaneously process many hmmscans at once.
from __future__ import print_function
import os, multiprocessing
from subprocess import call
from util import *

DEVNULL = open('/dev/null', 'w')

def run_hmmscan(f):
    fname = f.split('.')[0]
    faafile = fname+'.faa'

    # For faa-Complete (and faa-Scaffold if we ever delete our extracted ones)
    # We need to unzip the files from the ZIP_PATH
    #faa_file = open(os.path.join(TEMP_PATH, DB_NAME, faafile), 'w')
    #unzip_code = call(['gunzip', os.path.join(ZIP_PATH, f), '-c'],
    #    stdout=faa_file)
    #faa_file.close()

    command = "/usr/local/biotools/bin/hmmscan"
    options = ['--cut_ga', '--domtblout', os.path.join(TEMP_PATH, OUTPUT_DIR, fname+'.out')]
    command_array = [command] + options + [HMM_FILE, os.path.join(TEMP_PATH, DB_NAME, faafile)]
    code = call(command_array, stdout=DEVNULL)
    return code


def process_hmmscan():
    completed_files = [f for f in os.listdir(os.path.join(TEMP_PATH, OUTPUT_DIR))]
    print(len(completed_files))

    total_files = [f for f in os.listdir(ZIP_PATH) if f.endswith('.faa.gz')]
    # total_files = [f for f in os.listdir(os.path.join(TEMP_PATH, DB_NAME)) if f.endswith('.faa')]
    print(len(total_files))

    files_to_process = [f for f in total_files if f.split('.')[0]+'.out' not in completed_files]
    print(len(files_to_process))
    print("There are {} files to process... beginning processing with parallelism of {} at a time.".format(
            len(files_to_process), MULTIPROCESSING_FACTOR))
    p = multiprocessing.Pool(MULTIPROCESSING_FACTOR)
    return_codes = p.map(run_hmmscan, files_to_process)


def main():
    process_hmmscan()

if __name__ == '__main__':
    main()
    DEVNULL.close()
