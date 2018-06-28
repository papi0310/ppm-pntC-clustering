# coding: utf-8
import os
from util import *

OUTPUTPATH = '/Vagabundo/monica/temp/pngc+ppm/CUTGA-OUT-faa-Contig'
count = 0

def determine_empty(filename):
    global count
    empty = False
    with open(filename, 'r') as f:
        for _ in range(3):
            line = f.readline()
            if line == '' or not line.startswith('#'):
                empty = True
                break
        if not empty:
            # Read fourth line
            line = f.readline()
            if line.startswith('#'):
                empty=True
    if empty:
        # Comment this line for a dry run
        os.remove(filename)
        count += 1

def main():
    for filename in os.listdir(OUTPUTPATH.format(SETNAME)):
        filepath = os.path.join(OUTPUTPATH.format(SETNAME), filename)
        determine_empty(filepath)
    print ("Deleted {} empty files, for set {}".format(count, SETNAME))

if __name__ == '__main__':
    main()
