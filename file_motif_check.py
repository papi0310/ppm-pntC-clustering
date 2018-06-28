import sys
import re

MOTIF_REGEX = 'EDK\w{3,7}NS'
#MOTIF_REGEX = 'Atopobium'

def has_motif(content_string, motif_regex):
  return bool(re.search(motif_regex, content_string))

with open(sys.argv[1], 'r') as fin:
  with open(sys.argv[2], 'w') as fout:
    seqs = fin.read().split('#')
    for seq in seqs:
      if has_motif(seq, MOTIF_REGEX):
        fout.write('#')
        fout.write(seq)
