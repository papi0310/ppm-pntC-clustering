import sys

good_gcfs = set()
with open(sys.argv[1], 'r') as f: # NRGenomes.list
  for line in f:
    good_gcfs.add(line)

def out(gcf, lines, f):
  if gcf is not None and gcf in good_gcfs:
    f.write('# ')
    f.write(gcf)
    for line in lines:
      f.write(line)

with open(sys.argv[2], 'r') as fin: # ppm_sequences_gcf_complete
  with open(sys.argv[3], 'w') as fout: # output
    gcf = None
    lines = []
    for line in fin:
      if line[0] == '#':
        out(gcf, lines, fout)
        gcf = line[2:]
        lines = []
      else:
        lines.append(line)
    out(gcf, lines, fout)

