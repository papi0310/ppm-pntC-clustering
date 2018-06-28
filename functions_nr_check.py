import sys

good_gcfs = set()
with open(sys.argv[1], 'r') as f: # nrGenomes file
  for line in f:
    if line[0] != '#':
      good_gcfs.add(line.strip())
good_gcfs.add('GCF')
with open(sys.argv[2], 'r') as fin: # functions_output file
  with open(sys.argv[3], 'w') as fout: # new file
    for line in fin:
      gcf = line.split()[0]
      if gcf in good_gcfs:
        fout.write(line)
