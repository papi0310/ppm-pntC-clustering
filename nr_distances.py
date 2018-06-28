import sys

good_gcfs = set()
with open(sys.argv[1], 'r') as f: # nrGenomes file
  for line in f:
    if line[0] != '#':
      good_gcfs.add(' ' + line.strip() + '.out')
with open(sys.argv[2], 'r') as fin: # distances_output file
  with open(sys.argv[3], 'w') as fout: # new file
    # for line in fin:
    #   fout.write(line)
    #   break

    groups = fin.read().split('#')
    fout.write('#' + groups[1])
    for group in groups[2:]:
      parts = group.split('\n')
      gcf = parts[0]
      if gcf in good_gcfs:
        fout.write('#' + group)
