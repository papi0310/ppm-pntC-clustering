import sys

good_gcfs = {}
with open(sys.argv[1], 'r') as f: # ppm_sequences_gcf_complete
  for line in f:
    if line[0] == '#':
      gcf = line[2:-1]
      if gcf not in good_gcfs:
        good_gcfs[gcf] = set()
    elif line[0] == '>':
      wp = line.split()[0][1:]
      good_gcfs[gcf].add(wp)

nhoods = []
with open(sys.argv[2], 'r') as fin: # functions_output_10_update_motif_nr_all
  lines = fin.read().split('\n')
  line = lines[1]
  gcf = line.split()[0]
  wp = line.split()[5].split('+')[0]
  nhoods.append(((gcf, wp), []))
  nhoods[-1][1].append(line)
  for line in lines[2:-1]:
    gcf = line.split()[0]
    wp = line.split()[5].split('+')[0]
    if (gcf, wp) != nhoods[-1][0]:
      nhoods.append(((gcf, wp), []))
    nhoods[-1][1].append(line)

with open(sys.argv[3], 'w') as fout: # output
  fout.write('GCF\tPOS\tWP_ID\tFAMILY\tSOURCE\tMATCH_NEIGHBORHOOD\n')
  for nhood in nhoods:
    gcf = nhood[0][0]
    wp = nhood[0][1]
    if gcf in good_gcfs and wp in good_gcfs[gcf]:
      for line in nhood[1]:
        fout.write(line)
