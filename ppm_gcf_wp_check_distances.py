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

with open(sys.argv[2], 'r') as fin: # shortest_10_ppm_pngc_unique_com_update_motif_nr
  with open(sys.argv[3], 'w') as fout: # output
    data = fin.read()
    gcfs = data.split('#')
    for gcf_data in gcfs:
      lines = gcf_data.split('\n')
      gcf = lines[0][1:-4]
      wps = lines[1:-1]
      if gcf in good_gcfs:
        good_wp_lines = []
        for wp_line in wps:
          wp = wp_line.split()[0]
          if wp in good_gcfs[gcf]:
            good_wp_lines.append(wp_line)
        if good_wp_lines:
          fout.write('# ' + gcf + '.out\n')
          for line in good_wp_lines:
            fout.write(line + '\n')
