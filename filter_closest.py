import sys

MAX_DIST = 10

with open(sys.argv[1], 'r') as fin:
  min_dists = {}
  order = {}
  order_count = 0
  for i, line in enumerate(fin):
    if line and line[0] == '#':
      comment = line
    else:
      ppm, pngc, dist = line.split()
      dist = float(dist)
      if dist > MAX_DIST:
        continue
      if (comment, ppm) not in min_dists:
        min_dists[(comment, ppm)] = (dist, line)
        order[(comment, ppm)] = order_count
        order_count += 1
      elif dist < min_dists[(comment, ppm)][0]:
        min_dists[(comment, ppm)] = (dist, line)
        order[(comment, ppm)] = order_count
        order_count += 1

lines = []
for (comment, ppm), (dist, line) in min_dists.iteritems():
  lines.append((comment, ppm, line))
# import pdb; pdb.set_trace()
lines = sorted(lines, key=lambda line: order[(line[0], line[1])])
with open(sys.argv[2], 'w') as fout:
  for line in lines:
    fout.write(line[0])
    fout.write(line[2])
