import sys

MAX_DIST = 10

with open(sys.argv[1], 'r') as fin:
  with open(sys.argv[2], 'w') as fout:
    for line in fin:
      if line and line[0] == '#':
        comment = line
      else:
        ppm, pngc, dist = line.split()
        dist = float(dist)
        if dist <= MAX_DIST:
          fout.write(comment)
          fout.write(line)
