
import sys, re, os
from pathlib import Path
_, nlines, prefix, suffix, file_name = sys.argv
nlines = int(nlines)
n = 0
out = None
with open(file_name) as f:	
	for i, line in enumerate(f):	
		if not i % nlines:
			if out is not None:
				out.close()
			out = open("{}.{}{}".format(prefix, n, suffix), 'w')
			n += 1
		out.write(line)
