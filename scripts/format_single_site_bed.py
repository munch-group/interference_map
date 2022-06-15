from contextlib import redirect_stdout
import sys

_, orig_file, bed_file = sys.argv

with open(orig_file) as f:
    with open(bed_file, 'w') as out:
        with redirect_stdout(out):
            for line in f:
                cols = line.split()
                if not cols[0].startswith('chr'):
                	cols[0] = 'chr' + cols[0]
                print(cols[0], cols[1], int(cols[1])+1, sep='\t')

    
