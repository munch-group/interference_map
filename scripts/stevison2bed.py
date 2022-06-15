
from contextlib import redirect_stdout
import sys

_, map_file, bed_file = sys.argv

with open(map_file) as f:
    next(f) # skip header
    with open(bed_file, 'w') as out:
        with redirect_stdout(out):
            for line in f:
                cols = line.split()
                #cols = cols[9:] + cols[:8]
                cols = cols[1:] + cols[:1]
                if cols[0].endswith('a') or cols[0].endswith('b'): # lowercase chr2A to chr2a 
                    cols[0] = cols[0].lower()
                cols[1] = str(int(float(cols[1]) * 1000000))
                cols[2] = str(int(float(cols[2]) * 1000000))
                cols[3] = str(int(float(cols[3]) * 1000000))
                print('\t'.join(cols))
    
