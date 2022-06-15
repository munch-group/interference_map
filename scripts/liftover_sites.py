
from contextlib import redirect_stdout
import sys
_, intervals_file_name, output_file_name = sys.argv

with open(output_file_name, 'w') as out:
	with redirect_stdout(out):
		with open(intervals_file_name) as f:
			for line in f:
				chrom, start, end, *rest = line.split()
				start, end = int(start), int(end)
				for i in range(end-start):
					print(chrom, start+i, start+i+1, '\t'.join(rest), sep='\t')
