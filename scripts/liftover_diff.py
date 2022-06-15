
from contextlib import redirect_stdout
import sys
# _, mapped_bed_file, unmapped_bed_file, output_bed_file = sys.argv
_, mapped_bed_file, unmapped_bed_file = sys.argv

unmapped_sites = set()

#read file with unmapped sites
with open(unmapped_bed_file) as f:
	# read them into a set for quick lookup
	for line in f:
		if not line.startswith('#'):
			chrom, start, end, *rest = line.split()
			unmapped_sites.add((chrom, start, end))
	

# open output file and redirect stdout here
# with open(output_bed_file, 'w') as out:
# 	with redirect_stdout(out):
# open bed file with mapped sites
with open(mapped_bed_file) as f:
	# read mapped sites
	for line in f:
		chrom, start, end, *rest = line.split()
		if (chrom, start, end) not in unmapped_sites:
			print(line, end='')
