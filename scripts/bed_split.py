
import sys, re, os
from pathlib import Path
_, output_dir, *output_base_names = sys.argv

chrom_regex = re.compile(r'(chr[a-zA-Z0-9]+)')
chromosomes = [chrom_regex.search(x).group(1) for x in output_base_names]

output_dir = Path(output_dir)

if not output_dir.exists():
	os.makedirs(str(output_dir))

output_files = dict()
for chrom, output_base_name in zip(chromosomes, output_base_names):
	output_path = output_dir / output_base_name
	output_files[chrom] = open(str(output_path), 'w')

for line in sys.stdin:	
	chrom = line.split()[0]
	if chrom not in output_files:
		print(line, end='', file=sys.stderr)
	else:
		output_files[chrom].write(line)
