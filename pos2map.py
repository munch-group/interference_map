import argparse
from pathlib import Path
from bisect import bisect
from collections import defaultdict


parser = argparse.ArgumentParser()

parser.add_argument("--skip-header",
                    dest="skip_header",
                    action="store_true",
                    help="Skip first line of input file")
parser.add_argument("--map-scale",
                    dest="map_scale",
                    default=1e6,
                    type=float,
                    help="Scaling of genetic map")
parser.add_argument("-rate-scale",
                    dest="rate_scale",
                    default=1,
                    type=float,
                    help="Scaling of input recombination rate")
parser.add_argument("map_file",
                    type=Path,
                    help="File name with genetic map")
parser.add_argument("positions_file",
                    type=Path,
                    help="File name with physical positions to turn into map positions")

args = parser.parse_args()

gen_map = defaultdict(list)
phys_map = defaultdict(list)
map_pos = 0.0
prev_pos, prev_rate = 0.0, 0.0
with open(str(args.map_file)) as f:
    if args.skip_header:
        next(f)
    for line in f:
        chrom, pos, rate, *_ = line.split()
        pos, rate = int(pos), float(rate)
        map_step = (pos - prev_pos) * (prev_rate / args.rate_scale) / args.map_scale
        prev_pos, prev_rate = pos, rate
        map_pos += map_step
        gen_map[chrom].append(map_pos)
        phys_map[chrom].append(pos)

for p, g in zip(phys_map['chr1'], gen_map['chr1']):
    print(p, g)

with open(str(args.positions_file)) as f:
    for line in f:
        chrom, pos, *_ =  line.split()
        pos = int(pos)
        idx = bisect(phys_map[chrom], pos)
        if 0 < idx < len(gen_map[chrom]):
            gen_left, gen_right = gen_map[chrom][idx-1], gen_map[chrom][idx]
            phys_offset = pos - phys_map[chrom][idx-1]
            rate = (gen_right - gen_left) * args.rate_scale
            gen_pos = gen_left + phys_offset * rate / args.map_scale
        else:
            gen_pos = gen_map[chrom][idx]
        print(pos, gen_pos)

