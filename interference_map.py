
import sys
import argparse
from pathlib import Path
from math import exp, log
from bisect import bisect_left, bisect
from collections import defaultdict
from contextlib import redirect_stdout

assert sys.version_info > (3, 0)

def read_tsv_file(path, skip_first=False):
    """
    Generator reading a tsv file yielding the fields
    of each line read.
    """
    with open(str(path)) as f:
        if skip_first:
            next(f)
        for line in f:
            yield line.split()


def cumprod(iterator):
    """
    Generator yielding the cumulative producuct
    of the iterable argument.
    """
    tot = None
    for x in iterator:
        if tot is None:
            tot = x
        else:
            tot *= x
        yield tot


def compute_scores(positions):
    """
    Compute interference score for each non-synonymous
    genetic map position.
    """

    def compute(positions):
        """
        Compute scores for each position in one direction
        (either upstream or downstream)
        """

        def rates(lst):
            """
            Generator yielding the interval between 
            sucessive values an iterator
            """
            for i in range(len(lst)-1):
                yield abs(lst[i] - lst[i+1]) 

        score = sum(cumprod(1 - r for r in rates(positions)))
        scores = [score]
        for rate in rates(positions):
            assert rate < 1
            score = score / (1-rate) - 1
            scores.append(score)
        return scores

    upst_scores = compute(list(reversed(positions)))
    downst_scores = compute(positions)

    scores_lookup = list()
    pos_lookup = list()
    for i, pos in enumerate(positions):
        scores_lookup.append(downst_scores[i] + upst_scores[-i-1])
        pos_lookup.append(pos)

    return pos_lookup, scores_lookup


def get_genetic_map(recombination_map_file, rate_scale=1, map_scale=1e6, skip_first=False):
    """
    Convert a recombination map to a genetic map.
    """
    gen_map = defaultdict(list)
    phys_map = defaultdict(list)
    map_pos = 0.0
    prev_pos, prev_rate = 0.0, 0.0
    for chrom, pos, rate, *_ in read_tsv_file(recombination_map_file, skip_first=skip_first):
        pos, rate = int(pos), float(rate)
        map_step = (pos - prev_pos) * (prev_rate / rate_scale) / map_scale
        prev_pos, prev_rate = pos, rate
        map_pos += map_step
        gen_map[chrom].append(map_pos)
        phys_map[chrom].append(pos)
    return phys_map, gen_map


def genetic_coord(chrom, pos, phys_map, gen_map, rate_scale=1, map_scale=1e6):
    """
    Turn a physical coordinate into a genetic coordinate.
    """
    idx = bisect(phys_map[chrom], pos)
    if 0 < idx < len(gen_map[chrom]):
        gen_left, gen_right = gen_map[chrom][idx-1], gen_map[chrom][idx]
        phys_offset = pos - phys_map[chrom][idx-1]
        rate = (gen_right - gen_left) * rate_scale
        gen_pos = gen_left + phys_offset * rate / map_scale
    else:
        gen_pos = gen_map[chrom][idx]
    return gen_pos


def get_score(pos, pos_lookup, scores_lookup):
    """
    Look up the score for a synonymous site.
    """
    idx = bisect(pos_lookup, pos)
    if idx == 0:
        # leftmost intefering pos
        right_pos, right_score = pos_lookup[idx], scores_lookup[idx]
        score = (1 - abs(pos - right_pos)) * right_score
    elif idx == len(scores_lookup):
        # rightmost interfering pos
        left_pos, left_score = pos_lookup[idx-1], scores_lookup[idx-1]
        score = (1 - abs(pos - left_pos)) * left_score
    else:
        left_pos, left_score = pos_lookup[idx-1], scores_lookup[idx-1]
        right_pos, right_score = pos_lookup[idx], scores_lookup[idx]
        if left_pos == pos:
            score = left_score
        else:
            left_rate = (1 - abs(pos - left_pos))
            right_rate = (1 - abs(pos - right_pos))
            tot = left_rate + right_rate
            score = left_score * left_rate/tot + right_score * right_rate /tot
    return score


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--output",
                        type=Path,
                        help="Output file name")
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
    parser.add_argument("--skip-first",
                        dest="skip_first",
                        action="store_true",
                        help="Skip first line of recombination map file")
    parser.add_argument("recombination_map_file",
                        type=Path,
                        help="File name with recombination map")
    parser.add_argument("selected_sites_file",
                        type=Path,
                        help="File name with physical map positions for selected sites")
    parser.add_argument("all_sites_file",
                        type=Path,
                        help="File name with physical map positions for all sites")

    args = parser.parse_args()

    if args.output:
        out_stream = open(str(args.output), 'w')
    else:
        out_stream = sys.stdout

    # get the genetic map
    phys_coords, gen_coords = get_genetic_map(args.recombination_map_file, skip_first=True)

    # for p, g in zip(phys_coords['chr1'], gen_coords['chr1']):        
    #     print(p, g)
    # print()

    # compute genetic coordinates of all nsyn sites
    nsyn_gen_map_coords = list()
    for chrom, pos, *_ in read_tsv_file(args.selected_sites_file):
        gen_pos = genetic_coord(chrom, int(pos), phys_coords, gen_coords)
        nsyn_gen_map_coords.append(gen_pos)

    #     print(chrom, pos, gen_pos)
    # print()

    # compute scores for nsyn sites
    pos_lookup, scores_lookup = compute_scores(nsyn_gen_map_coords)

    # for p, s in zip(pos_lookup, scores_lookup):
    #     print(p, s)
    # print()

    with redirect_stdout(out_stream):

        # compute scores for syn sites by interpolating between scores for nsyn sites
        for chrom, pos, *_ in read_tsv_file(args.all_sites_file):
            gen_pos = genetic_coord(chrom, int(pos), phys_coords, gen_coords)

            print(pos, get_score(gen_pos, pos_lookup, scores_lookup))

if __name__ == "__main__":
    main()

    # python interference_map.py gen.txt phys.txt all.txt
