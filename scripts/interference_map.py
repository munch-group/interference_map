
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


def compute_scores(positions, generations=1):
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
                # prop of recomb in one generation
                rate = abs(lst[i] - lst[i+1])

                # rate = 1 - exp(-rate) # transform to 0-1 range...
                if rate > 1:
                    rate = 1

                assert rate <= 1, rate
                # yield prop of at least one recomb in a number of generations
                yield 1 - (1 - rate)**generations

        scores = [0]
        for rate in rates(positions):
            assert 0 <= rate <= 1, rate
            score = (1-rate) + (1-rate)*scores[-1]
            assert score >= 0
            scores.append(score)
        return scores

    upst_scores = compute(positions)
    downst_scores = compute(list(reversed(positions)))

    # def compute(positions):
    #     """
    #     Compute scores for each position in one direction
    #     (either upstream or downstream)
    #     """

    #     def rates(lst):
    #         """
    #         Generator yielding the interval between 
    #         sucessive values an iterator
    #         """
    #         for i in range(len(lst)-1):
    #             # prop of recomb in one generation
    #             rate = abs(lst[i] - lst[i+1])
    #             print(rate)
    #             if rate > 0.5:
    #                 rate = 0.5
    #             assert rate < 1, rate
    #             # yield prop of at least one recomb in a number of generations
    #             yield 1 - (1 - rate)**generations

    #     score = sum(cumprod(1 - r for r in rates(positions)))
    #     scores = [score]
    #     for rate in rates(positions):
    #         score = score / (1-rate) - 1
    #         scores.append(score)
    #     return scores

    # upst_scores = compute(list(reversed(positions)))
    # downst_scores = compute(positions)


    scores_lookup = list()
    pos_lookup = list()
    for i, pos in enumerate(positions):
        scores_lookup.append(downst_scores[i] + upst_scores[-i-1])
        pos_lookup.append(pos)

    return pos_lookup, scores_lookup

# def get_genetic_map(recombination_map_file, rate_scale=1, map_scale=1e6, col_idx=(0,1,2), skip_first=False):
#     """
#     Convert a recombination map to a genetic map.
#     """
#     gen_map = defaultdict(list)
#     phys_map = defaultdict(list)
#     map_pos = 0.0
#     prev_pos, prev_rate = 0.0, 0.0

#     for row in read_tsv_file(recombination_map_file, skip_first=skip_first):
#         chrom, start, end, rate = [row[i] for i in col_idx]

#         try:
#             start, rate = int(start), float(rate)
#         except ValueError:
#             continue

#         # interval in megabases ((end-start)/map_scale) times rate in cM/Mb (prev_rate)
#         # (rate can be divided by a constant say Ne if it is a rho value). The whole thing is 
#         # divided by a megabase (map_scale) to get genetic distance as probability between two markers.
#         #map_step = ((start - prev_start) / map_scale) * (prev_rate / rate_scale) / map_scale

#         # span between physical positions on bases times rate (or rho) per base
#         map_step = (start - prev_start) * (prev_rate / rate_scale) / map_scale
        
#         prev_start, prev_rate = start, rate
#         assert map_step >= 0
#         assert start >= prev_start
#         map_pos += map_step
#         gen_map[chrom].append(map_pos)
#         phys_map[chrom].append(start)

#     # assert(sorted(gen_map[chrom]) == gen_map[chrom])
#     # assert(sorted(phys_map[chrom]) == phys_map[chrom])

#     return phys_map, gen_map

def get_genetic_map(recombination_map_file, rate_scale, map_scale, dummy_rate,
                    col_idx=(0,1,2,3), skip_first=False):
    """
    Convert a recombination map to a genetic map.
    """
    gen_map = defaultdict(list)
    phys_map = defaultdict(list)
    map_pos = 0.0
    prev_start, prev_end, prev_rate = None, None, None

    for row in read_tsv_file(recombination_map_file, skip_first=skip_first):
        chrom, start, end, rate = [row[i] for i in col_idx]

        try:
            start, end, rate = int(start), int(end), float(rate)
        except ValueError:
            # skip rows with na
            continue

        assert prev_start is None or start >= prev_start, (start, prev_start)
        assert prev_end is None or start >= prev_end, (start, prev_end, start >= prev_end)

        if prev_end is not None and start > prev_end:
            # there is a gap in the map so we need to compute a dummy map step
            map_step = (start - prev_end) * (dummy_rate / rate_scale) / map_scale
            assert map_step >= 0
            map_pos += map_step
            gen_map[chrom].append(map_pos)
            phys_map[chrom].append(start) # start not end here

        # compute map step from recombination rate in interval
        map_step = (end - start) * (rate / rate_scale) / map_scale        
        assert map_step >= 0
        map_pos += map_step
        gen_map[chrom].append(map_pos)
        phys_map[chrom].append(end)

        prev_start, prev_end, prev_rate = start, end, rate


    assert(sorted(gen_map[chrom]) == gen_map[chrom])
    assert(sorted(phys_map[chrom]) == phys_map[chrom])

    return phys_map, gen_map

# def get_genetic_map(recombination_map_file, rate_scale=1, map_scale=1e6, col_idx=(0,1,2), skip_first=False):
#     """
#     Convert a recombination map to a genetic map.
#     """
#     gen_map = defaultdict(list)
#     phys_map = defaultdict(list)
#     map_pos = 0.0
#     prev_pos, prev_rate = 0.0, 0.0

#     for row in read_tsv_file(recombination_map_file, skip_first=skip_first):
#         chrom, pos, end, rate = [row[i] for i in col_idx]

#         try:
#             pos, rate = int(pos), float(rate)
#         except ValueError:
#             continue

#         # interval in megabases ((end-start)/map_scale) times rate in cM/Mb (prev_rate)
#         # (rate can be divided by a constant say Ne if it is a rho value). The whole thing is 
#         # divided by a megabase (map_scale) to get genetic distance as probability between two markers.
#         #map_step = ((pos - prev_pos) / map_scale) * (prev_rate / rate_scale) / map_scale

#         # span between physical positions on bases times rate (or rho) per base
#         map_step = (pos - prev_pos) * (prev_rate / rate_scale) / map_scale
        
#         prev_pos, prev_rate = pos, rate
#         assert map_step >= 0
#         assert pos >= prev_pos
#         map_pos += map_step
#         gen_map[chrom].append(map_pos)
#         phys_map[chrom].append(pos)

#     # assert(sorted(gen_map[chrom]) == gen_map[chrom])
#     # assert(sorted(phys_map[chrom]) == phys_map[chrom])

#     return phys_map, gen_map


#def genetic_coord(chrom, pos, phys_map, gen_map, rate_scale=1, map_scale=1e6):
def genetic_coord(chrom, phys_pos, phys_map, gen_map):
    """
    Turn a physical coordinate into a genetic coordinate.
    """
    idx = bisect(phys_map[chrom], phys_pos)
    if 0 < idx < len(phys_map[chrom]):
        phys_left, phys_right = phys_map[chrom][idx-1], phys_map[chrom][idx]
        gen_left, gen_right = gen_map[chrom][idx-1], gen_map[chrom][idx]

#        phys_offset = phys_pos - phys_map[chrom][idx-1]
        # rate = (gen_right - gen_left) * rate_scale
        # gen_pos = gen_left + phys_offset * rate / map_scale

        per_base_interval_rate = (gen_right - gen_left) / (phys_right - phys_left)
        assert per_base_interval_rate >= 0, per_base_interval_rate
        assert phys_pos - phys_left >= 0, (phys_pos, phys_left)
        gen_pos = gen_left + (phys_pos - phys_left) * per_base_interval_rate

#        rate = (gen_right - gen_left)
#        gen_pos = gen_left + phys_offset * rate
    else:
        #gen_pos = gen_map[chrom][idx]
        raise ValueError('position outside map: {} {}'.format(chrom, phys_pos))
    return gen_pos


def get_score(pos, pos_lookup, scores_lookup):
    """
    Look up the score for a synonymous site.
    """
    idx = bisect(pos_lookup, pos)
    if idx == 0:
        # leftmost interfering pos
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
    parser.add_argument("--rate-scale",
                        dest="rate_scale",
                        default=1,
                        type=float,
                        help="Scaling of input recombination rate")
    parser.add_argument("--generations",
                        dest="generations",
                        default=1,
                        type=int,
                        help="Number of generations relevant to interference (Ne is a good number).")
    parser.add_argument("--map-col-idx",
                        dest="map_col_idx",
                        default='0,1,2,3',
                        type=str,
                        help="chrom, start, end, rate indexes in map file")
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
    map_col_idx = [int(x) for x in args.map_col_idx.split(',')]
    # phys_coords, gen_coords = get_genetic_map(args.recombination_map_file, col_idx=map_col_idx,
    #     rate_scale=args.rate_scale, map_scale=args.map_scale, skip_first=args.skip_first)
    phys_coords, gen_coords = get_genetic_map(args.recombination_map_file,
        rate_scale=args.rate_scale, map_scale=args.map_scale, dummy_rate=1.0,
        col_idx=map_col_idx, skip_first=args.skip_first)

    with open('phys_gen.txt', 'w') as f:
        for p, g in zip(phys_coords['chr10'], gen_coords['chr10']):
            print(p, g, sep='\t', file=f)

    # print(len(phys_coords['chr22']))

    # print(phys_coords['chr22'][:10])
    # print(gen_coords['chr22'][:10])
    # assert 0

    # with open(str(args.recombination_map_file)) as f:
    #     print(f.read())
    # for p, g in zip(phys_coords['chr1'], gen_coords['chr1']):        
    #     print(p, g, g*args.map_scale)
    # print()


    # import pandas as pd
    # import GenomicIntervals
    # df = pd.read_table(args.selected_sites_file, header=None, names=['chrom', 'start', 'end', 'score'])
    # df_collapsed = GenomicIntervals.interval_collapse(df)
    # raw = (df.end - df.start).sum()
    # collapsed = (df_collapsed.end - df_collapsed.start).sum()
    # print("###", raw, collapsed)
    # assert 0

    # compute genetic coordinates of all nsyn sites
    nsyn_gen_map_coords = list()
#    for chrom, pos, *_ in read_tsv_file(args.selected_sites_file):
    for chrom, start, end, *_ in read_tsv_file(args.selected_sites_file):
        for phys_pos in range(int(start), int(end)):
            try:
                # gen_pos = genetic_coord(chrom, int(pos), phys_coords, gen_coords, 
                #     rate_scale=args.rate_scale, map_scale=args.map_scale)
                gen_pos = genetic_coord(chrom, int(phys_pos), phys_coords, gen_coords)

                nsyn_gen_map_coords.append(gen_pos)
            except ValueError as err:
                print(err, file=sys.stderr)

    with open('nsyn.txt', 'w') as f:
        for n in nsyn_gen_map_coords:
            print(n, file=f)


    print(len(nsyn_gen_map_coords))

    # compute scores for nsyn sites
    pos_lookup, scores_lookup = compute_scores(nsyn_gen_map_coords, generations=args.generations)


    with open('gen_scores.txt', 'w') as f:
        for p, s in zip(pos_lookup, scores_lookup):
            print(p, s, sep='\t', file=f)


    with redirect_stdout(out_stream):

        # compute scores for syn sites by interpolating between scores for nsyn sites
        for chrom, pos, *_ in read_tsv_file(args.all_sites_file):
            try:
                # gen_pos = genetic_coord(chrom, int(pos), phys_coords, gen_coords, 
                #     rate_scale=args.rate_scale, map_scale=args.map_scale)
                gen_pos = genetic_coord(chrom, int(pos), phys_coords, gen_coords)

                print(chrom, pos, get_score(gen_pos, pos_lookup, scores_lookup))
            except ValueError as err:
                print(err, file=sys.stderr)

if __name__ == "__main__":
    main()

    ##  python ./scripts/interference_map.py --generations 10000 --map-col-idx 0,1,4 --map-scale 1000 --rate-scale 100000  --output test.scores  /faststorage/project/hri/people/kmt/steps/interference_score/chimp/map_species_coord/sorted/chr22.bed /faststorage/project/hri/people/kmt/steps/interference_score/chimp/gerp/annotation_species_coord/sorted/chr22.bed /faststorage/project/hri/people/kmt/steps/interference_score/chimp/sites_species_coord/sorted/chr22.bed
    #python ./scripts/interference_map.py --map-col-idx 0,1,2,4 --map-scale 1000 --rate-scale 1  --output test.scores  /faststorage/project/hri/people/kmt/steps/interference_score/chimp/converted_map_files/chr10.bed /faststorage/project/hri/people/kmt/steps/interference_score/chimp/gerp/annotation_species_coord/sorted/chr10.bed /faststorage/project/hri/people/kmt/steps/interference_score/chimp/sites_species_coord/sorted/chr10.bed
