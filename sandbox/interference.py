
from itertools import islice
from math import exp, log
import argparse
from bisect import bisect_left, bisect
#from numba import jit
from pathlib import Path


def cumsum(iterator):
    tot = 0
    for x in iterator:
        tot += x
        yield tot


def cumprod(iterator):
    tot = None
    for x in iterator:
        if tot is None:
            tot = x
        else:
            tot *= x
        yield tot


def differences(pos, pos_itr):
    prev = pos
    for p in pos_itr:
        diff = abs(p - prev)
        prev = p
        yield diff


def score1(pos, pos_itr):
    # logs of distances between positions in rec prob (incl. focus pos)
    return sum(cumprod(1 - d for d in differences(pos, pos_itr)))


# same as dependent (stupid implementation, but maybe more numerically stable?)
def score2(pos, pos_itr):
    # logs of distances between positions in rec prob (incl. focus pos)
    log_probs = (log(1 - d) for d in differences(pos, pos_itr))
    # return sum of exp of each cumsum
    return sum(exp(c) for c in cumsum(log_probs))


# sanity check: should be the same as compute_stats if genetic dist between all sites is the same.
def exponential(pos, pos_itr):
    return sum(exp(-abs(pos-p)) for p in pos_itr)


def compute_scores_slow(positions, statistics, start_idx=None, end_idx=None):

    # defaults
    if start_idx is None:
        start_idx = 0
    if end_idx is None:
        end_idx = len(positions)

    length = len(positions)

    # iterator over specified records to use as focus position
    pos_iter = islice(positions, start_idx, end_idx)

    # loop focus positions
    for focus_pos in pos_iter:

        # find index dividing interfering positions at focus_pos
        idx = bisect_left(positions, focus_pos)

        # compute skip to exclude the focus site in the interfering set
        if idx < length:
            skip = int(positions[idx] == focus_pos)
        else:
            # unless it is the last site
            skip = 0

        stats = list()
        for f in statistics:
            # get iterators for downstream and upstream interfering sites
            downst_iter = reversed(positions[:idx])
            # FIGURE OUT HOW TO MAKE A REVERSED GENERATOR HERE WITHOUT 
            # SLICING THE BIG LIST...
            upst_iter = islice(positions, idx+skip, length)
            # compute scores for downstream and upstream sites
            downst_score = f(focus_pos, downst_iter)
            upst_score = f(focus_pos, upst_iter)
            stats.append(downst_score + upst_score)

        yield focus_pos, stats


def compute_scores(positions):

    def compute(positions):

        def rates(lst):
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

    scores = list()
    for i, pos in enumerate(positions):
        scores.append((pos, downst_scores[i] + upst_scores[-i-1]))
    return scores

def inject_positions_iter(scores, positions):

    pos_lookup = [pos for pos, score in scores]

    for pos in positions:
        idx = bisect(pos_lookup, pos)
        if idx == 0:
            # leftmost intefering pos
            right_pos, right_score = scores[idx]
            score = (1 - abs(pos - right_pos)) * right_score
        elif idx == len(scores):
            # rightmost interfering pos
            left_pos, left_score = scores[idx-1]
            score = (1 - abs(pos - left_pos)) * left_score
        else:
            left_pos, left_score = scores[idx-1]
            right_pos, right_score = scores[idx]
            if left_pos == pos:
                score = left_score
            else:
                left_rate = (1 - abs(pos - left_pos))
                right_rate = (1 - abs(pos - right_pos))
                tot = left_rate + right_rate
                score = left_score * left_rate/tot + right_score * right_rate /tot

        yield (pos, score)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--start",
                        dest="start",
                        type=int,
                        default=None,
                        help="Start index")
    parser.add_argument("-e", "--end",
                        dest="end",
                        type=int,
                        default=None,
                        help="End index")
    parser.add_argument("--neutral",
                        dest="neutral_sites",
                        type=Path,
                        help="File name with genetic map positions for neutral sites")
    parser.add_argument("selected_sites",
                        type=Path,
                        help="File name with genetic map positions for selected sites")

    args = parser.parse_args()

    # all_positions = list()
    # interfering_positions = list()
    # with open(str(args.selected_file)) as f:
    #     for line in f:
    #         chrom, pos, kind = line.split()
    #         all_positions.append(pos)
    #         if kind == 'nsyn':
    #             interfering_positions.append(pos)

    # DUMMY: #####################
    import numpy
    all_genetic_map_coord = list(numpy.linspace(0, 1, num=20000000))
#    all_genetic_map_coord = list(sorted(numpy.abs(numpy.random.normal(0.01, 0.01, size=20))))
#    all_genetic_map_coord = [0.01 * x**2 for x in range(20)]
    interfering_genetic_map_coord = all_genetic_map_coord[0:len(all_genetic_map_coord):10]
    ##############################

    if args.start is None: args.start = 0
    if args.end is None: args.end = len(all_genetic_map_coord)

    # statistics = [score1]
    # scores1 = compute_scores_slow(interfering_genetic_map_coord, statistics,
    #                               start_idx=args.start, end_idx=args.end)
    # for pos, score in scores1:
    #     print(*([pos] + score), sep='\t')
    # print()

    scores2 = compute_scores(interfering_genetic_map_coord)
    for pos, score in scores2:
        print(pos, score, sep='\t')
    print()

    for pos, score in inject_positions_iter(scores2, all_genetic_map_coord):
        print(pos, score, sep='\t')


if __name__ == "__main__":
    main()
