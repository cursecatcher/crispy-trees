#!/usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import collections
import csv
import os, sys

from timeit import default_timer as timer
import wekatree



def format_time(time):
    strtime = "{:03.5f} m".format(time/60) if time >= 60 else "{:03.5f} s".format(time)
    return strtime



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--in", "-i", dest="input", action="store", type=str, required=True)
    parser.add_argument("--out", "-o", dest="output", action="store", type=str, required=True)
    parser.add_argument("--depth", "-d", dest="max_depth", action="store", type=int, default=None)
    parser.add_argument("--threshold", "-t", dest="threshold", action="store", type=float, default="1.0")
    parser.add_argument("--verbose", "-v", dest="verbose", action="store_true", default=False)

    args = parser.parse_args()
    max_depth = args.max_depth


    #obtaining files to parse
    filelist = [args.input]         #single file

    if os.path.isdir(args.input):   #entire directory
        filelist = [
            os.path.join(args.input, filename)
                for filename in os.listdir(args.input)
        ]


    start = timer()         #3, 2, 1....
    explorations = list()

    print("- Parsing and visiting trees...")

    #parsing
    for filename in filelist:
        try:
            if args.verbose:
                print("Parsing and visiting tree in {}... ".format(filename))

            tree = wekatree.DecisionTree.parse(filename, args.verbose)
            explorations.append(tree.bfs(max_depth))

        except wekatree.UnparsableTreeException:
            print("\033[01;31mCannot parse {}\033[00m".format(filename))
    else:
        #calculating max_depth if the user has not provided the proper parameter
        if max_depth is None:
            max_depth = max([bfs_result[-1][1] for bfs_result in explorations])
            print("- Max depth is {}".format(max_depth))

    end = timer()

    print("- Parsing of {} trees completed in {}".format(len(explorations), format_time(end-start)))

    #obtaining greta's stuff
    print("- Doing stuff with BFS data...", end="")

    start = timer()

    counters = [collections.Counter() for depth in range(max_depth + 1)]
    distributions = list()

    for bfs_result in explorations:
        for node, level in bfs_result:
            counters[level][node.node_id] += 1
    else:
        for level, counter in enumerate(counters):
            tot = sum(counter.values())
            distr = sorted([(k, v, v/tot) for k, v in counter.items()], key=lambda t: t[1], reverse=True)
            distributions.append((level, distr))

    end = timer()

    print("completed in {}".format(format_time(end-start)))


    print("- Output data will be written in {} folder...".format(args.output))

    try:
        os.mkdir(args.output)
    except FileExistsError:
        pass

    with open("{}/all.csv".format(args.output), "w") as fo:
        csvfo = csv.writer(fo, delimiter="\t")
        csvfo.writerow(["depth", "gene_id", "num", "freq"])

        for level, current in distributions:
            with open("{}/level_{}.csv".format(args.output, level), "w") as fo_level:
                csvlevel = csv.writer(fo_level, delimiter="\t")
                csvlevel.writerow(["gene_id", "num", "freq"])

                partial_sum = 0

                for gene, num, freq in current:
                    curr = [level, gene, num, "{:.4f}".format(freq)]

                    if partial_sum <= args.threshold:
                        csvfo.writerow(curr)
                        csvlevel.writerow(curr[1:])

                        partial_sum += freq
                    else:
                        break

    print("Done")
