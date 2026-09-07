#!/usr/bin/env python

# Parses votu accounting table

# Modules
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='accounting table (tsv-format)')
parser.add_argument('-a', help='downsize Algorithm option')
args = parser.parse_args()

tablefile = args.i
downsizeAlgorithm = args.a
readcounts = []
with open(tablefile,'r') as infile:
    header = infile.readline().strip().split('\t')
    columnidx = header.index('vOTU_mapped')
    for line in infile:
        row = line.strip().split('\t')
        readcount = row[columnidx]
        readcounts.append(int(readcount))

if int(downsizeAlgorithm) == 1:
    print(min(readcounts))
else:
    print("You have to select an eligible downsize Algorithm")
    sys.exit(1)
