#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, required=True, help="MethylDackel output file")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file")

args = parser.parse_args()

outfile = open(args.outfile, "w")
with open(args.infile) as infile:
	infile.readline()
	for line in infile:
		data = line.rstrip().split('\t')
		chrom, start ,end , perc, m , u = data
		b  = float(m) / float(int(m) + int(u))
		b_str = "%.4f"%(b)
		outfile.write(chrom + '\t' + start + '\t' + end + '\t' + b_str + '\n')
		
