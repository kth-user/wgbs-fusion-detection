#!/usr/bin/python

"""
Bs_seq_conversion_rate.py

Conversion rate of BS-seq
"""
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, required=True, help="CHH bedGraph using MethylDackel")
parser.add_argument("-r", "--reference", type=str, required=True, help="Reference FASTA")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Outpfulename ")

args = parser.parse_args()

bedgraph = open(args.infile)
bedgraph.readline()
fasta = pysam.FastaFile(args.reference)

chroms = [str(x) for x in range(1,23)]
chroms.append('X')
chroms.append('Y')
chroms.append('pUC19')
chroms.append('Lambda_NEB')
outfile = open(args.outfile, 'w')
outfile.write('chrom\tnumBases\tnumConverted\tratio\n')
chr_dict = {}
for chrom in chroms:
	chr_dict[chrom] = {'ALL' : 0, 'CONV' : 0}
for line in bedgraph:
	data = line.rstrip().split('\t')
	chrom, start, end , perc, m , u = data
	try:
		chr_dict[chrom]
	except KeyError:
		continue
	start = int(start)
	end= int(end)
	base = fasta.fetch(chrom, start, end) # single base
	if base == 'C':
		m = int(m)
		u = int(u)
		t = m + u
		chr_dict[chrom]['ALL'] += t
		chr_dict[chrom]['CONV'] += u

for chrom in chroms:
	tmp = chr_dict[chrom]
	print(chrom, "is summarizing...")
	T = str(tmp['ALL'])
	M = str(tmp['CONV'])
	try:
		r = float(M) / float(T)
	except ZeroDivisionError:
		r = 'NA'
	outfile.write(chrom + '\t' + T + '\t' + M + '\t' + str(r) + '\n')
outfile.close()
