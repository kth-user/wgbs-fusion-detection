#!/usr/bin/env python

import sys

text=sys.argv[1]
cigar_list=[]

for n in text:
	if n in ['0','1','2','3','4','5','6','7','8','9']:
		continue
	if n in ['I','P']:
		continue
	cigar_list.append(n)

first_cigar = cigar_list[0]
if first_cigar == 'S' or first_cigar == 'H':
	answer='S'
else:
	answer='M'
print(answer)
