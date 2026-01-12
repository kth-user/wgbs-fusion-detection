#!/usr/bin/env python
import sys

if len(sys.argv) < 2:
    print(0); sys.exit(0)

cigar = sys.argv[1]
matched = ['M', 'D', 'X', 'N', '=']
total = 0

import re

def parse_cigar(cigar):
    numbers = [int(x) for x in re.findall(r'\d+', cigar)]
    letters = re.findall(r'[A-Z]', cigar)
    return numbers, letters

numbers, letters = parse_cigar(cigar)
#print("Numbers:", numbers)  # [50, 10, 25, 5, 15]
#print("Letters:", letters)  # ['M', 'D', 'M', 'I', 'S']

for i in range(len(letters)):
	if letters[i] in matched:
		total += int(numbers[i])

print(total)
