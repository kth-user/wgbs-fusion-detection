#!/usr/bin/env python
# find intersect between given text (line-by-line)
import sys

if len(sys.argv) != 4:
   print("Usage: python script.py <file1> <file2>")
   sys.exit(1)

infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile = sys.argv[3]

try:
   # Read lines from both files
   with open(infile1, 'r') as f1:
       lines1 = set(line.strip() for line in f1)
   
   with open(infile2, 'r') as f2:
       lines2 = set(line.strip() for line in f2)
   
   # Find intersection
   intersection = lines1 & lines2
   
   # Print results
   with open(outfile, 'w') as f_out:
       for line in sorted(intersection):
           f_out.write(line + '\n')
       
except FileNotFoundError as e:
   print(f"Error: {e}")
   sys.exit(1)
except Exception as e:
   print(f"Error: {e}")
   sys.exit(1)



