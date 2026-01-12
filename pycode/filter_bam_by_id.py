#!/usr/bin/env python

import pysam
import sys
id_file = sys.argv[1]
input_bam = sys.argv[2]
output_bam = sys.argv[3]

# Load Read ID set
with open(id_file, 'r') as f:
    target_ids = set(line.strip() for line in f)

# Filter BAM file 
with pysam.AlignmentFile(input_bam, 'rb') as infile:
    with pysam.AlignmentFile(output_bam, 'wb', 
                           template=infile, header=infile.header) as outfile:
        for read in infile:
            if read.query_name in target_ids:
                outfile.write(read)
