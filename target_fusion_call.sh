#$ -S /bin/bash

# path setting
# Edit these paths before use
PYTHON=/path/to/python3
SAMTOOLS=/path/to/samtools
INTERSECT_PY=/path/to/pycode/intersect.py
FILTER_ID_PY=/path/to/pycode/filter_bam_by_id.py
MATCHED_LEN_PY=/path/to/pycode/calc_matched_len_cigar.py
GET_FIRST_CIGAR_PY=/path/to/pycode/get_first_cigar.py

# argument setting
input_bam=$1
out_dir=$2
gene_A_name=$3
gene_A_chr=$4
gene_A_start=$5
gene_A_end=$6
gene_B_name=$7
gene_B_chr=$8
gene_B_start=$9
gene_B_end=${10}
FROM=${11}
TO=${12}
gene_A_coord=${gene_A_chr}:${gene_A_start}-${gene_A_end}
gene_B_coord=${gene_B_chr}:${gene_B_start}-${gene_B_end}
fusion_name=${gene_A_name}_${gene_B_name}
echo "#######################################################################"
echo "fusion call pipeline"
echo "input_bam: $input_bam"
echo "out_dir: $out_dir"
echo "gene_A_name: ${gene_A_name}"
echo "gene_A coordinate: ${gene_A_coord}"
echo "gene_B_name: ${gene_B_name}"
echo "gene_B coordinate: ${gene_B_coord}"
echo "fusion_name: ${fusion_name}"
echo "from step ${FROM}"
echo "to step ${TO}"
echo "#######################################################################"
# filename setting
output_base_name=$(basename $input_bam .bam)
output_A=${out_dir}"/"${output_base_name}_${gene_A_name}.bam
output_B=${out_dir}"/"${output_base_name}_${gene_B_name}.bam

# output
summary=${out_dir}"/"summary.txt
breakpoint_output=${out_dir}"/"${output_base_name}_breakpoints.txt

# check existing files and softwares
if [ ! -f $summary ]; then
  touch $summary
  echo -e "filename\tfusion_name\tprimary\tsupple\tsplit_pair\tdiscordant_pair" >> $summary
fi
if [ ! -f $breakpoint_output ]; then
  touch $breakpoint_output
  echo -e "fusion_name\tread_id\tbreak_chr\tbreak_pos\tsa_break_chr\tsa_break_pos" > $breakpoint_output
fi
if [ ! -f $SAMTOOLS ]; then
  SAMTOOLS=/usr/local/bin/samtools
fi


## Preprocess
if [ $FROM -le 1 ] && [ $TO -ge 1 ]
then
echo "step 1. Preprocess bam files & find reads mapped to both target genes"
# extract target region bam
$SAMTOOLS view -h $input_bam ${gene_A_coord} > $output_A
$SAMTOOLS view -h $input_bam ${gene_B_coord} > $output_B

# extract target region ids (before MAPQ filtering; get intersection of reads that mapped to both genes)
$SAMTOOLS view $output_A | cut -f1 | sort | uniq > ${out_dir}"/"${output_base_name}_${gene_A_name}.ids
$SAMTOOLS view $output_B | cut -f1 | sort | uniq > ${out_dir}"/"${output_base_name}_${gene_B_name}.ids
cat ${out_dir}"/"${output_base_name}_${gene_A_name}.ids ${out_dir}"/"${output_base_name}_${gene_B_name}.ids | sort | uniq > ${out_dir}"/"${output_base_name}_${fusion_name}.ids

## Extract Fusion reads
$PYTHON $INTERSECT_PY \
		${out_dir}"/"${output_base_name}_${gene_A_name}.ids \
		${out_dir}"/"${output_base_name}_${gene_B_name}.ids \
		${out_dir}"/"${output_base_name}_${fusion_name}.ids

$PYTHON $FILTER_ID_PY \
		${out_dir}"/"${output_base_name}_${fusion_name}.ids \
		$input_bam \
		${out_dir}"/"${output_base_name}_${fusion_name}.bam

rm $output_A
rm $output_B
rm ${out_dir}"/"${output_base_name}_${gene_A_name}.ids
rm ${out_dir}"/"${output_base_name}_${gene_B_name}.ids
rm ${out_dir}"/"${output_base_name}_${fusion_name}.ids

echo "step 1 finished"
fi


## Sort and Index fusion bam file
if [ $FROM -le 2 ] && [ $TO -ge 2 ]
then
echo "step 2. Sort and Index fusion bam"

$SAMTOOLS sort -o ${out_dir}"/"${output_base_name}_${fusion_name}.sorted.bam ${out_dir}"/"${output_base_name}_${fusion_name}.bam
mv ${out_dir}"/"${output_base_name}_${fusion_name}.sorted.bam ${out_dir}"/"${output_base_name}_${fusion_name}.bam
$SAMTOOLS index ${out_dir}"/"${output_base_name}_${fusion_name}.bam

echo "step 2 finished"
fi


## make summary fusion reads
# filter pairs by MAPQ >= 20 in both reads
if [ $FROM -le 3 ] && [ $TO -ge 3 ]
then
echo "step 3. filter mapped read pairs by MAPQ>=20 & make summary for fusion reads"
# FLAG list
# -F  256: remove secondary alignment
# -F 2048: remove supplementary alignment
# -F 2304: remove secondary & supplementary alignment
# -F 4: remove unmapped
# -F 8: remove pair unmapped
# -F 12: both reads mapped
# -F 2316: mapped & primary alignments
# get only primary alignments; since supplementary alignments are marked within primary alignment with SA:Z tag
n_primary=$($SAMTOOLS view -F 2316 ${out_dir}"/"${output_base_name}_${fusion_name}.bam | cut -f1 | sort | uniq | wc -l)
n_supple=$($SAMTOOLS view -F 2316 ${out_dir}"/"${output_base_name}_${fusion_name}.bam | grep SA:Z | cut -f1 | sort | uniq | wc -l)
$SAMTOOLS view -F 2316 ${out_dir}"/"${output_base_name}_${fusion_name}.bam | grep SA:Z | cut -f1 | sort | uniq > ${out_dir}"/"${output_base_name}_${fusion_name}_supple.ids
echo "n_primary_alignments: $n_primary"
echo "n_supplementary_alignments: $n_supple"
# read IDs with both R1 & R2 satisfying MAPQ≥20
$SAMTOOLS view -F 2316 -q 20 -f 64 ${out_dir}"/"${output_base_name}_${fusion_name}.bam | cut -f1 | sort | uniq > ${out_dir}"/"${output_base_name}_${fusion_name}_R1_mapq20.ids
$SAMTOOLS view -F 2316 -q 20 -f 128 ${out_dir}"/"${output_base_name}_${fusion_name}.bam | cut -f1 | sort | uniq > ${out_dir}"/"${output_base_name}_${fusion_name}_R2_mapq20.ids

# 1. R1 & R2 MAPQ≥20 read IDs
$PYTHON $INTERSECT_PY \
		${out_dir}"/"${output_base_name}_${fusion_name}_R1_mapq20.ids \
		${out_dir}"/"${output_base_name}_${fusion_name}_R2_mapq20.ids \
		${out_dir}"/"${output_base_name}_${fusion_name}_R1R2_mapq20.ids

# 2. Split reads: both_mapq20 ∩ split
$PYTHON $INTERSECT_PY \
		${out_dir}"/"${output_base_name}_${fusion_name}_R1R2_mapq20.ids \
		${out_dir}"/"${output_base_name}_${fusion_name}_supple.ids \
		${out_dir}"/"${output_base_name}_${fusion_name}_split.ids

# 3. Count
n_both_mapq20=$(cat ${out_dir}"/"${output_base_name}_${fusion_name}_R1R2_mapq20.ids | wc -l | tr -d ' \n')
n_split=$(cat ${out_dir}"/"${output_base_name}_${fusion_name}_split.ids | wc -l | tr -d ' \n')
n_discordant=$((n_both_mapq20 - n_split))

single_summary="${output_base_name}	${fusion_name}	${n_primary}	${n_supple}	${n_split}	${n_discordant}"
echo $single_summary >> $summary

rm ${out_dir}"/"${output_base_name}_${fusion_name}_R1_mapq20.ids
rm ${out_dir}"/"${output_base_name}_${fusion_name}_R2_mapq20.ids
rm ${out_dir}"/"${output_base_name}_${fusion_name}_R1R2_mapq20.ids
rm ${out_dir}"/"${output_base_name}_${fusion_name}_supple.ids

fi


## spot breakpoint coordinates and organize them by genes
if [ $FROM -le 4 ] && [ $TO -ge 4 ]
then
echo "step 4. spot breakpoint coordinates"

for read_id in $(cat ${out_dir}"/"${output_base_name}_${fusion_name}_split.ids);do
  echo $read_id
  line=$($SAMTOOLS view -F 2316 -q 20 ${out_dir}"/"${output_base_name}_${fusion_name}.bam |\
	  grep $read_id | grep SA:Z | head -1) # randomly select 1 read (if two reads within the pair passes the filter)
  echo $line
  echo $line | wc -l

    if [ -n "$line" ]; then
    chr=$(echo "$line" | cut -f3)
	pos=$(echo "$line" | cut -f4)
	flag=$(echo "$line" | cut -f2)
	cigar=$(echo "$line" | cut -f6)
	sa_tag=$(echo "$line" | grep -o 'SA:Z:[^[:space:]]*' | cut -d: -f3)
	sa_chr=$(echo $sa_tag | cut -d, -f1 | cut -c2-) # wgbs: remove f/r one letter generated by bwa-meth
	#sa_chr=$(echo $sa_tag | cut -d, -f1) # wgs: bwa does not add f/r prefix to chromosome name 
	sa_pos=$(echo $sa_tag | cut -d, -f2)
	sa_strand=$(echo $sa_tag | cut -d, -f3)
	sa_cigar=$(echo $sa_tag | cut -d, -f4)

    #echo "Debug: chr=$chr, pos=$pos, flag=$flag"
    #echo "Debug: sa_chr=$sa_chr, sa_pos=$sa_pos, sa_strand=$sa_strand"

	if (( (flag & 16) == 0 )); then
	  strand="+"
	else
	  strand="-"
	fi

    # count matched length (M,D,X,N,=)
	matched_length=$($PYTHON $MATCHED_LEN_PY "$cigar")
	sa_matched_length=$($PYTHON $MATCHED_LEN_PY "$sa_cigar")
    #echo "Debug: matched_length=$matched_length, sa_matched_length=$sa_matched_length"
	# get first cigar (S:SH, M:MDXN=, ignore IP)
	first_cigar=$($PYTHON $GET_FIRST_CIGAR_PY "$cigar")
	sa_first_cigar=$($PYTHON $GET_FIRST_CIGAR_PY "$sa_cigar")
    #echo "Debug: first_cigar=$first_cigar, sa_first_cigar=$sa_first_cigar"
	# end coordinates
	end_pos=$((pos + matched_length - 1))
	sa_end_pos=$((sa_pos + sa_matched_length - 1))
    #echo "Debug: end_pos=$end_pos, sa_end_pos=$sa_end_pos"

	if [ "$strand" = "+" ]; then
	  if [ "$first_cigar" = "S" ]; then
	    breakpoint_pos=$pos
	  else
	    breakpoint_pos=$end_pos
	  fi
	else	          
	  if [ "$first_cigar" = "S" ]; then
	    breakpoint_pos=$end_pos
	  else 
	    breakpoint_pos=$pos
	  fi
	fi

	if [ "$sa_strand" = "+" ]; then
	  if [ "$sa_first_cigar" = "S" ]; then
	    sa_breakpoint_pos=$sa_pos
	  else
	    sa_breakpoint_pos=$sa_end_pos
	  fi
	else	          
	  if [ "$sa_first_cigar" = "S" ]; then
	    sa_breakpoint_pos=$sa_end_pos
	  else 
	    sa_breakpoint_pos=$sa_pos
	  fi
	fi
								        
	# Print result
	echo -e "${fusion_name}\t${read_id}\t${chr}\t${breakpoint_pos}\t${sa_chr}\t${sa_breakpoint_pos}" >> $breakpoint_output
  fi
done

echo "step 4 finished"
fi




