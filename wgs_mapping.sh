#$ -S /bin/bash

#################################################################################################
#
# WGS mapping script
#
#################################################################################################

### Configure global arguments
# Edit these paths before use
GRCH37=/path/to/Homo_sapiens.GRCh37.dna.primary_assembly.fa
PYTHON=/path/to/python # Python 3.9.5 
GATK=/path/to/gatk # GATK 4.0.5.1
FASTP=/path/to/fastp # fastp 0.20.1
SAMTOOLS=/path/to/samtools # samtools 1.10
BWA=/path/to/bwa # bwa 0.7.15-r1142-dirty
PYCODE=/path/to/pycode

#source path_to_conda.sh
#conda activate env_name

s=$1 # unique sampleID
s1=$2
s2=$3
WD=$4 # path for analysis
NGS=$5 # path to raw FASTQ
FROM=$6
TO=$7

echo 's: '$s
echo 'read1: '$s1
echo 'read2: '$s2
echo 'dir: '$WD
echo 'ngs_dir: '$NGS
echo 'from step: '$FROM
echo 'to step: '$TO

## change extension appropriately!
READ1=$NGS/${s1}.fastq.gz
READ2=$NGS/${s2}.fastq.gz

mkdir $WD/analysis

DST=$WD/analysis/${s}
if [ ! -d $DST ]
then
	mkdir $DST
	mkdir $DST/trimmed
	mkdir $DST/bam
fi
echo "directory"


if [ $FROM -le 1 ] && [ $TO -ge 1 ]
then
echo "Step 1 : Adapter trimming"
### Adapter trimming ####
### Input  : Raw FASTQ reads (paired end reads)
### Output : Trimmed FASTQ reads (paired end reads)

$FASTP \
-i ${READ1} -I ${READ2} \
-o ${DST}/trimmed/${s}_1.trimmed.fastq.gz \
-O ${DST}/trimmed/${s}_2.trimmed.fastq.gz \
-h ${DST}/trimmed/${s}.html \
-j ${DST}/trimmed/${s}.json \
-q 20 -u 20 -y -3 -p -g -t 1 -T 1 \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-w 8

fi


if [ $FROM -le 2 ] && [ $TO -ge 2 ]
then
echo "Step 2 : Alignment by BWA"
### BWA mapping ###
### Input : Trimmed FASTQ reads
### Output : Sorted BAM (.bwa.sorted.bam)

# Can use file directly or use links
TRIMMED_R1_LINK=${DST}/trimmed/${s}_1.trimmed.fastq.gz
TRIMMED_R2_LINK=${DST}/trimmed/${s}_2.trimmed.fastq.gz
TRIMMED_FASTQ_R1=$(readlink -f $TRIMMED_R1_LINK)
TRIMMED_FASTQ_R2=$(readlink -f $TRIMMED_R2_LINK)

BWA_SAM=${DST}/bam/${s}.bwa.sam
BWA_BAM=${DST}/bam/${s}.bwa.bam

$BWA mem -t 8 \
$GRCH37 \
$TRIMMED_FASTQ_R1 \
$TRIMMED_FASTQ_R2 \
> $BWA_SAM

$SAMTOOLS view -b -T $GRCH37 \
-S $BWA_SAM \
-o $BWA_BAM

rm $BWA_SAM

# Sort and Indexing
$SAMTOOLS sort -o ${DST}/bam/${s}.bwa.sorted.bam $BWA_BAM
$SAMTOOLS index ${DST}/bam/${s}.bwa.sorted.bam
rm $BWA_BAM
echo "BWA alignment completed: $BWA_BAM"
fi


if [ $FROM -le 3 ] && [ $TO -ge 3 ]
then
echo "Step 3 : Markduplicate after BWA"
### Remove duplicates ###
### Input : Sorted BAM (.bwa.sorted.bam)
### Output : Markduplicated BAM (.bwa.markdup.bam)

$GATK --java-options "-Xmx8G -Djava.io.tmpdir=./ -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
MarkDuplicates \
-I ${DST}/bam/${s}.bwa.sorted.bam \
-O ${DST}/bam/${s}.bwa.markdup.bam \
--REMOVE_DUPLICATES \
-M ${DST}/bam/${s}.bwa.sorted.metrics

$SAMTOOLS index ${DST}/bam/${s}.bwa.markdup.bam
fi



