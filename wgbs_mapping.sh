#$ -S /bin/bash

#################################################################################################
#
# WGBS mapping script
#
#################################################################################################

### Configure global arguments
#GRCH37=/path/to/Homo_sapiens.GRCh37.dna.primary_assembly.fa
GRCH37=/path/to/hg19_pUC19_lambda.fa # pUC19 and lambda merged hg19 genome
PYTHON=/path/to/python # Python 3.9.5
GATK=/path/to/gatk # GATK 4.0.5.1
FASTP=/path/to/fastp # fastp 0.20.1
SAMTOOLS=/path/to/samtools # samtools 1.10
METHYLDACKEL=/path/to/MethylDackel # MethylDackel 0.3.0
BWAMETH=/path/to/bwameth.py # bwa-meth 0.2.7
PYCODE=/path/to/pycode

#source /path/to/conda.sh
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
	mkdir $DST/CpG
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
echo "Step 2 : Alignment by BWA-meth"
### BWA-meth mapping ###
### Input : Trimmed FASTQ reads
### Output : Sorted BAM (.bwameth.sorted.bam)

# bwa-meth reference indexing
# this should be done once for first-time running
#PYTHON $BWAMETH index $GRCH37 

# Can use file directly or use links
TRIMMED_R1_LINK=${DST}/trimmed/${s}_1.trimmed.fastq.gz
TRIMMED_R2_LINK=${DST}/trimmed/${s}_2.trimmed.fastq.gz
TRIMMED_FASTQ_R1=$(readlink -f $TRIMMED_R1_LINK)
TRIMMED_FASTQ_R2=$(readlink -f $TRIMMED_R2_LINK)

BWMETH_BAM=${DST}/bam/${s}.bwameth.bam

# bwa-meth does not require additional reference indexing

$PYTHON $BWAMETH \
  --reference $GRCH37 \
  --do-not-penalize-chimeras \
  $TRIMMED_FASTQ_R1 $TRIMMED_FASTQ_R2 \
  > $BWMETH_BAM

# Sort and Index
$SAMTOOLS sort -o ${DST}/bam/${s}.bwameth.sorted.bam $BWMETH_BAM
$SAMTOOLS index ${DST}/bam/${s}.bwameth.sorted.bam
rm $BWMETH_BAM
echo "BWA-meth alignment completed: $BWMETH_BAM"
fi


if [ $FROM -le 3 ] && [ $TO -ge 3 ]
then
echo "Step 3 : Markduplicate after BWA-meth"
### Remove duplicates ###
### Input : Sorted BAM (.bwameth.sorted.bam)
### Output : Markduplicated BAM (.bwameth.markdup.bam)

$GATK --java-options "-Xmx8G -Djava.io.tmpdir=./ -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
MarkDuplicates \
-I ${DST}/bam/${s}.bwameth.sorted.bam \
-O ${DST}/bam/${s}.bwameth.markdup.bam \
--REMOVE_DUPLICATES \
-M ${DST}/bam/${s}.bwameth.sorted.metrics

$SAMTOOLS index ${DST}/bam/${s}.bwameth.markdup.bam
fi


if [ $FROM -le 4 ] && [ $TO -ge 4 ]
then
echo "Step 4 : Extract CpG from mark duplicated BAM"
### Extract CpG ###
### Input : Markduplicated BAM (.bwameth.markdup.bam)                
### Output : CpG context count (.markdup.mergeContext.bedGraph)
$METHYLDACKEL extract \
$GRCH37 \
${DST}/bam/${s}.bwameth.markdup.bam \
-o ${DST}/CpG/${s}_markdup

$METHYLDACKEL mergeContext \
$GRCH37 \
${DST}/CpG/${s}_markdup_CpG.bedGraph \
-o $DST/CpG/${s}.markdup.mergeContext.bedGraph

$PYTHON $PYCODE/per_sample_convert_to_beta.py \
--infile $DST/CpG/${s}.markdup.mergeContext.bedGraph \
--outfile $DST/CpG/${s}.markdup.beta.bed
fi


if [ $FROM -le 5 ] && [ $TO -ge 5 ]
then
echo "Step 5 : Bisulfite conversion efficency"
### Bisulfite conversion efficency ###
### Input : Markduplicated BAM (.bwameth.markdup.bam)                
### Output : Efficiency (.rate)

$METHYLDACKEL extract \
$GRCH37 ${DST}/bam/${s}.bwameth.markdup.bam \
--CHH -o $DST/CpG/${s}_markdup

$PYTHON /Data-node02/KTH/scripts/BS_seq_conversion_rate_v2.py \
--infile ${DST}/CpG/${s}_markdup_CHH.bedGraph \
--reference $GRCH37 \
--outfile ${DST}/CpG/${s}.conversion.rate
fi


if [ $FROM -le 6 ] && [ $TO -ge 6 ]           
then                                          
echo "Step 6 : Bisulfite conversion efficency in lambda"
### Bisulfite conversion efficency ###            
### Input : Markduplicated BAM (.bwameth.markdup.bam)                
### Output : Efficiency (.rate)                       

$METHYLDACKEL extract \
$GRCH37 ${DST}/bam/${s}.bwameth.markdup.bam \
--CHG -o $DST/CpG/${s}_markdup                       

cat ${DST}/CpG/${s}_markdup_CpG.bedGraph \
${DST}/CpG/${s}_markdup_CHG.bedGraph \
${DST}/CpG/${s}_markdup_CHH.bedGraph > ${DST}/CpG/${s}_markdup.merged
cat ${DST}/CpG/${s}_markdup.merged | grep Lambda_NEB > ${DST}/CpG/${s}_markdup_lambda.bedGraph
rm ${DST}/CpG/${s}_markdup.merged

$PYTHON /Data-node02/KTH/scripts/BS_seq_conversion_rate_v2.py \
--infile ${DST}/CpG/${s}_markdup_lambda.bedGraph \
--reference $GRCH37 \
--outfile ${DST}/CpG/${s}.lambda.conversion.rate                      
fi                                                             

