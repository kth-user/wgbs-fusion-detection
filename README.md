# WGBS-based Fusion Detection Pipeline

Pipeline for detecting gene fusions from whole-genome bisulfite sequencing data.

## Requirements
- Python 3.9.5
- GATK 4.0.5.1
- fastp 0.20.1
- samtools 1.10
- MethylDackel 0.3.0
- bwa-meth 0.2.7
- bwa 0.7.15


## Usage

### 1. Read Mapping

#### WGBS Mapping
```bash
bash wgbs_mapping.sh
```

**Arguments:**
- `sampleID`: Unique sample identifier
- `read1`: Read 1 FASTQ filename
- `read2`: Read 2 FASTQ filename
- `work_dir`: Working directory for analysis outputs
- `fastq_dir`: Path to raw FASTQ files
- `start_step`: Starting step number
- `end_step`: Ending step number

**Example:**
```bash
bash wgbs_mapping.sh sample01 sample01_R1.fq.gz sample01_R2.fq.gz \
    /path/to/workdir /path/to/fastq 1 6
```

#### WGS Mapping
```bash
bash wgs_mapping.sh
```

**Arguments:** Same as WGBS mapping

**Example:**
```bash
bash wgs_mapping.sh sample01 sample01_R1.fq.gz sample01_R2.fq.gz \
    /path/to/workdir /path/to/fastq 1 3
```

### 2. Fusion Detection
```bash
bash target_fusion_call.sh  
```

**Arguments:**
- `input_bam`: Path to mapped BAM file
- `output_dir`: Output directory for fusion call results
- `geneA_name`: Gene A symbol
- `geneA_chr`: Gene A chromosome (e.g., chr9)
- `geneA_start`: Gene A start position
- `geneA_end`: Gene A end position
- `geneB_name`: Gene B symbol
- `geneB_chr`: Gene B chromosome
- `geneB_start`: Gene B start position
- `geneB_end`: Gene B end position
- `start_step`: Starting step number
- `end_step`: Ending step number

**Example:**
```bash
bash target_fusion_call.sh sample01.bam ./fusion_output \
    BCR 22 23521891 23660224 \
    ABL1 9 133589333 133763062 \
    1 4
```

**Note:** This script supports both WGBS and WGS data. WGS-specific steps are marked with comments in the script.

## Pipeline Overview

1. **Mapping**: Align sequencing reads to reference genome
   - WGBS: BWA-meth-based bisulfite-aware alignment
   - WGS: BWA-based alignment

2. **Fusion Detection**: Identify fusion junction reads spanning two genes
   - Extract reads mapping to target gene regions
   - Identify split and discordant read pairs
   - Call fusion breakpoints


