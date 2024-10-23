# ucircfull: a UMI-guided tool to identify and quantify circRNAs from ucircFL-seq data

## Introduction

The ucircfull package is a UMI-guided tool for identifying and quantifying circRNAs from UMI-tagged full-length circRNA sequencing (ucircFL-seq) data. It is designed to handle the challenges of low read depth and high sequencing error rates in ucircFL-seq data, and to provide accurate quantification of circRNA abundance.

## Installation

### Install ucircfull from source code:

#### Dependencies:

- cmake >= 3.4
- make
- g++ >= 11.0
- gcc >= 11.0
- libboost-all-dev
- libseqan3-dev >= 3.4.0
- libseqan2-dev
- minimap2
- rust
- porechop
- seqkit
- samtools

#### Build ucircfull from source code:

```bash
git clone https://github.com/yangence/ucircfull.git
cd ucircfull
mkdir build & cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

### Install ucircfull from apptainer:

```bash
wget https://github.com/yangence/ucircfull/releases/download/v1.0.0/ucircfull-1.0.0.sif
```

## Required files

Users can prepare the external files under the following instructions:

1. Indexed genome fasta file

```bash
samtools faidx $genome
```

2. gene annotation GTF file

## Usage

```bash
Usage: ucircfull [--help] [--version] {circ_call,clust_umi,extract_umi}

Optional arguments:
  -h, --help    shows help message and exits
  -v, --version prints version information and exits

Subcommands:
  circ_call     circRNA identification.
  clust_umi     UMI clusting guided consensus generation.
  extract_umi   Extract UMI sequence and identify strand from ucircFL-seq raw fastq
```

### Step 1: UMI extraction

```bash
Usage: ucircfull extract_umi [--help] [--version] --input FQ --anchorx SEQ --umi SEQ [--noumi] --outdir DIR --prefix PREFIX --thread INT [--seqkit PATH] [--porechop PATH]

Extract UMI sequence and identify strand from circFL-seq2 raw fastq

Optional arguments:
  -h, --help            shows help message and exits
  -v, --version         prints version information and exits
  -i, --input FQ        ucircFL-seq raw fastq file. [required]
  -x, --anchorx SEQ     anchor sequence used in 1st strand cDNA synthesis. [required]
  -u, --umi SEQ         umi pattern. [default: "CTCNNNYRNNNYRNNNYRNNNGAG"]
  -n, --noumi           no UMIs were added to 1st strand anchor.
  -o, --outdir DIR      output directory. [default: "."]
  -p, --prefix PREFIX   output prefix. [default: "circFL"]
  -t, --thread INT      number of threads used. [default: 4]
  --seqkit PATH         path to seqkit. [default: "seqkit"]
  --porechop PATH       path to porechop. [default: "porechop"]
```

For ucircFL-seq data in default library preparation data (`$rawfastq`), run:

```bash
ucircfull extract_umi -i $rawfastq -x CTACACGACGCTCTTCCGATCT -o . -p $sample -t $thread
```

#### output:

- `$sample`\_strand.fastq
- `$sample`\_umi.fasta

### Step 2: UMI clustering

```bash
Usage: ucircfull clust_umi [--help] [--version] --input FQ --umi FA [--seq VAR] --outdir DIR --prefix PREFIX --thread INT [--seqkit PATH]

UMI clusting guided consensus generation.

Optional arguments:
  -h, --help            shows help message and exits
  -v, --version         prints version information and exits
  -i, --input FQ        stranded fastq file. [required]
  -u, --umi FA          umi fasta file. [required]
  -s, --seq             generate sequence cluster. [default: false]
  -o, --outdir DIR      output directory. [default: "."]
  -p, --prefix PREFIX   output prefix. [default: "circFL"]
  -t, --thread INT      number of threads used. [default: 4]
  --seqkit PATH         path to seqkit. [default: "seqkit"]
```

```bash
ucircfull clust_umi -i ${sample}_strand.fastq -u ${sample}_umi.fasta -o . -p $sample -t $thread
```

#### output:

- `$sample`\_umi.clstr: UMI clustering result
- `$sample`\_seq.clstr (optional): sequence clustering result

### Step3: circRNA identification

```bash
Usage: ucircfull circ_call [--help] [--version] --mode STR [--notstranded VAR] --input FQ [--bam BAM] --ref REF --anno GTF [--umi CLSTR] --splice MOTIF --outdir DIR --prefix PREFIX --thread INT [--minimap2 PATH] [--samtools PATH] [--debug]

circRNA identification.

Optional arguments:
  -h, --help            shows help message and exits
  -v, --version         prints version information and exits
  -m, --mode STR        circRNA calling mode. (RG, ucRG) [required]
  -i, --input FQ        stranded fastq file. [required]
  --bam BAM             mapped bam file as input.
  -r, --ref REF         CIRI-long reference directory. [required]
  -a, --anno GTF        reference annotation GTF file. [required]
  -u, --umi CLSTR       umi clust results file. [default: "-"]
  --splice MOTIF        output directory. [default: "AGGT,AGGC,ACAT,ACGT,AGAT"]
  -o, --outdir DIR      output directory. [default: "."]
  -p, --prefix PREFIX   output prefix. [default: "circFL"]
  -t, --thread INT      number of threads used. [default: 4]
  --minimap2 PATH       path to minimap2. [default: "minimap2"]
  --samtools PATH       path to samtools. [default: "samtools"]
  --debug               enable debug output.
```

```bash
ucircfull circ_call -m ucRG -i ${sample}_strand.fastq -r $genome -a $gtfFile -u ${sample}_umi.clstr -o ./ucRG -p $sample -t $thread
```

#### output:

- `$sample`.circ.gtf: identification and quantification results of circRNAs
