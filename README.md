<p align="center">
  <img src="logo.png" alt="ucircfull logo" width="90%" />
</p>

## About

The ucircfull package is a UMI-guided tool for identifying and quantifying circRNAs from UMI-tagged full-length circRNA sequencing (ucircFL-seq) data. It is designed to handle the challenges of low read depth and high sequencing error rates in ucircFL-seq data, and to provide accurate quantification of circRNA abundance.

## Installation

### Install ucircfull binary release:

```bash
wget https://github.com/YueqiJin/ucircfull/releases/download/v1.1.2/ucircfull-1.1.2-Linux.tar.gz
tar -xzf ucircfull-1.1.2-Linux.tar.gz
/path/to/ucircfull-1.1.2-Linux/bin/ucircfull --help
```

### Install ucircfull from source code:

#### Dependencies:

- cmake >= 3.4
- make
- g++ with C++23 support (GCC >= 13 recommended)
- gcc with C++23 support (GCC >= 13 recommended)
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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

#### Install ucircfull to a deployment prefix:

```bash
cmake --install build --prefix /opt/ucircfull
```

This installs executables to `/opt/ucircfull/bin` and shared libraries to `/opt/ucircfull/lib`.

### Run ucircfull in apptainer:

```bash
apptainer pull docker://jinyueqi/ucircfull
apptainer exec /path/to/ucircfull-latest.sif ucircfull --help
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

Extract UMI sequence and identify strand from ucircFL-seq raw fastq

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
  --debug               enable debug output.
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
Usage: ucircfull circ_call [--help] [--version] --mode STR [--notstranded VAR] --input FQ [--bam BAM] --ref REF --anno GTF [--umi CLSTR] --splice MOTIF --outdir DIR --prefix PREFIX --thread INT [--threshold INT] [--minimap2 PATH] [--samtools PATH] [--debug]

circRNA identification.

Optional arguments:
  -h, --help            shows help message and exits
  -v, --version         prints version information and exits
  -m, --mode STR        circRNA calling mode. (RG, ucRG, cRG) [required]
  -sn, --notstranded    find splicing signal in both strand. [default: false]
  -i, --input FQ        stranded fastq file. [required]
  --bam BAM             mapped bam file as input.
  -r, --ref REF         CIRI-long reference directory. [required]
  -a, --anno GTF        reference annotation GTF file. [required]
  -u, --umi CLSTR       umi clust results file. [default: "-"]
  --splice MOTIF        splice motifs. [default: "AGGT,AGGC,ACAT,ACGT,AGAT"]
  -o, --outdir DIR      output directory. [default: "."]
  -p, --prefix PREFIX   output prefix. [default: "circFL"]
  -t, --thread INT      number of threads used. [default: 4]
  --threshold INT       minimum supporting reads for a circRNA transcript. [default: 2]
  --minimap2 PATH       path to minimap2. [default: "minimap2"]
  --samtools PATH       path to samtools. [default: "samtools"]
  --debug               enable debug output.
```

```bash
ucircfull circ_call -m ucRG -i ${sample}_strand.fastq -r $genome -a $gtfFile -u ${sample}_umi.clstr -o ./ucRG -p $sample -t $thread
```

#### output:

- `$sample`.circ.gtf: identification and quantification results of circRNAs
- `$sample`.fusion.txt: detected fusion circRNAs

## Output files

### `$prefix`.circ.gtf

`$prefix`.circ.gtf is a standard GTF file generated by `ucircfull circ_call`. It contains three feature types for each reported circRNA locus:

- `BSJ`: back-splice junction locus, grouping all isoforms from the same BSJ site
- `transcript`: one full-length circRNA isoform (exon composition + splice-site variants)
- `exon`: exon structure of the isoform

The file uses the standard 9-column GTF layout:

| Column | Description |
| :--- | :--- |
| 1 `seqname` | Chromosome name |
| 2 `source` | Always `circfull` |
| 3 `feature` | `BSJ`, `transcript`, or `exon` |
| 4 `start` | 1-based start coordinate |
| 5 `end` | 1-based end coordinate |
| 6 `score` | Always `.` |
| 7 `strand` | Strand of the circRNA isoform |
| 8 `frame` | Always `.` |
| 9 `attribute` | Feature-specific attributes described below |

Attributes written by `ucircfull`:

- `BSJ` records: `gene_id`, `circ_type`, `host_gene_id`, `host_gene_name`, `bsj`
- `transcript` records: `gene_id`, `transcript_id`, `uniform_id`, `circ_type`, `host_gene_id`, `host_gene_name`, `bsj` (per-isoform read count)
- `exon` records: `gene_id`, `transcript_id`, `uniform_id`, `exon_number`

Field meanings:

- `gene_id`: BSJ position ID in the format `chr:start-end`
- `transcript_id`: isoform ID in the format `chr:start-end:strand|exon1_start-exon1_end,exon2_start-exon2_end,...`
- `bsj` on `BSJ` records: total number of supporting reads summed across all isoforms at the BSJ locus
- `bsj` on `transcript` records: number of supporting reads for that isoform
- `circ_type`: circRNA classification (`exon`, `intron`, or `intergenic_region`)
- `host_gene_id`: best-matching host gene identifier from the annotation
- `host_gene_name`: gene symbol of the host gene
- `uniform_id`: standardized circRNA name following the proposed naming scheme; includes exon composition, splice-site variants (L/S), retained introns (RI), novel exons (NE), and an ordered `.N` suffix per distinct BSJ site (e.g. `circAKT3(2,3).1`, `circMCU(2,L3).2`)
- `exon_number`: exon order within the isoform, starting from 1

Only circRNA isoforms with supporting reads at or above the `--threshold` value are output to this file.

### `$prefix`.fusion.txt

`$prefix`.fusion.txt is a tab-separated table describing fusion circRNAs detected by `ucircfull circ_call`.

| Column | Description |
| :--- | :--- |
| `circID` | Fusion circRNA locus ID in the format `chr_first|start_first|end_first|chr_second|start_second|end_second` |
| `isoID` | Fusion isoform ID, including exon composition and strand for both loci |
| `chr_first` | Chromosome of the first fusion locus |
| `start_first` | Start coordinate of the first fusion locus |
| `end_first` | End coordinate of the first fusion locus |
| `len_first` | Length of the first fusion locus |
| `exonNum_first` | Number of exons in the first fusion locus |
| `exon_start_first` | Comma-separated exon start coordinates for the first locus |
| `exon_end_first` | Comma-separated exon end coordinates for the first locus |
| `exon_leftSeq_first` | Left splice-site sequence for each exon in the first locus |
| `exon_rightSeq_first` | Right splice-site sequence for each exon in the first locus |
| `strand_first` | Strand of the first fusion locus |
| `geneName_first` | Gene annotation overlapping the first locus, if available |
| `chr_second` | Chromosome of the second fusion locus |
| `start_second` | Start coordinate of the second fusion locus |
| `end_second` | End coordinate of the second fusion locus |
| `len_second` | Length of the second fusion locus |
| `exonNum_second` | Number of exons in the second fusion locus |
| `exon_start_second` | Comma-separated exon start coordinates for the second locus |
| `exon_end_second` | Comma-separated exon end coordinates for the second locus |
| `exon_leftSeq_second` | Left splice-site sequence for each exon in the second locus |
| `exon_rightSeq_second` | Right splice-site sequence for each exon in the second locus |
| `strand_second` | Strand of the second fusion locus |
| `geneName_second` | Gene annotation overlapping the second locus, if available |
| `readCount` | Number of reads supporting the fusion isoform |
| `readID` | Comma-separated read IDs supporting the fusion isoform |

## Citation
Yueqi Jin, Xueyan Hu, Yun Zhang, Jianqi She, Changyu Tao, Ence Yang, Amplification Optimized and Unique Molecular Identifier Guided High Accuracy Full-length CircRNA Sequencing, Genomics, Proteomics & Bioinformatics, 2026; qzag009, https://doi.org/10.1093/gpbjnl/qzag009