# kcov

A k-mer coverage analysis tool for variant detection and genome comparison.

## Overview

kcov is a bioinformatics tool that utilizes k-mer based approaches to analyze sequencing data and detect genomic variations. The tool leverages efficient data structures and multi-threading to process large-scale next-generation sequencing (NGS) data.

## Tools

### kcov

A k-mer coverage analysis program that counts k-mers in NGS files and compares them against a reference genome.

**Features:**
- Fast k-mer counting from FASTA/FASTQ files (supports gzipped formats)
- Multi-threaded k-mer counting using pipeline parallelism
- K-mer coverage calculation for reference genome positions
- Output in WIG format for visualization in IGV (Integrative Genomics Viewer)
- Efficient hash table-based k-mer storage with canonical k-mer representation (forward and reverse complement)

**Usage:**
```bash
kcov [options] Reference NGS_files > kmer-coverage.wig
```

**Supported formats:** `*.fq`, `*.fa`, `*.fq.gz`, `*.fa.gz`

**Options:**
- `-k INT`: k-mer size (default: 31)
- `-t INT`: number of threads for k-mer counting (default: 3)
- `-c INT`: minimal k-mer coverage for variant calling (default: 8)
- `-l INT`: maximal assembly length (default: 500)
- `-e FLOAT`: sequencing error rate for P-value calculation (default: 0.025)

**Example commands:**
```bash
# Single-end reads
kcov ref_genome.fa reads.fq > kmer-coverage.wig

# Paired-end reads
kcov ref_genome.fa reads1.fq reads2.fq > kmer-coverage.wig

# Multiple files
kcov ref_genome.fa reads1.fq reads2.fq reads3.fq > kmer-coverage.wig
```

**Technical details:**
- Uses canonical k-mers (minimum of forward and reverse complement)
- Hash table partitioning for parallel processing
- Greedy path search algorithm for variant detection
- Supports insertion, deletion, and substitution variant calling

## Dependencies

### C Program (kcov)
- GCC or compatible C compiler (C11 standard)
- zlib library (for gzip file support)
- pthread library (for multi-threading)
- klib headers (ketopt.h, kthread.h, khashl.h, kseq.h)

## Building

Compile the C program using the provided Makefile:
```bash
make
```

This will generate the `kcov` executable.

## Algorithm

The k-mer analysis pipeline consists of:
1. **K-mer counting**: Extract and count k-mers from NGS reads using multi-threaded pipeline
2. **Hash table storage**: Store k-mers in partitioned hash tables with coverage counts
3. **Reference analysis**: Scan reference genome and query k-mer coverages
4. **Variant detection**: Identify regions with coverage discrepancies indicating variations
5. **Path searching**: Greedy assembly-based approach to reconstruct variant sequences
6. **Output generation**: Produce WIG format coverage tracks and VCF-style variant calls

## Author

Lifu Song (songlf@tib.cas.cn)

## Version

Version 20230301

## License

See LICENSE file for details.
