Juan F. Santana, Ph.D. (juan-santana@uiowa.edu), University of Iowa, Iowa City, I.A.
As published in [Santana et al., 2024](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad1253/7513806).

# UMI deduplication

This tool is designed to eliminate duplicate of non-stranded paired-end reads in BED files that contain identical unique molecular identifiers (UMIs). It identifies and flags instances where two or more reads share identical start and end positions. If flagged reads also share the same UMI, they are subsequently removed from the BED file.

## Requirements

### System

- Linux OS

### Data files

- Three FASTQ files per sample: R1, R2, R3 - where R2 correspond to the UMI sequence in the following format:
```
@A00876:427:H3WF7DRX3:1:2101:1181:1000 2:N:0:CACCGG
GTAGGATA
+
FFFFFFFF
```

- A `sampleskey.csv` with the following format (keep the header)
```
name,fastq_1,fastq_2,fastq_3
DMSO-TBP-1,Sample1_R1.fastq.gz,Sample1_R2.fastq.gz,Sample1_R3.fastq.gz
DMSO-TBP-2,Sample1_R1.fastq.gz,Sample1_R2.fastq.gz,Sample1_R3.fastq.gz
```

### Usage

In the working directory, save `sampleskey.csv` and `.fastq.gz` files. 
The bed file name should correspond to a name present in `sampleskey.csv`.
Run:

```
./dedup <path to bed file> <path to samplekey.csv> > <path to output>
```

### Output

- A deduplicated bed file
- Log file containing deduplication statistics as shown below:

| Total reads | Total duplicates | Percentage of reads removed | 
|:-----------:|:----------------:|:---------------------------:|
| 23212949    |  1892565         |                    8.2      |


| Number of duplications for a given UMI | Number of UMIs | Total | 
|:-----------:|:----------------:|:---------------------------:|
| 1           |     1631         |                   1631      |
| 2           |     2162         |                   4324      |
| 3           |     2542         |                   7626      |
| 4           |     2627         |                  10508      |
