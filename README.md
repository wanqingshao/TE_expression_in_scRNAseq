# Transcript assembly improves expression quantification of transposable elements in single cell RNA-seq data

## Introduction

This document contains a step by step tutorial for analyzing the expression of transposable elements (TEs) using scRNA-seq data.

Comparing to bulk RNA-seq signals, scRNA-seq signals are sparse and noisy. As a result, counting scRNA-seq reads at individual TEs or TE subfamilies/families can aggregate noise and create substantial challenges for downstream analyses. Furthermore, counting reads at individual TEs or TE subfamilies/families fails to take into account the structures of the full-length transcripts, which can consist of multiple TEs from different subfamilies/families. Consequently, different expression values can be assigned to individual TEs within the same transcript. This caveat is especially obvious when dealing with scRNA-seq datasets with signal enrichment at either the 5’ or 3’ end of the RNA (for instance, data generated using Drop-seq or 10x protocols).

To overcome these limitations, we constructed a computational pipeline that quantifies the expression of TE containing non-coding RNAs at the transcript level. In brief, our pipeline uses high quality bulk RNA-seq datasets to perform de novo transcript assembly, selects transcripts that overlap with TEs but not the exons of protein-coding genes (termed TE transcripts), then quantifies scRNA-seq signal at both the protein-coding genes and the assembled TE transcripts.



## Software requirement

* FeatureCounts from subread http://subread.sourceforge.net/
* pigz https://github.com/madler/pigz
* R (>=3.5) https://www.r-project.org/
  packages from bioconductor: GenomicRanges
  packages from CRAN: datatable, optparse, magrittr, parallel, ggplot2, ggrepel, ggpubr
* samtools (>=1.9) http://www.htslib.org/
* STAR (>=2.7) https://github.com/alexdobin/STAR
* StringTie2  https://github.com/skovaka/stringtie2
* TACO https://tacorna.github.io/
* zUMIs original scripts: https://github.com/sdparekh/zUMIs  Updated scripts were provided under the `scripts` folder

## Annotations

* STAR index
* RefSeq annotation
* repeatmasker output
* chromosome length

## Data processing for bulk samples

Bulk RNA-seq data can be used for constructing TE transcript references (Note: in theory, scRNA-seq generated with Smart-seq derived protocols could also be used for transcript assembly. However, from our experience, transcript assembly using scRNA-seq data is challenging and often results in low quality candidates).

### Bulk RNA-seq alignment

Bulk RNA-seq data can be aligned using STAR. We recommend the following parameters (example is for pair end samples):

```
STAR --genomeDir path_to_genome_index \
  --runThreadN 6 \
  --readFilesIn fastq_1.gz fastq_2.gz \
  --outFileNamePrefix output_name \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --outSAMattributes NH HI NM MD XS AS \
  --readFilesCommand zcat
```

`--outFilterMultimapNmax` was set to `500` to increase the number of multiple aligned reads, this is to preserve reads that originated from repetitive regions

`--outSAMattributes` was set to "NH HI NM MD XS AS" to facilitate the downstream transcript assembly step


### Transcript assembly

Transcript assembly was performed using StringTie2. We recommend the following parameters:

```
stringtie input.bam -j 2 -s 5 -f 0.05 -c 2 -p 5 \
  -o output.gtf
```
This set of parameters requires a minimum of 2 reads for junction assembly, 5 reads per bp for single-exon transcript assembly, 2 reads per bp for multi-exon transcript assembly, 5% reads for alternative transcript assembly.

### Checking the quality of bulk RNA-seq samples

Since bulk RNA-seq samples are essential for transcript assembly, we recommend checking the quality of bulk RNA-seq samples before proceeding.

Here are a few heuristic rules that we found useful in determining the quality of bulk RNA-seq samples:

1. Check the percentage of reads that are mapped to the exons of protein-coding genes. High quality RNA-seq samples should have high percentages of protein-coding gene exon derived reads. For a regular polyA selected RNA-seq sample, this percentage should be higher than 70%.

2. Generate bigwig or bedgraph files and examine the reads distribution in genome browser. PolyA selected bulk RNA-seq samples should show clear exon-intron patterns without severe 3' bias.

3. Check the assembly results. High quality RNA-seq samples should yield full-length transcripts with clear exon-intron structures. Transcripts should largely match known annotations at protein-coding genes.


### Merging assembled transcripts

Assembled transcripts can be merged using TACO with the following command:

```
cat *gtf > gtf_to_merge.txt
taco_run -o merged_stringtie -p 10 gtf_to_merge.txt
```

### Downloading and processing annotation files

First download the gene annotation files for the genome you are working on, the example listed below uses the mm10 refGene annotation file hosted by UCSC (http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz).

Download the repeatmasker output file, the example listed below uses the mm10 repeatmasker file hosted by UCSC (http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/rmsk.txt.gz)

Collect chromosome length. This file can be found in the STAR-index folder: chrNameLength.txt

The following script can be used to process downloaded annotations.

```
Rscript export_saf_file.r -g mm10.refGene.txt -r mm10.rmsk.txt -l mm10.chrNameLength.txt -n mm10 -m chrM -p 10
```
The script above takes in the annotation files and outputs processed files in SAF format (The SAF annotation format has five required columns, including GeneID, Chr, Start, End and Strand. These columns can be in any order. More columns can be included in the annotation. Columns are tab-delimited. Column names are case insensitive. See featureCounts document for details. https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts).

`-g` path to refGene.txt file

`-r` path to rmsk.txt file

`-l` path to chrNameLength.txt file

`-m` mitochondria chromosome patterns, chrM or chrMT is the most common mitochondria chromosome names. To specify more than two patterns, separate the two patterns with `|`. For example: `chrM1|chrM2`. If `-m` or `-l` is not set, mitochondria saf files will not be generated.

`-n` basename of the output files. If basename is set to `mm10`, the output files will be `mm10_chrM.saf`, `mm10_nc_exon.saf`, `mm10_pc_exon.saf`, `mm10_te.saf` and `mm10_te.gr.rds`

`-p` number of threads to use, default is set to 5

### Selecting TE transcripts

For this project, we focus on transcripts that contain exonized TE sequences but do not overlap with exons of the protein coding genes. The following script will perform preliminary QC analysis and export the selected TE transcripts. Preliminary QC analysis generates figures showing the types of assembled transcripts and the number of TE transcripts that overlap with known ncRNAs.

```
Rscript processing_assembled_transcripts.r -f assembly.gtf -g mm10_pc_exon.saf -a mm10_nc_exon.saf \
  -m mm10_chrM.saf -r mm10_te.saf -n mm10 -p 10
```

`-f` path to assembled transcript gtf, this is one of TACO's outputs

`-g` path to mRNA exon saf file, file was generated in the previous step

`-a` path to ncRNA exon saf file, file was generated in the previous step

`-m` path to chrM saf file, file was generated in the previous step, default is set as unknown

`-r` path to TE saf file, file was generated in the previous step

`-n` basename for the output files. If basename is set as `mm10`, the outputs will be `mm10_pc_exon_te_tx.saf`(saf file containing mRNA exons, TE transcript exons and mitochondria chromosome if specified), `mm10_te_tx.rds`(rds file containing TE transcripts info as a GRanges object) and mm10_transcript_qc.pdf (pdf file showing the types of assembled transcripts and the number of TE transcripts that overlap with known ncRNAs)

`-p` number of threads to use, default is 5


## Processing scRNA-seq data

Here we provide procedures to process two of the most popular scRNA-seq strategies: scRNA-seq using Smart-seq derived protocols and scRNA-seq using UMI based protocols

### Samples generated with Smart-seq derived protocols

scRNA-seq datasets generated with Smart-seq derived protocols can be processed similarly as bulk RNA-seq data. The example below works for paired-end reads.

#### Sample alignment
```
STAR --genomeDir path_to_genome_index \
  --runThreadN 6 \
  --readFilesIn fastq_1.gz fastq_2.gz \
  --outFileNamePrefix output_name \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterMultimapNmax 500 \
  --outSAMattributes NH HI NM MD XS AS \
  --readFilesCommand zcat \
```

To speed up the process, we recommend loading the STAR index into the memory and process multiple cells in parallel.
To load genome into the memory `STAR --genomeDir path_to_genome_index -- genomeLoad LoadAndExit`. To use the loaded genome `STAR --genomeDir ath_to_genome_index -- genomeLoad LoadAndKeep --limitBAMsortRAM ram_limit`(See STAR document for details).

#### Quantifying signal at TE transcripts

We first use featureCounts to assign each read to its corresponding feature(s).  Reads that overlap with multiple features are then redistribute using an EM-algorithm with uniquely mapped reads as priors.

```
## Assigning reads to features
featureCounts  -F SAF -O --fracOverlap 0.1  -p -B -M --fraction -T 5 \
  -a mm10_pc_exon_te_tx.saf -o Cell1_te_tx_featurecount -R BAM Cell1.bam

## Reformating bam files
samtools view -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x XS \
  Cell1.bam.featureCounts.bam |  grep -e 'XT:Z' > Cell1_featurecount.sam.txt

## Reassigning multi-mapped reads using EM algorithm
Rscripts  redistribute_multiple_aligned_reads.r -f Cell1_featurecount.sam.txt \
  -r mm10_pc_exon_te_tx.saf -n Cell1 -s 50 -m 1 -p 10
```   
The final output for the code above is Cell1_count_dt.txt. This file contains four columns: "feature" is the feature name (genes and TE transcripts), "uniq" quantifies the number of reads that are mapped to single features, "EM_distri" quantifies the number of reads at each feature after distributing reads that are mapped to multiple features using EM algorithm, "even_distri" quantifies the number of reads at each feature after evenly distributing reads that are mapped to multiple features.

For script `redistribute_multiple_aligned_reads.r`:

`-f` path to reformatted alignment files

`-r` path to features in saf format

`-n` basename of the final output

`-s` maximum number of cycles for the EM algorithm, default value is 50

`-m` stop the EM algorithm if the maximum number of reads changes per feature is less than m, default value is set as 1

`-p` number of threads, default value is 5

### Samples generated using UMI based protocols

scRNA-seq datasets generated using UMI based protocols (10x, Drop-Seq, CEL-seq, MARS-seq...) are processed using zUMIs to generate the reformatted bam files. We have modified the original zUMIs scripts to facilitate TE expression quantification. The modified zUMI package can be found under the `scripts` folder.

#### Sample processing

```
zUMIs-master-TE.sh -y zUMIs_sample.yaml -d path_to_modified_zUMIs_package

```
The scripts above follows the design of the original zUMIs. Key parameters are set using the yaml file. The outputs contain `Sample_1_featurecount.sam.txt`, which can be used for signal quantification at the TE transcripts.


#### Quantifying signal at TE transcripts

```
scRNA_UMI_counting.R -f Sample_1_featurecount.sam.txt -r mm10_pc_exon_te_tx.saf -n Sample_1 -s 50 -m 1 -p 10

```
The final output for is Sample_1_count_dt.txt. This file contains five columns:  "feature" is the feature name (genes and TE transcripts), "uniq" quantifies the number of reads that are mapped to single features, "EM_distri" quantifies the number of reads at each feature after distributing reads that are mapped to multiple features using EM algorithm, "even_distri" quantifies the number of reads at each feature after evenly distributing reads that are mapped to multiple features, "RG" specifies the cell barcode.
