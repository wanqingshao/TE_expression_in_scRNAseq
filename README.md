# Leveraging transcript assembly for transposable elements expression quantification

## Introduction

This document contains a step by step tutorial for analyzing the expression of transposable elements (TE) by incorporating.

Comparing to bulk RNA-seq, scRNA-seq data are often sparse and noisy. As a result, directly counting scRNA-seq reads at individual TEs or TE subfamilies/families can aggregate noise and create substantial challenges for downstream analysis. Furthermore, counting reads at individual TEs or subfamilies/families fails to take into account the structures of the full-length transcripts, which can consist of multiple TEs from different subfamilies/families. Consequently, different expression values will be assigned to individual TEs within the same transcript. This caveat is especially obvious when dealing with scRNA-seq datasets where sequencing reads are enriched at either the 5’ or 3’ end of the RNA.

To overcome these limitations, we constructed a computational pipeline that quantifies the expression of non-coding RNAs that contain TE derived sequences. Instead of counting reads of individual TEs, our pipeline assembles transcripts that contains TE derived sequences and quantify their expression at the transcript level. In brief, our pipeline constructs transcripts using high quality bulk RNA-seq datasets, selects transcripts that overlap TEs but not the exons of protein-coding genes. These transcripts are termed TE transcripts. Our pipeline then quantifies scRNA-seq signal at both the protein-coding genes and the assembled TE transcripts.



## Software requirement

* FeatureCount from subread http://subread.sourceforge.net/
* pigz https://github.com/madler/pigz
* R (>=3.5) https://www.r-project.org/
  packages from bioconductor: GenomicRanges
  packages from CRAN: datatable, optparse, magrittr, parallel, ggplot2, ggrepel, ggpubr
* samtools (>=1.9) http://www.htslib.org/
* STAR (>=2.7) https://github.com/alexdobin/STAR
* Stringtie2  https://github.com/skovaka/stringtie2
* Taco https://tacorna.github.io/
* zUMIs original scripts: https://github.com/sdparekh/zUMIs  Updated scripts were provided under the `scripts` folder

## Annotations

* STAR index
* RefSeq annotation
* repeatmasker output
* chromosome length

## Data processing for bulk samples

Bulk RNA-seq data can be used for constructing TE transcript references (Note: in theory, scRNA-seq generated with Smart-seq related protocols could also be used for transcript assembly, however, from our experience, transcript assembly using scRNA-seq data is challenging and often results in low quality candidates).

### Bulk RNA-seq alignment

Bulk RNA-seq can be aligned with STAR, we recommend the following parameters (example is for pair end samples):

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

`--outFilterMultimapNmax` was set to "500" to increase the number of multiple aligned reads, this is to preserve reads that originated from repetitive regions

`--outSAMattributes` was set to "NH HI NM MD XS AS" facility the downstream transcript assembly step


### Transcript assembly

Transcript assembly can be performed with StringTie2, we recommend the following parameters:

```
stringtie input.bam -j 2 -s 5 -f 0.05 -c 2 -p 5 \
  -o output.gtf
```
This set of parameters requires 2 reads for junction assembly, 5 reads per bp for single-exon transcript assembly, 2 reads per bp for multi-exon transcript assembly, 1% reads for alternative transcript assembly.

### Checking the quality of bulk RNA-seq samples

Since the quality of bulk RNA-seq samples has a significant influence on the quality of assembled transcripts, we recommend checking the quality of bulk RNA-seq samples before moving forwards.

Here are a few heuristic rules we found useful in determining the quality of bulk RNA-seq samples

1. The percentage reads that mapped to the exons of protein-coding genes. High quality RNA-seq samples often have high percentage of protein-coding gene exon derived reads. For a regular polyA selected RNA-seq sample, this percentage should be > 70%
2. Generate bigwig or bedgraph files and examine the reads distribution in genome browser. For a regular polyA selected RNA-seq sample, reads should show clear exon-intron pattern without severe 3' signal bias.
3. Check the assembly results. High quality RNA-seq samples should yield full-length transcripts with clear exon-intron structures, transcripts should largely match known annotations at protein-coding genes.


### Merging assembled transcripts

Assembled transcripts can be merged with Taco using the following command

```
cat *gtf > gtf_to_merge.txt
taco_run -o merged_stringtie -p 10 gtf_to_merge.txt
```

### Downloading and processing annotation files

First download the gene annotation files for the genome you are working on, the example listed below uses the mm10 refGene annotation file hosted by the UCSC (http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz).

Then download the repeatmasker frile, the example listed below uses the mm10 repeatmasker file hosted by UCSC (http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/rmsk.txt.gz)

Collect chromosome length, this file can be found in the STAR-index folder: chrNameLength.txt

The following script can be used to process downloaded annotations

```
Rscript export_saf_file.r -g mm10.refGene.txt -r mm10.rmsk.txt -l mm10.chrNameLength.txt -n mm10 -m chrM -p 10
```
The script above takes in the annotation files and output processed files in SAF format (The SAF annotation format has five required columns, including GeneID, Chr, Start, End and Strand. These columns can be in any order. More columns can be included in the annotation. Columns are tab-delimited. Column names are case insensitive, see featurecount document for details https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts).

`-g` path to refGene.txt file

`-r` path to rmsk.txt file

`-l` path to chrNameLength.txt file

`-m` mitochondria chromosome patterns, chrM or chrMT is the most common mitochondria chromosome names, if specify more than two pattern, separate the two pattern with `|`  for example `chrM1|chrM2`. If `-m` or `-l` is not set, mitochondria saf files will not be generated.

`-n` basename for the output files. If basename is set to `mm10`, the outputs will be `mm10_chrM.saf`, `mm10_nc_exon.saf`, `mm10_pc_exon.saf`, `mm10_te.saf` and `mm10_te.gr.rds`

`-p` number of threads to use, default is 5

### Selecting TE transcripts

For the scope of this project, we focus on transcripts that contain exonized TE sequences but do not overlap the exons the protein coding genes. The following script will perform some preliminary QC analysis and export the selected TE transcripts. Preliminary QC analysis including figures showing the type of assembled transcripts and the number of TE transcripts that overlap with known ncRNA genes.

```
Rscript processing_assembled_transcripts.r -f assembly.gtf -g mm10_pc_exon.saf -a mm10_nc_exon.saf \
  -m mm10_chrM.saf -r mm10_te.saf -n mm10 -p 10
```

`-f` path to assembled transcript gtf, this is one of taco's outputs

`-g` path to mRNA exon saf file, file was generated in the previous step

`-a` path to ncRNA exon saf file, file was generated in the previous step

`-m` path to chrM saf file, file was generated in the previous step, default unknown

`-r` path to TE saf file, file was generated in the previous step

`-n` basename for the output files. If basename is set as `mm10`, the outputs will be `mm10_pc_exon_te_tx.saf`(saf file containing mRNA exons, TE transcript exons and mitochondria chromosome if specified), `mm10_te_tx.rds`(rds file containing TE transcripts info as a GRanges object) and mm10_transcript_qc.pdf (pdf file showing the types of assembled transcripts and the number of TE transcripts that overlap with known ncRNA genes)

`-p` number of threads to use, default is 5


## Processing scRNA-seq data

We provide processing procedures for two of the most popular scRNA-seq strategies: scRNA-seq using Smart-seq derived protocols and scRNA-seq using UMI based protocols

### Samples generated with Smart-seq derived protocols

scRNA-seq datasets generated with Smart-seq derived protocols can be processed similarly as bulk RNA-seq data, the example below works for paired-end reads.

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

To speed up the processing process, we recommend loading the STAR index into the memory and process multiple cells in parallel.
To load genome into the memory `STAR --genomeDir path_to_genome_index -- genomeLoad LoadAndExit`. To use the loaded genome `STAR --genomeDir ath_to_genome_index -- genomeLoad LoadAndKeep --limitBAMsortRAM ram_limit`(See STAR document for details).

#### Quantifying signal at TE transcripts

We first use featureCount to assign each read to its corresponding features, multiple mapped reads were then redistribute with an EM-algorithm using uniquely mapped reads as priors.

```
## Assign reads to features
featureCounts  -F SAF -O --fracOverlap 0.1  -p -B -M --fraction -T 5 \
  -a mm10_pc_exon_te_tx.saf -o Cell1_te_tx_featurecount -R BAM Cell1.bam

## reformat bam files
samtools view -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x XS \
  Cell1.bam.featureCounts.bam |  grep -e 'XT:Z' > Cell1_featurecount.sam.txt

## reassign multi-mapped reads using EM algorithm
Rscripts  redistribute_multiple_aligned_reads.r -f Cell1_featurecount.sam.txt \
  -r mm10_pc_exon_te_tx.saf -n Cell1 -s 50 -m 1 -p 10
```   
The final output for the code above is Cell1_count_dt.txt, this file contains four columns, "feature" is name of feature, "uniq" quantifies the number of reads that were only mapped to single features, "EM_distri" quantifies the number of reads at features using EM algorithm, "even_distri" quantifies the number of reads at feature by evenly distribute reads that were mapped to multiple features.

For script `redistribute_multiple_aligned_reads.r`:

`-f` path to reformatted alignment files

`-r` path to features in saf format

`-n` basename of the final output

`-s` maximum number of cycles for the EM algorithm, default value 50

`-m` stop the EM algorithm if the maximum number of reads changes per feature is less than m, default value 1

`-p` number of threads, default value 5

### Samples generated using UMI based protocols

scRNA-seq datasets generated with UMI based protocols (10x, Drop-Seq, CEL-seq, MARS-seq...) were processed using zUMIs to generate the reformatted bam files. We have modified the original zUMIs scripts to facilitate TE expression quantification. The modified zUMI package can be found under the `scripts` folder.

#### Sample processing

```
zUMIs-master-TE.sh -y zUMIs_sample.yaml -d path_to_modified_zUMIs_package

```
The scripts above follows the design of the original zUMIs, key parameters were set using the yaml file. The outputs contains `Sample_1_featurecount.sam.txt`, which could be used for signal quantification at the TE transcripts.


#### Quantifying signal at TE transcripts

```
scRNA_UMI_counting.R -f Sample_1_featurecount.sam.txt \
-r mm10_pc_exon_te_tx.saf -n Sample_1 -s 50 -m 1 -p 3

```

`-f` path to reformatted alignment files

`-r` path to features in saf format

`-n` basename of the final output

`-s` maximum number of cycles for the EM algorithm, default value 50

`-m` stop the EM algorithm if the maximum number of reads changes per feature is less than m, default value 1

`-p` number of threads, default value 5

The final output for the code above is Sample_1_count_dt.txt. This file contains five columns, "feature" is name of feature, "uniq" quantifies the number of reads that were only mapped to single features, "EM_distri" quantifies the number of reads at features using EM algorithm, "even_distri" quantifies the number of reads at feature by evenly distribute reads that were mapped to multiple features, "RG" specifies the cell barcode.
