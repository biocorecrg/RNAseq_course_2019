# RNAseq_course_2019
Material for the RNAseq course

## Prior to course
Send survey to participants regarding cluster access and usage.

## The plan (draft)

### Day 1
* Introduction (lecture)
** goals of course
** goals of RNA-seq experiments
** overview of platforms (Illumina, Nanopore, PacBio ?)
** polyA, ribo0, total RNA, microRNA
** stranded, unstranded
** adaptors
** read type: single, paired end; read length.
** experimental / sequencing design
** sequencing depth
* Get raw data from public repository (Drosophila ?): check fastq file
* Adaptor trimming
* FastQC (compare what should be seen in genome versus transcriptome), FastqScreen, MultiQC

### Day 2
* Get reference genome: ENSEMBL, Gencode, UCSC
* Map data to reference genome:
** STAR
** SALMON
* SAM / BAM formats (play with samtools)

### Day 3
* DESeq2: import data from STAR and SALMON
* online tool that integrates DESeq2, edgeR, limma, and more: http://52.90.192.24:3838/rnaseq2g/
* RSEM after SALMON?
* Gene selection
* PCA, heatmap

### Day 4
* Genome browser: ENSEMBL, UCSC or IGV
* bigwig?
* Conversion from UCSC chromosome naming convention to ENSEMBL's ?
* Gene Ontology analysis:
** enrichR
** DAVID ?
** GSEA

### Day 5
* Mini project on small genome: provide link to public data


