# RNAseq_course_2019
Material for the RNAseq course

## Prior to course
Send survey to participants regarding cluster access and usage.

## The plan (draft)
The webpage is at https://biocorecrg.github.io/RNAseq_course_2019/ 

### Day 1
* Introduction (lecture)
  * goals of course
  * goals of RNA-seq experiments
  * overview of platforms (Illumina, Nanopore, PacBio ?)
  * polyA, ribo0, total RNA, microRNA
  * stranded, unstranded
  * adaptors
  * read type: single, paired end; read length.
  * experimental / sequencing design
  * sequencing depth
* Get raw data from public repository (Drosophila ?): check fastq file
* Adaptor trimming
* FastQC (compare what should be seen in genome versus transcriptome), FastqScreen, MultiQC

### Day 2
* Get reference genome or transcriptome: ENSEMBL, Gencode, UCSC
* GTF format; a bit of exploration and excercises 
* Map reads to genome or transcriptome? 
* Existing approaches/methods to read mapping in an RNA-seq experiment: pros and cons, which to choose?
* Map data to reference genome:
  * STAR
  * SALMON
* SAM / BAM formats (play with samtools)
* How to explore BAM files (e.g., using NCBI Genome Workbench https://www.ncbi.nlm.nih.gov/tools/gbench/tutorial6/; UCSC browser https://genome.ucsc.edu/goldenPath/help/bam.html - which http server can be used in the class for this?) other tools to view BAM file SeqMonk, RNAseqViewer, IGB,...
* Make bigwig-files from BAMs and load into GenomeBrowser

I suggest to move Genome Browser here to show its major functionality and to look at BAM and bigwig files. It can be ok to expand it to Day 3. Maybe, move SALMON to Day 3 and compare alignments from STAR ans SALMON using statistics on mapping and looking at bigwig-files in the borwser? 

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
  * enrichR
  * DAVID ?
  * GSEA

### Day 5
* Mini project on small genome: provide link to public data

link to CRG course:  https://public-docs.crg.es/rguigo/courses/rnaseq/2017

