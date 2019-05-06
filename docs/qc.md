---
layout: page
title: Quality control
navigation: 9
---

# QC of sequencing reads
To assess the quality of sequencing data, we will use the programms [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [**Fastq Screen**](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/). 

FastQC calculates statistics about the composition and quality of raw sequences, while Fastq Screen looks for possible contaminations. 

```{bash}
mkdir QC
$RUN fastqc resources/A549_25_3chr10_*.fastq.gz -o ./QC/

Started analysis of A549_25_3chr10_1.fastq.gz
Approx 5% complete for A549_25_3chr10_1.fastq.gz
Approx 10% complete for A549_25_3chr10_1.fastq.gz
Approx 15% complete for A549_25_3chr10_1.fastq.gz
Approx 20% complete for A549_25_3chr10_1.fastq.gz
...
Approx 85% complete for A549_25_3chr10_2.fastq.gz
Approx 90% complete for A549_25_3chr10_2.fastq.gz
Approx 95% complete for A549_25_3chr10_2.fastq.gz
Analysis complete for A549_25_3chr10_2.fastq.gz
```

We can display the results with a browser; e.g., Firefox, for each file individually or all files with one command:
```{bash}
firefox QC/A549_25_3chr10_1_fastqc.html

firefox QC/*.html
```

<img src="images/fastqc.png" width="800"/>

<br/>

Below is an example of a poor quality dataset. As you can see, the average quality drops dramatically towards the 3'-end.

<img src="images/bad_fastqc.png" width="800"/>


**Fastq Screen** requires a number of databases to be installed. A number of prepared by the authors of the program databases can be downloaded by using the following command (**DO NOT LAUNCH IT NOW!!!**):

```{bash}
$RUN fastq_screen --get_genomes
``` 

This will download 11 genomes (arabidopsis, drosophila, E. coli, human, lambda, mouse, mitochondria, phiX, rat, worm and yeast) and 3 collection of sequences (adapters, vectors, rRNA) indexed with bowtie2. This step is quite slow so we are not going to launch it now. The files will be stored inside **FastQ_Screen_Genomes** folder.
We have to modify the fastq_screen.conf file using the real path of the executable **/usr/local/bin/bowtie2** and the one of the genome.

To execute Fastq Screen: 

```{bash}
$RUN fastq_screen --conf FastQ_Screen_Genomes/fastq_screen.conf \
           resources/A549_0_1chr10_1.fastq.gz \
           --outdir ./QC/A549_0_1

Using fastq_screen v0.13.0
Reading configuration from 'fastq_screen.conf'
Aligner (--aligner) not specified, but Bowtie2 path and index files found: mapping with Bowtie2
Adding database Human
Adding database Mouse
Adding database Rat
Adding database Drosophila
Adding database Worm
Adding database Yeast
Adding database Arabidopsis
Adding database Ecoli
Adding database rRNA
Adding database MT
Adding database PhiX
Adding database Lambda
Adding database Vectors
Adding database Adapters
Using 7 threads for searches
Option --subset set to 100000 reads
Processing A549_0_1chr10_1.fastq.gz
Counting sequences in A549_0_1chr10_1.fastq.gz
Making reduced sequence file with ratio 711:1
...
```

Below is shown an example of the result aligning reads from fastq_screen.conf A549_0_1_1.fastq.gz to different databases to detect contaminations, failure of ribosomal depletion, and other artifacts of sequencing.  

<img src="images/A549_0_1_1_screen_2.png" />
<img src="images/A549_0_1_1_screen_1.png" />

<br/>

## Initial processing of sequencing reads

Before mapping reads to the genome/transcriptome or performing a de novo assembly, the reads must be pre-processed as follows: 
* Demultiplex by index or barcode (it is usually done in the sequencing facility)
* Remove adapter sequences
* Trim reads by quality
* Discard reads by quality/ambiguity
* Filter reads by k-mer coverage (recommended for the de novo assembly)
* Normalize k-mer coverage (recommended for the de novo assembly)

As shown before, both the presence of low quality reads and adapters are reported in the **fastqc** output. 

Adapters are usually expected in small RNA-Seq because the molecules are tipically shorter than the reads, and that makes an adapter to be present at 3'-end. Let's run FastQC on a fastq file for small RNA-Seq.


```{bash}
$RUN fastqc resources/subsample_to_trim.fq.gz -o ./QC

Started analysis of subsample_to_trim.fq.gz
Approx 5% complete for subsample_to_trim.fq.gz
Approx 10% complete for subsample_to_trim.fq.gz
Approx 15% complete for subsample_to_trim.fq.gz
Approx 20% complete for subsample_to_trim.fq.gz
Approx 25% complete for subsample_to_trim.fq.gz
...
```


<img src="images/fastqc_small_rnas.png" width="800"/>

<br/>

## EXERCISE

Once you got the FastQC report (above), how to find out the sequence(s) of the adapter(s) that needs to be trimmed?



There are many tools for trimming reads and removing adapters, such as **Trim Galore!**, **Trimmomatic**, **Cutadapt**, [**skewer**](https://github.com/relipmoc/skewer), **AlienTrimmer**, **BBDuk**, and the most recent **SOAPnuke** and [**fastp**](https://www.ncbi.nlm.nih.gov/pubmed/30423086). 

Let's use **skewer7** to trim the Illumina small RNA 3' adapter.  

```{bash}
$RUN skewer resources/subsample_to_trim.fq.gz -x TGGAATTCTCGGGTGCCAAGG -o QC/subsample_to_trim

.--. .-.
: .--': :.-.
`. `. : `'.' .--. .-..-..-. .--. .--.
_`, :: . `.' '_.': `; `; :' '_.': ..'
`.__.':_;:_;`.__.'`.__.__.'`.__.':_;
skewer v0.2.2 [April 4, 2016]
Parameters used:
-- 3' end adapter sequence (-x):	TGGAATTCTCGGGTGCCAAGG
-- maximum error ratio allowed (-r):	0.100
-- maximum indel error ratio allowed (-d):	0.030
-- minimum read length allowed after trimming (-l):	18
-- file format (-f):		Sanger/Illumina 1.8+ FASTQ (auto detected)
-- minimum overlap length for adapter detection (-k):	3
Thu Apr 18 17:51:18 2019 >> started
|=================================================>| (100.00%)
Thu Apr 18 17:51:25 2019 >> done (6.789s)
1000000 reads processed; of these:
  30171 ( 3.02%) short reads filtered out after trimming by size control
   2220 ( 0.22%) empty reads filtered out after trimming by size control
 967609 (96.76%) reads available; of these:
 958360 (99.04%) trimmed reads available after processing
   9249 ( 0.96%) untrimmed reads available after processing
log has been saved to "QC/subsample_to_trim-trimmed.log".
```

We can look at the read distribution after the trimming of the adapter by inspecting the log-file or re-launching FastQC.

```{bash}
$RUN fastqc QC/subsample_to_trim-trimmed.fastq -o QC

Started analysis of subsample_to_trim-trimmed.fastq
Approx 5% complete for subsample_to_trim-trimmed.fastq
Approx 10% complete for subsample_to_trim-trimmed.fastq
Approx 15% complete for subsample_to_trim-trimmed.fastq
Approx 20% complete for subsample_to_trim-trimmed.fastq
Approx 25% complete for subsample_to_trim-trimmed.fastq
Approx 30% complete for subsample_to_trim-trimmed.fastq
...

```
<img src="images/size_dist_small.png" width="800"/>
<img src="images/adapter_removed.png" width="800"/>



## EXERCISE

Let's explore the tool **skewer** in more detail, using "skewer --help" command.
* Which parameter indicates the minimum read length allowed after trimming? And what is its default value?
* Which parameter indicates the threshold on the average read quality to be filtered out?
* Using skewer filter out reads in "subsample_to_trim.fq-trimmed.fastq" that have average quality below 30 and trim them on 3'-end until the base quality is reached 30. How many reads were filtered out and how many remained?



