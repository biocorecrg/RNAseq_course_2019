---
layout: page
title: Module 2
navigation: 3
---

# Module 2

## Read alignment to reference transcriptome

|Mapping of short reads|
| :---:  |
|<img src="images/800px-Mapping_Reads.png" width="350" align="middle" />|
|from Wikipedia|

Once we tested the quality of our input reads we can proceed with the alignment to the reference transcriptome. In last 10 years, several tools for mapping have been developed for this purpose. All of them requires to build an index of the reference genome / transcriptome using different algorithms.

### Fast aligners 
These tools can be used for aligning short reads to a trancriptome or genome reference, but in the latter case they cannot align in splicing junctions. They can be much faster than traditional aligners like [**Blast**](https://blast.ncbi.nlm.nih.gov/Blast.cgi) but less sensitive and have limitations due to the read size. 

1. [**Bowtie**](http://bowtie-bio.sourceforge.net/index.shtml) is an ultrafast, memory-efficient short read aligner geared toward quickly aligning large sets of short DNA sequences (reads) to large genomes. Bowtie indexes the genome with a Burrows-Wheeler index. 
2. [**Bowtie2**](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index. 
3. [**BWA**](http://bio-bwa.sourceforge.net/) is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. BWA indexes the genome with an FM Index.
4. [**GEM**](https://github.com/smarco/gem3-mapper) is a high-performance mapping tool for aligning sequenced reads against large reference genomes. In particular, it is designed to obtain best results when mapping sequences up to 1K bases long. GEM3 indexes the reference genome using a custom FM-Index design and performs an adaptive gapped search based on the characteristics of the input and the user settings. 

### Splice-aware aligners 
These aligners are able to map to the splicing junctions described in the annotation and even to detect novel ones. Some of them can detect gene fusions and SNPs and RNA editing.

1. [**Tophat**](https://ccb.jhu.edu/software/tophat/index.shtml) is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons.
2. [**Hisat2**](http://ccb.jhu.edu/software/hisat2/index.shtml) is a fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes (as well as to a single reference genome). The reference transcriptome is indexed using a Hierarchical Graph FM index (HGFM). 
3. [**STAR**](https://github.com/alexdobin/STAR) is an ultrafast universal RNA-seq aligner. It uses sequential maximum mappable seed search in uncompressed suffix arrays followed by seed clustering and stitching procedure. It is also able to search for gene fusions.

### Quasi-mapping aligners 
These tools are way faster than the previous ones because they don't need to report the resulting alignments but only  associate a read to a given transcript for quantification.  

1. [**Salmon**](https://salmon.readthedocs.io/en/latest/index.html) is a tool for wicked-fast transcript quantification from RNA-seq data. It requires a set of target transcripts to quantify and a K-mer parameter to make the index (i.e. minimum acceptable alignment). 
2. [**Kallisto**](https://pachterlab.github.io/kallisto/) is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.


