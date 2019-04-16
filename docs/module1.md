---
layout: page
title: Module 1
navigation: 2
---

# Module 1
## Introduction: what is RNA-Seq? 

|RNA molecule|
| :---:  |
|<img src="images/Double-stranded_RNA.gif" width="300" align="middle" />|
|from Wikipedia|


RNA sequencing is a technique that allows to detect and quantify RNA molecules within biological samples by using next-generation sequencing (NGS). In few words this technology is used to analyze the transcriptome by revealing:
* gene expression
* alternative spliced transcripts 
* gene fusion and SNPs
* RNA editing

Other experiments used for assessing the RNA expression are Northern Blot [1], Real time PCR [2] and hybridization-based  microarray [3].

RNA-Seq can be performed on:
* Only messenger RNAs by performing polyA enrichment. 
* total RNA, after a step of ribosomal depletion.
* small RNAs, selecting the size of RNA molecule (usually < 100 nt)
* Only RNA molecules that are being transcribed in that moment (ribosomal profiling).

Depending on the kind of sequencing the RNAseq can produce:
* Single short reads: particularly used for gene quantification
* Paired end reads: useful for splicing detection and annotation refinement
* Stranded or unstranded: the former allows detection of antisense molecules or genes on both directions, the latter is sometimes needed when very little amount of RNA is available. 
* Single long reads (PACBio or Nanopore): used for de novo identification of new transcripts / improving transcriptome. 

## mRNA sequencing
RNA is extracted and then converted to cDNA by using a polyT adapter that binds the polyA tail. In this way you can exclude non poly-adenylated transcripts like rRNA, tRNA and the majority of long ncRNAs. Then     

|RNASeq protocol|
| :---:  |
|<img src="images/tileshop.jpeg" align="middle" />|
|from Wang et al 2009 [4]|



1) https://en.wikipedia.org/wiki/Northern_blot
2) https://en.wikipedia.org/wiki/Real-time_polymerase_chain_reaction
3) https://en.wikipedia.org/wiki/DNA_microarray
4) Wang Z, Gerstein M, Snyder M. RNA-Seq: a revolutionary tool for transcriptomics. Nat Rev Genet. 2009 Jan;10(1):57-63. doi: 10.1038/nrg2484.
