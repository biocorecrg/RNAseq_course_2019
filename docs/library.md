---
layout: page
title: Library preparation
navigation: 4
---

# Library preparation

* **mRNA purification**. 

mRNA is isolated using a polyT adapter that binds to the polyA tail of RNA. In the result, non poly-adenylated transcripts - rRNA, tRNA, long ncRNAs, miRNA, histone mRNA, degraded RNA, bacterial transcripts, and many viral trascripts - are excluded from the reaction (washed away). 

* **RNA quantification, quality control and fragmentation**


<img src="images/illumina1.png" width="500" align="middle" />

* **RNA is converted to cDNA**


<img src="images/transcription.jpg" width="500" align="middle" />

During a typical RNAseq experiment the information about DNA strands is lost after both strands of cDNA are synthesized. There is however a number of methods for making stranded RNAseq libraries that preserve the strand information. (see for detail, [https://galaxyproject.org/tutorials/rb_rnaseq/](https://galaxyproject.org/tutorials/rb_rnaseq/))

One of such methods (shown below) is implemented in the Illumina's TruSeq Stranded mRNA protocol that uses the introduction of dUTP instead of dTTP during the amplification. The incorporation of dUTP in the second strand synthesis quenches the second strand during amplification, because the polymerase used in the assay is not incorporated past this nucleotide.  
Stranded protocol allows detection of antisense molecules or genes in both 5' and 3' direction. 

<img src="images/stranded_protocol.png" width="500" align="middle" />

Note that in a stranded protocol shown here (in other protocols it can be different), Read 1 is mapped to the antisense strand (this is also true for single-end reads), while Read 2, to the sense strand.

|Read mapping in a stranded vs. unstranded sequencing|
| :---:  |
|<img src="images/stranded_vs_unstranded.jpg" width="900" align="middle" />|
|from [https://galaxyproject.org/tutorials/rb_rnaseq/](https://galaxyproject.org/tutorials/rb_rnaseq/)|

<br/>

* **cDNA multiplexing** 

Fragmented cDNA is indexed with a hexamer or octamer barcode (so that cDNA from different samples can be pooled into a single lane for multiplexed sequencing).

|cDNA multiplexing|
| :---:  |
|<img src="images/multiplexing.jpg" width="600" align="middle" />|
|from [https://github.com/hbctraining/rnaseq_overview](https://github.com/hbctraining/rnaseq_overview)|


* **cDNA amplification** <br/>

<br/>

* **cDNA library quality control and fragment selection** <br/>

<br/>

* **Sequencing**

<br/>

The output of RNA-seq is then demultiplexed yielding either one fastq-file per sample (for single-end reads protocol) or two fastq-files per sample (for paired-end reads protocol).

<br/>

### Experimental design

| Things to consider|
| :---:  |
|<img src="images/exp_design.jpg" width="700" align="middle" />|
|from [https://galaxyproject.org/tutorials/rb_rnaseq/](https://galaxyproject.org/tutorials/rb_rnaseq/)|

<br/>

**HOMEWORK (until tomorrow)**

* Read and do an exercise from [https://github.com/hbctraining/rnaseq_overview/blob/master/lessons/experimental_planning_considerations.md](https://github.com/hbctraining/rnaseq_overview/blob/master/lessons/experimental_planning_considerations.md).
* Read [https://rawgit.com/bioinformatics-core-shared-training/experimental-design/master/ExperimentalDesignManual.pdf](https://rawgit.com/bioinformatics-core-shared-training/experimental-design/master/ExperimentalDesignManual.pdf).

<br/>




