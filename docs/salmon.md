---
layout: page
title: Mapping with Salmon
navigation: 14
---


# Mapping using Salmon

<img src="images/RNAseq_workflow.png" width="1000"/>

[**Salmon**](https://combine-lab.github.io/salmon/) is a tool for quantifying the expression of transcripts using RNA-seq data. It is a quasi-mapper as it doesn't produce the read alignment (and doesn't output BAM/SAM files).

For indexing with **Salmon** we need to use transcripts sequences in a fasta file. **Salmon** does not need any decompression of the input so we can index by using this command:

```{bash}
$RUN salmon index -t annotations/gencode.v29.transcripts.fa.gz -i indexes/transcripts

Version Info: This is the most recent version of salmon.
index ["transcripts"] did not previously exist  . . . creating it
[2019-04-30 18:12:59.272] [jLog] [info] building index
[2019-04-30 18:12:59.275] [jointLog] [info] [Step 1 of 4] : counting k-mers

[....]
[2019-04-30 18:18:07.251] [jLog] [info] done building index
```

For aligning with **Salmon** we need to specify the strandess of the library (**Fragment Library Types**). In brief you have to specify three letters:

**The first:**

|Symbol |Meaning | Reads|  
| :---: | :----: |:----: |
|I|inward|-> ... <- |
|O|outward|<- ... ->|
|M|matching|> ... ->|

**The second:**

|Symbol |Meaning |
| :---: | :----: |
|S|stranded|
|U|unstranded|

**The third:**

|Symbol |Meaning |
| :---: | :----: |
|F|read 1 (or single-end read) comes from the forward strand|
|R|read 1 (or single-end read) comes from the reverse strand|

In our case we have **Inward**, **Stranded** and **Reverse**. Moreover if we want to assign the reads to the genes too we need to provide a GTF file with correlation between transcripts and genes (option **-g**).

```{bash}
$RUN salmon quant -i indexes/transcripts -l ISR \
    -1 resources/A549_0_1chr10_1.fastq.gz \
    -2 resources/A549_0_1chr10_2.fastq.gz \
    -o alignments/salmon_A549_0_1 \
    -g annotations/gencode.v29.annotation_chr10.gtf \
    --seqBias \
    --validateMappings 

Version Info: This is the most recent version of salmon.
### salmon (mapping-based) v0.13.1
### [ program ] => salmon 
### [ command ] => quant 
....
```

We can check the results inside the folder 

```{bash}
ls alignments/salmon_A549_0_1/

aux_info  cmd_info.json  lib_format_counts.json  libParams  logs  quant.genes.sf  quant.sf
```

And in particular the file **quant.genes.sf**, that is a tabular file with the following information:


|Column |Meaning |   
| :----: | :---- |
|Name| Gene name|
|Length| Gene length|
|EffectiveLength| Effective length after considering biases|
|TPM|Transcripts Per Million|
|NumReads|Estimated number of reads considering both univocally and multimapping reads|

```{bash}
head -n 5 alignments/salmon_A549_0_1/quant.genes.sf 

Name	Length	EffectiveLength	TPM	NumReads
ENSG00000285803.1	1152	1116.21	7.96961	15.764
ENSG00000285712.1	1590	1545.58	2.19064	6
ENSG00000285824.1	1120	860.855	6.91601	10.551
ENSG00000285884.1	790	515.683	3.28285	3

head -n 5 alignments/salmon_A549_0_1/quant.sf 
Name	Length	EffectiveLength	TPM	NumReads
ENST00000016171.5	2356	1970.742	659.861626	2304.468
ENST00000020673.5	4183	5925.497	0.000000	0.000
ENST00000173785.4	925	868.802	0.000000	0.000
ENST00000181796.6	3785	3216.057	0.000000	0.000
```
We will use this information for calculating differential expression (DE) analysis. 

# Combining reports
At this point we can summarize the work done by using the tool [**multiqc**](https://multiqc.info/). First we link our mapping results to QC.

```{bash}
cd QC/
ln -s ../alignments/* .
```

Then we join the different analyses:

```{bash}
$RUN multiqc .
[INFO   ]         multiqc : This is MultiQC v1.7 (7d02b24)
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching 'QC/'
Searching 70 files..  [####################################]  100% 
...

firefox multiqc_report.html
```

Here the result:

<img src="images/multiqc.png"  align="middle" />

<br/>
