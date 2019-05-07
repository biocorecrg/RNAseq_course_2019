---
layout: page
title: Data repositories
navigation: 8
---

# RNA-Seq data repositories


The major repositories for gene expression data:
* [**GEO**](https://www.ncbi.nlm.nih.gov/geo/) 
* [**Array-express**](https://www.ebi.ac.uk/arrayexpress/)
* [**ENCODE**](https://www.encodeproject.org/)

These repositories  are linked to the repositories of NGS raw data (Fastq files):
* [**SRA**](https://www.ncbi.nlm.nih.gov/sra) (Sequence Read Archive) 
* [**ENA**](https://www.ebi.ac.uk/ena) (European Nucleotide Archive) 
* [**DDBJ-DRA**](https://www.ddbj.nig.ac.jp/dra/index-e.html) 

<br/>

## EXERCISE
Let's explore [one of the GEO records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126535)
* Which platform and protocol were used for sequencing?
* What type of RNA was sequenced?
* How many samples were sequenced?

<br/>

**NOTE: You will need to download data from SRA for an independent project after this week!**
To download raw data from **SRA**, it is possible to use **fastq-dump program** from [**SRA toolkit**](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) or to download files from the NCBI ftp website using wget (for detail, see [https://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.when_to_use_a_command_line](https://www.ncbi.nlm.nih.gov/books/NBK158899/#SRA_download.when_to_use_a_command_line).

To download data, use the SRA identifier specifying whether reads are single or paired-end, otherwise paired-end reads will be downloaded as a single interleaved file; for paired-end reads use the parameter **--split-files in fastq-dump**. Fastq-dump adds SRA ID to each read in the file, to avoid it, use the parameter **--origfmt**. The parameter **--gzip** compresses fastq files immediately after download. The command below will download fastq-file(s) for one sample only (for example, using SRR identifier SRR8571764 from the exercise above; it is slow - it might take up to 30-40 minutes):

```{bash}
fastq-dump --gzip --origfmt --split-files SRR-IDENTIFIER
```

To download all samples for a specific GEO experiment, use the SRA study identifier (e.g., for the GEO experiment considered above, it is SRP185848) and follow the steps:
* First, download a list of SRR identifiers for all samples in the study by going to [the NCBI SRA page for this study](https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=522280) and clicking on the right top "Send" --> "File" --> "Accession List" --> "Save to file". That will give you the text file with all SRR identifiers for this study; save it for example to the file "sra_ids.txt". 
* Second, run the following command:

```{bash}

fastq-dump --gzip --origfmt --split-files $(<sra_ids.txt)

```
<br/>

Another source of high quality data on gene expression in human and mouse is [The Encyclopedia of DNA Elements (ENCODE)](https://www.encodeproject.org/). Using the ENCODE portal one can access data produced by the ENCODE Consortium.


<br/>
