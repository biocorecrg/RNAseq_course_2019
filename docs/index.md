---
layout: page
title: Home
navigation: 1
---

# Biocore's RNAseq course 2019

## Dates, time & location

Dates:

| Day  | Date  | Time  |
| :---:  | :---  | ---:  |
| 1 | Tuesday 13th of May 2019|9:30-13:30|
| 2 | Thursday 14th of May 2019|9:30-13:30| 
| 3 | Thursday 15th of May 2019|14:00-17:00| 
| 4 | Thursday 16th of May 2019|14:00-17:00| 
| 5 | Thursday 27th of May 2019|9:30-13:00| 


Location:
Carrer del Dr. Aiguader, 88, 08003 Barcelona.
CRG Training center, PRBB building ground floor. 

<iframe src="https://www.google.com/maps/embed?pb=!1m14!1m8!1m3!1d11973.94726186489!2d2.1942455!3d41.3852331!3m2!1i1024!2i768!4f13.1!3m3!1m2!1s0x0%3A0x81e449abea5aae0e!2sPRBB+Parc+de+Recerca+Biom%C3%A8dica+de+Barcelona!5e0!3m2!1sit!2ses!4v1551808726678" width="500" height="500" align="middle" frameborder="0" style="border:0" allowfullscreen></iframe>

<br/>
<br/>

## Instructors

|[Luca Cozzuto](mailto:luca.cozzuto@crg.eu)| [Julia Ponomarenko](mailto:julia.ponomarenko@crg.eu)  | [Sarah Bonnin](mailto:sarah.bonnin@crg.eu) |
| :---:  | :---:  | :---:  |
|<a href="https://biocore.crg.eu/wiki/User:Lcozzuto"><img src="pics/lcozzuto.jpg" width="100"/> </a> |<a href="https://biocore.crg.eu/wiki/User:Jponomarenko"><img src="pics/ponomarenko.jpg" width="100"/> </a> |<a href="https://biocore.crg.eu/wiki/User:SBonnin"><img src="pics/sbonnin.jpg" width="100"/></a> | 


from the CRG [Bioinformatics core facility](https://biocore.crg.eu/) (office 460, 4th floor hotel side)

Material available at https://biocorecrg.github.io/RNAseq_course_2019/

##  Learning objectives
At the end of the course, the participants will be able to:
*	Understand the steps from raw reads to expression counts, differential expression and interpretation of gene lists using enrichment analysis
*	Define a good experimental design, including experimental design, sequencing design, and quality control steps)
*	Perform quality assessment of RNA-seq data, raw and processed
*	Understand file formats commonly used in RNA-seq data analysis
*	Gain an overview on common software tools for RNA-seq data analysis and their limitations
*	Run RNA-seq pipeline to perform differential expression analysis


##  Course Program
*	Day 1: Introduction, download raw data, quality control of data
*	Day 2: Read mapping to reference genome
*	Day 3: Differential expression analysis
*	Day 4: Genome browser, gene ontology enrichment analysis
*	Day 5: Troubleshooting of a mini-project


## Prerequisites
[Basic experience of command line Linux](https://biocorecrg.github.io/advanced_linux_2019/)

## Basic concepts on Linux Containers
During this course we will use a number of tools that have been stored within a [**Linux Container**](https://en.wikipedia.org/wiki/LXC). 

A Linux Container can be seen as minimal virtual environment that can be used in any Linux compatible machine, this allows us to save a lot of time / resources related to installation of tools and libraries and improves the reproducbility of an analysis. In brief you will nedd to install (or ask your IT to have installed) just one program.

In particular we created a [**Docker**](https://www.docker.com/) image and stored it in [DockerHub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/rnaseq2019). 

This image can be downloaded and used in computers running Linux or Mac OS or can be converted into a another Linux Container called [**Singularity**](https://www.sylabs.io/docs/) that is the one that will be used in this course. 

The Singularity image is basically a file that can be accessed by the program Singulairty for executing any software is installed within that image. You can create this image to use the tools we will show in this course at your convenience by using the following command:

```{bash}
singularity pull docker://biocorecrg/rnaseq2019:X.X
```

And access it in the following way 

```{bash}
export RUN="singularity exec ./rnaseq2019-1.0.simg"
$RUN STAR --help
```

If you are a CRG user with access to our cluster you can just access it at this path:

```{bash}
export RUN="singularity exec $PATH/rnaseq2019-1.0.simg"
$RUN STAR --help
```






