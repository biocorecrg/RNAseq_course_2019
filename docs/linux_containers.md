---
layout: page
title: Linux Containers
navigation: 2
---

# Basic concepts on Linux Containers
During this course we will use a number of tools that are stored within a [**Linux Container**](https://en.wikipedia.org/wiki/LXC). 

A Linux Container can be seen as a minimal virtual environment that can be used in any Linux-compatible machine. Using this container allows researchers to save on time and resources related to installation of tools and libraries and improves the reproducibility of the analysis. 
In order for all of us to use exactly the same tools (and their versions) during the course and when you will do your final project, we made up a [**Docker**](https://www.docker.com/) image from [this Dockerfile](https://github.com/biocorecrg/RNAseq_course_2019/blob/master/Dockerfile) and uploaded it in [DockerHub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/rnaseq2019). 

This image can be downloaded and used on computers running Linux/Mac OS or can be converted into another Linux Container called [**Singularity**](https://www.sylabs.io/docs/), which we will be using in this course. 

The Singularity image is a file that can be accessed by the program Singularity for executing any software installed within that image. This image can be created using the following command:

```{bash}
singularity pull docker://biocorecrg/rnaseq2019:1.5
```

And access it as

```{bash}
export RUN="singularity exec $PWD/rnaseq2019-1.5.simg"
$RUN STAR --help
```

If you are a CRG user with an access to the CRG cluster you can access the image as 

```{bash}
export RUN="singularity exec /db/images/rnaseq2019-1.5.simg"
$RUN STAR --help
```

<br/>

If you are not a CRG user ask your IT to install this Singularity image on your cluster.

<br/>

If you want to run the analysis on your own computer (and if your computer has enough RAM), you have the following options:
1. (difficult, time-consuming and might be impossible to install some software on some OSs) Install all sofware on your own one-by-one; for detail, see [this Dockerfile](https://github.com/biocorecrg/RNAseq_course_2019/blob/master/Dockerfile).
2. (should work on any OS) [Install Virtual Box](https://www.virtualbox.org/wiki/Downloads) and [download the virtual machine for this course](https://public-docs.crg.es/biocore/projects/training/vm/2019/). It contains Docker, Singularity and all programms to run the analysis, but it doesn't contain the Singularity image; you will need to make it:
```{bash}
singularity pull docker://biocorecrg/rnaseq2019:1.5
```

3. (should work on any OS) [Install Docker](https://docs.docker.com/v17.12/install/) (you will have to register at Dockerhub) and [install Singularity](https://www.sylabs.io/guides/3.0/user-guide/installation.html) on your computer and make the Singularity image: 
```{bash}
singularity pull docker://biocorecrg/rnaseq2019:1.5
```

<br/>

