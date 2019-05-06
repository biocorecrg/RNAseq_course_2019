---
layout: page
title: Linux Containers
navigation: 2
---

# Basic concepts on Linux Containers
During this course we will use a number of tools that have been stored within a [**Linux Container**](https://en.wikipedia.org/wiki/LXC). 

A Linux Container can be seen as minimal virtual environment that can be used in any Linux compatible machine, this allows us to save a lot of time / resources related to installation of tools and libraries and improves the reproducbility of an analysis. In brief you will nedd to install (or ask your IT to have installed) just one program.

In particular we created a [**Docker**](https://www.docker.com/) image and stored it in [DockerHub](https://cloud.docker.com/u/biocorecrg/repository/docker/biocorecrg/rnaseq2019). 

This image can be downloaded and used in computers running Linux or Mac OS or can be converted into a another Linux Container called [**Singularity**](https://www.sylabs.io/docs/) that is the one that will be used in this course. 

The Singularity image is basically a file that can be accessed by the program Singulairty for executing any software is installed within that image. You can create this image to use the tools we will show in this course at your convenience by using the following command:

```{bash}
singularity pull docker://biocorecrg/rnaseq2019:1.5
```

And access it in the following way 

```{bash}
export RUN="singularity exec $PWD/rnaseq2019-1.5.simg"
$RUN STAR --help
```

If you are a CRG user with access to our cluster you can just access it at this path:

```{bash}
export RUN="singularity exec $PATH/rnaseq2019-1.5.simg"
$RUN STAR --help
```
