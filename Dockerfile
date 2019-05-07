FROM biocorecrg/centos-perlbrew-pyenv-java


# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu>
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

ARG FASTQC_VERSION=0.11.8
ARG STAR_VERSION=2.7.0f
ARG QUALIMAP_VERSION=2.2.1
ARG SKEWER_VERSION=0.2.2
ARG MULTIQC_VERSION=1.7
ARG SAMTOOLS_VERSION=1.9
ARG BCFTOOLS_VERSION=1.9
ARG FASTQSCREEN_VERSION=0.13.0
ARG BOWTIE2_VERSION=2.3.5.1
ARG SALMON_VERSION=0.13.1
ARG R_VERSION=3.5.2
ARG SRATOOLKIT_VERSION=2.9.6

#upgrading pip
RUN pip install --upgrade pip

#INSTALLING SRATOOLKIT
RUN bash -c 'curl -k -L https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-centos_linux64.tar.gz > sratoolkit.tar.gz'
RUN tar -zvxf sratoolkit.tar.gz; cd sratoolkit.${SRATOOLKIT_VERSION}-centos_linux64/bin; ln -s $PWD/fastq-dump /usr/local/bin/fastq-dump; cd ../../

#INSTALLING FASTQC
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip > fastqc.zip'
RUN unzip fastqc.zip; chmod 775 FastQC/fastqc; ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

#INSTALLING FASTQ_SCREEN
RUN yum install -y perl-GD wget gd-devel.x86_64
RUN cpanm -f GD::Graph::bars
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v${FASTQSCREEN_VERSION}.tar.gz > fastqc_screen.tar.gz'
RUN tar -zvxf fastqc_screen.tar.gz
RUN cp fastq_screen_v${FASTQSCREEN_VERSION}/* /usr/local/bin/
RUN sed s/perl/env\ perl/g fastq_screen_v${FASTQSCREEN_VERSION}/fastq_screen > /usr/local/bin/fastq_screen 
RUN rm fastqc_screen.tar.gz

#INSTALLING BOWTIE2
RUN yum install -y yum install tbb.x86_64
RUN bash -c 'curl -k -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip/download > bowtie2.zip'
RUN unzip bowtie2.zip; mv bowtie2-${BOWTIE2_VERSION}*/bowtie2* /usr/local/bin/


# Installing Skewer
RUN bash -c 'curl -k -L https://downloads.sourceforge.net/project/skewer/Binaries/skewer-${SKEWER_VERSION}-linux-x86_64 > /usr/local/bin/skewer'
RUN chmod +x /usr/local/bin/skewer

# Installing STAR
RUN bash -c 'curl -k -L https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz > STAR.tar.gz'
RUN tar -zvxf STAR.tar.gz
RUN cp STAR-${STAR_VERSION}/bin/Linux_x86_64/* /usr/local/bin/
RUN rm STAR.tar.gz

# Installing QUALIMAP
RUN bash -c 'curl -k -L https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v${QUALIMAP_VERSION}.zip > qualimap.zip'
RUN unzip qualimap.zip
RUN rm qualimap.zip

RUN ln -s $PWD/qualimap_v${QUALIMAP_VERSION}/qualimap /usr/local/bin/

# Installing MULTIQC // Latest dev version is much better. 
RUN pip install -Iv https://github.com/ewels/MultiQC/archive/v${MULTIQC_VERSION}.tar.gz 

# Installing samtools
RUN yum install -y xz-devel.x86_64
RUN bash -c 'curl -k -L https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 > samtools.tar.bz2' 
RUN tar -jvxf samtools.tar.bz2
RUN cd samtools-${SAMTOOLS_VERSION}; ./configure; make; make install; cd ../ 

# Installing BCFtools
RUN bash -c 'curl -k -L https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 > bcftools.tar.bz2'
RUN tar -jvxf bcftools.tar.bz2
RUN cd bcftools-${BCFTOOLS_VERSION}; ./configure; make; make install; cd ../

# Installing salmon
RUN bash -c 'curl -k -L https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz > salmon.tar.gz'
RUN tar -zvxf salmon.tar.gz
RUN ln -s $PWD/salmon-latest_linux_x86_64/bin/salmon /usr/local/bin/salmon

# Install R and R packages
RUN yum install epel-release libxml2-devel libcurl-devel -y
RUN yum install R-${R_VERSION} -y
RUN mkdir -p /usr/share/doc/R-${R_VERSION}/html
RUN Rscript -e "install.packages('BiocManager', , repos='http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install(c('GenomicRanges', 'SummarizedExperiment', 'genefilter', 'geneplotter'))"
RUN Rscript -e "BiocManager::install(c('DESeq2', 'tximport', 'pheatmap'), dependencies=TRUE)"

#cleaning
RUN yum clean all
RUN rm -fr *.zip *.gz *.bz2 
RUN rm -fr STAR-* bedtools2  samtools-*
