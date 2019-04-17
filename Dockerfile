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

#INSTALLING FASTQC
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip > fastqc.zip'
RUN unzip fastqc.zip; chmod 775 FastQC/fastqc; ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc


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



#cleaning
RUN yum clean all
RUN rm -fr *.zip *.gz *.bz2 
RUN rm -fr STAR-* bedtools2  samtools-*
