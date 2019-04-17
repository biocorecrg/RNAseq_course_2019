FROM biocorecrg/centos-perlbrew-pyenv

# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu>
MAINTAINER Luca Cozzuto <lucacozzuto@gmail.com> 

ARG FASTQC_VERSION=0.11.7
ARG JAVA_VERSION=1.8.0
ARG RIBO_VERSION=0.4.3
ARG STAR_VERSION=2.6.0c
ARG QUALIMAP_VERSION=2.2.1
ARG SKEWER_VERSION=0.2.2
ARG MULTIQC_VERSION=1.6dev
ARG SAMTOOLS_VERSION=1.8
ARG BEDTOOLS_VERSION=2.27.1
ARG BOWTIE_VERSION=1.2.2
ARG HTSEQ_VERSION=0.9.1
ARG SHORTSTACK_VERSION=3.8.5
ARG TOOL_MULTIQC_VERSION=1.1
ARG GFFREAD_VERSION=0.9.12

#INSTALLING FASTQC
RUN bash -c 'curl -k -L https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip > fastqc.zip'
RUN unzip fastqc.zip; chmod 775 FastQC/fastqc; ln -s $PWD/FastQC/fastqc /usr/local/bin/fastqc

RUN yum install -y java-${JAVA_VERSION}-openjdk which
RUN yum install -y java-${JAVA_VERSION}-openjdk which
RUN yum install -y gzip gunzip

#Installing BOWTIE
RUN bash -c 'curl -k -L https://sourceforge.net/projects/bowtie-bio/files/bowtie/${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip/download > bowtie.zip'
RUN unzip bowtie.zip
RUN mv bowtie-${BOWTIE_VERSION}-linux-x86_64/bowtie* /usr/local/bin/
RUN rm bowtie.zip
RUN rm -fr bowtie-${BOWTIE_VERSION}-linux-x86_64

#Installing SHORTSTACK
RUN bash -c 'curl -k -L https://github.com/MikeAxtell/ShortStack/archive/v${SHORTSTACK_VERSION}.tar.gz > shortstack.tar.gz'
RUN tar -zvxf shortstack.tar.gz 
RUN mv ShortStack-${SHORTSTACK_VERSION}/ShortStack /usr/local/bin/
RUN rm shortstack.tar.gz 
RUN rm -fr ShortStack-${SHORTSTACK_VERSION}

# Installing RiboPicker
RUN bash -c 'curl -k -L https://sourceforge.net/projects/ribopicker/files/standalone/ribopicker-standalone-${RIBO_VERSION}.tar.gz > ribopicker.tar.gz'
RUN tar -zvxf ribopicker.tar.gz
RUN rm ribopicker.tar.gz

#Customization for improving ribopicker
COPY bin/ribopicker_par.pl ./ribopicker-standalone-${RIBO_VERSION}/ribopicker_par.pl
COPY bin/riboPickerConfig.pm ./ribopicker-standalone-${RIBO_VERSION}/
COPY bin/check_rRNA_contam.pl  /usr/local/bin/
COPY bin/ribo.sh  /usr/local/bin/
COPY bin/my_functions.pl  /usr/local/bin/
COPY bin/estim_read_size.sh   /usr/local/bin/
COPY bin/splitfastq.sh /usr/local/bin/

RUN chmod +x /usr/local/bin/check_rRNA_contam.pl
RUN chmod +x /usr/local/bin/ribo.sh
RUN chmod +x /usr/local/bin/splitfastq.sh

RUN mkdir /ribodb
COPY bin/hum_ribo.fa /ribodb

RUN ln -s $PWD/ribopicker-standalone-${RIBO_VERSION}/ribopicker_par.pl   /usr/local/bin/
RUN ln -s $PWD/ribopicker-standalone-${RIBO_VERSION}/riboPickerConfig.pm   /usr/local/bin/
RUN ln -s $PWD/ribopicker-standalone-${RIBO_VERSION}/bwa64   /usr/local/bin/


RUN bwa64 index /ribodb/hum_ribo.fa 

#installing perl modules
RUN cpanm Data::Dumper File::Basename Pod::Usage Getopt::Long Math::Round File::Temp Cwd Pod::Usage FindBin

# Installing RiboPicker
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
# RUN pip install -Iv https://github.com/ewels/MultiQC/archive/v${MULTIQC_VERSION}.tar.gz 
RUN pip install --upgrade --force-reinstall git+https://github.com/ewels/MultiQC.git

# Improving indexing of star
COPY bin/index_star.sh   /usr/local/bin/
RUN chmod +x /usr/local/bin/index_star.sh

# Adding coverage calc 
COPY bin/coverage.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/coverage.sh

# Adding a converter for mirBASE GFF3
COPY bin/fix_smallRNAgff3.pl /usr/local/bin/
RUN chmod +x /usr/local/bin/fix_smallRNAgff3.pl

# Installing samtools
RUN yum install -y xz-devel.x86_64
RUN bash -c 'curl -k -L https://downloads.sourceforge.net/project/samtools/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 > samtools.tar.bz2'
RUN tar -jvxf samtools.tar.bz2
RUN cd samtools-${SAMTOOLS_VERSION}; ./configure; make; make install; cd ../
RUN rm samtools.tar.bz2

# Installing bedtools 
RUN bash -c 'curl -k -L https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz > bedtools.tar.gz'
RUN tar -zvxf bedtools.tar.gz 
RUN cd bedtools2; make; cp ./bin/* /usr/local/bin/
RUN rm bedtools.tar.gz 

#Installing kenttools
RUN yum install -y libpng12
RUN bash -c 'curl -k -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig > /usr/local/bin/bedGraphToBigWig'
RUN bash -c 'curl -k -L http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedSort > /usr/local/bin/bedSort'
RUN bash -c 'chmod +x /usr/local/bin/bed*'

#Adding perl script for improving multiQC report
RUN bash -c 'curl -k -L  https://github.com/CRG-CNAG/make_tool_desc_for_multiqc/archive/v${TOOL_MULTIQC_VERSION}.tar.gz > tool_ver.tar.gz'
RUN tar -zvxf tool_ver.tar.gz
RUN mv make_tool_desc_for_multiqc-${TOOL_MULTIQC_VERSION}/make_tool_desc_for_multiqc.pl /usr/local/bin/ 
RUN chmod +x /usr/local/bin/make_tool_desc_for_multiqc.pl
RUN rm -fr make_tool_desc_for_multiqc-* v${TOOL_MULTIQC_VERSION}.tar.gz 

#Adding perl module
RUN cpanm List::MoreUtils

#Installing HTSEQ
RUN pip install numpy matplotlib  
RUN pip install HTSEQ==${HTSEQ_VERSION}

#Installing GFFREAD
# looking at versions does not work... hope to fix it
RUN git clone https://github.com/gpertea/gclib
RUN git clone https://github.com/gpertea/gffread
RUN cd gffread; make; cp gffread /usr/local/bin/; chmod +x /usr/local/bin/gffread

#cleaning
RUN yum clean all
RUN rm -fr *.zip *.gz *.bz2 
RUN rm -fr STAR-* bedtools2  samtools-*
