# This Dockerfile is based on https://github.com/ahmedmoustafa/bioinformatics-toolbox
FROM ubuntu:22.04

LABEL description="ECRRM Genomics Commons Docker Image"
LABEL maintainer="amoustafa@aucegypt.edu"
LABEL version="1.0"


RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

##########################################################################################
##########################################################################################

RUN apt-get update --fix-missing && \
apt-get -y upgrade && \
apt-get -y install apt-utils dialog software-properties-common
RUN add-apt-repository universe && \
add-apt-repository multiverse && \
add-apt-repository restricted

##########################################################################################
##########################################################################################

ARG SETUPDIR=/tmp/genomics-commons-setup/
RUN mkdir -p $SETUPDIR
WORKDIR $SETUPDIR

##########################################################################################
##########################################################################################

# Prerequisites
###############
###############

RUN apt-get -y install vim nano emacs rsync curl wget screen htop parallel gnupg lsof git locate unrar bc aptitude unzip bison flex \
build-essential libtool autotools-dev automake autoconf cmake \
libboost-dev libboost-all-dev libboost-system-dev libboost-program-options-dev libboost-iostreams-dev libboost-filesystem-dev \
gfortran libgfortran5 \
openjdk-17* ant \
python3 python3-dev python3-pip python3-venv \
libssl-dev libcurl4-openssl-dev \
libxml2-dev \
libmagic-dev \
hdf5-* libhdf5-* \
fuse libfuse-dev \
libtbb-dev \
liblzma-dev libbz2-dev \
libbison-dev \
libgmp3-dev \
libncurses5-dev libncursesw5-dev \
liblzma-dev \
# caffe-cpu \
cargo \
ffmpeg \
libmagick++-dev \
libavfilter-dev \
dos2unix \
git-lfs \
apt-transport-https \
autopoint po4a doxygen \
libreadline-dev

##########################################################################################
##########################################################################################

# Progamming
############
############

# BioPerl
#########
RUN apt-get -y install bioperl

# Biopython
###########
RUN pip3 install --no-cache-dir -U biopython numpy pandas matplotlib scipy seaborn statsmodels plotly bokeh scikit-learn tensorflow keras torch theano jupyterlab

# R
###
RUN apt-get update -qq && apt-get install --no-install-recommends software-properties-common dirmngr && \
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get -y install --no-install-recommends r-base r-base-dev

RUN R -e "install.packages (c('remotes, tidyverse', 'tidylog', 'readr', 'dplyr', 'knitr', 'printr', 'rmarkdown', 'shiny', \
'ggplot2', 'gplots', 'plotly', 'rbokeh', 'circlize', 'RColorBrewer', 'formattable', \
'reshape2', 'data.table', 'readxl', 'devtools', 'cowplot', 'tictoc', 'ggpubr', 'patchwork', 'reticulate', \
'rpart', 'rpart.plot', 'randomForest', 'randomForestExplainer', 'randomForestSRC', 'ggRandomForests', 'xgboost', 'gbm', 'iml', \
'gganimate', 'gifski', 'av', 'magick', 'ggvis', 'googleVis', \
'pheatmap', 'Rtsne', 'vsn', 'vegan', 'BiocManager'), ask = FALSE)"

RUN R -e "BiocManager::install(c('DESeq2', 'edgeR', 'dada2', 'phyloseq', 'metagenomeSeq', 'biomaRt'), ask = FALSE, update = TRUE)"

##########################################################################################
##########################################################################################

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

##########################################################################################
##########################################################################################

# NCBI Tools
############
############

RUN mkdir -p $SETUPDIR/ncbi && cd $SETUPDIR/ncbi && \
git clone https://github.com/ncbi/ncbi-vdb.git && \
git clone https://github.com/ncbi/ngs.git && \
git clone https://github.com/ncbi/ngs-tools.git && \
git clone https://github.com/ncbi/sra-tools.git && \
cd $SETUPDIR/ncbi/ncbi-vdb && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-sdk && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-python && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-java && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-bam && ./configure && make && make install && \
cd $SETUPDIR/ncbi/sra-tools && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs-tools && ./configure && make && make install

RUN cd $SETUPDIR/ncbi && \
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' && \
chmod +x datasets && \
mv datasets /usr/local/bin/

##########################################################################################
##########################################################################################

# Sequence Processing
#####################
#####################

# FASTX
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/agordon/libgtextutils.git && cd $SETUPDIR/libgtextutils/ && \
./reconf && ./configure && make && make install && \
cd $SETUPDIR/ && \
git clone https://github.com/agordon/fastx_toolkit.git && cd $SETUPDIR/fastx_toolkit && \
wget -t 0 https://github.com/agordon/fastx_toolkit/files/1182724/fastx-toolkit-gcc7-patch.txt && \
patch -p1 < fastx-toolkit-gcc7-patch.txt && \
./reconf && ./configure && make && make install

# Trimmomatic
#############
RUN cd $SETUPDIR/ && \
git clone https://github.com/timflutre/trimmomatic.git && \
cd $SETUPDIR/trimmomatic && \
make && make install INSTALL="/usr/local/"

# SeqKit
########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/shenwei356/seqkit/releases/download/v2.4.0/seqkit_linux_amd64.tar.gz && \
tar zxvf seqkit_linux_amd64.tar.gz && \
mv seqkit /usr/local/bin/

# fastp
#######
# RUN cd $SETUPDIR/ && \
# git clone https://github.com/OpenGene/fastp.git && \
# cd $SETUPDIR/fastp && \
# make && make install
RUN wget http://opengene.org/fastp/fastp && \
chmod a+x ./fastp && \
mv ./fastp /usr/local/bin/


# HTStream
##########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/s4hts/HTStream/releases/download/v1.3.3/HTStream_v1.3.3.tar.gz && \
tar zxvf HTStream_v1.3.3.tar.gz && \
mv hts_* /usr/local/bin/

# fqtrim
########
RUN cd $SETUPDIR/ && \
wget -t 0 http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.9.7.tar.gz && \
tar zxvf fqtrim-0.9.7.tar.gz && \
cd $SETUPDIR/fqtrim-0.9.7/ && \
make && mv fqtrim /usr/local/bin/

# seqmagick
###########
RUN pip3 install --no-cache-dir -U seqmagick

# seqtk
#######
RUN git clone https://github.com/lh3/seqtk.git && \
cd seqtk && \
make && \
mv seqtk /usr/local/bin/

# Parallel fastq-dump
#####################
RUN cd $SETUPDIR/ && \
git clone https://github.com/rvalieris/parallel-fastq-dump.git && \
cd $SETUPDIR/parallel-fastq-dump/ && \
mv parallel-fastq-dump /usr/local/bin/

##########################################################################################
##########################################################################################

# Sequence Search
#################
#################

# BLAST & HMMER
###############
RUN apt-get -y install ncbi-blast+ hmmer hmmer2

# Diamond
#########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz && \
tar zxvf diamond-linux64.tar.gz && \
mv diamond /usr/local/bin/

# CD-HIT
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/weizhongli/cdhit.git && \
cd $SETUPDIR/cdhit && \
make && make install

# vsearch
#########
RUN cd $SETUPDIR/ && \
git clone https://github.com/torognes/vsearch.git && \
cd $SETUPDIR/vsearch/ && \
./autogen.sh && ./configure && make && make install

##########################################################################################
##########################################################################################

# Alignment Tools
#################
#################

# JAligner
##########
RUN apt-get -y install jaligner

# MUSCLE
########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/rcedgar/muscle/archive/refs/tags/5.1.0.tar.gz && \
tar xvf 5.1.0.tar.gz && \
cd $SETUPDIR/muscle-5.1.0/src && \
make && mv Linux/muscle /usr/local/bin/

# MAFFT
#######
RUN cd $SETUPDIR/ && \
wget -t 0 https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz && \
tar zxvf mafft-7.505-with-extensions-src.tgz && \
cd $SETUPDIR/mafft-7.505-with-extensions/core && \
make clean && make && make install && \
cd $SETUPDIR/mafft-7.505-with-extensions/extensions/ && \
make clean && make && make install

# BWA
#####
RUN cd $SETUPDIR/ && \
git clone https://github.com/lh3/bwa.git && \
cd $SETUPDIR/bwa && \
make && mv bwa /usr/local/bin/

# TopHat
########
# RUN cd $SETUPDIR/ && \
# wget -t 0 https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
# tar zxvf tophat-2.1.1.Linux_x86_64.tar.gz && \
# cd $SETUPDIR/tophat-2.1.1.Linux_x86_64 && \
# mv tophat* /usr/local/bin/

# HISAT2
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/infphilo/hisat2.git && \
cd $SETUPDIR/hisat2 && \
make && mv hisat2-* /usr/local/bin/  &&  mv hisat2 /usr/local/bin/

# HISAT-3N
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n && \
cd hisat-3n && \
git checkout -b hisat-3n origin/hisat-3n && \
make && mv hisat-3n* /usr/local/bin/

# Bowtie2
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/BenLangmead/bowtie2.git && \
cd $SETUPDIR/bowtie2/ && \
make && make install

# STAR
######
RUN cd $SETUPDIR/ && \
git clone https://github.com/alexdobin/STAR.git && \
cd $SETUPDIR/STAR/source && \
make STAR && mv STAR /usr/local/bin/

# Salmon
########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
tar zxvf salmon-1.10.0_linux_x86_64.tar.gz && \
mv $SETUPDIR/salmon-latest_linux_x86_64/bin/* /usr/local/bin/ && \
mv $SETUPDIR/salmon-latest_linux_x86_64/lib/* /usr/local/lib/

# kallisto
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/pachterlab/kallisto.git && \
cd $SETUPDIR/kallisto/ext/htslib && \
autoheader && autoconf && \
make -j CFLAGS=-D_GNU_SOURCE lib-static && \
cd $SETUPDIR/kallisto/ && \
mkdir build && \
cd $SETUPDIR/kallisto/build && \
cmake .. && make && make install

RUN R -e "BiocManager::install('pachterlab/sleuth', ask = FALSE, update = TRUE)"

# BBMap
#######
RUN cd $SETUPDIR/ && \
wget -t 0 https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && \
tar zxvf BBMap_39.01.tar.gz && \
mv bbmap/* /usr/local/bin/

##########################################################################################
##########################################################################################

# BAM Processing
################
################

# HTSlib
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/samtools/htslib.git && \
cd $SETUPDIR/htslib && \
autoreconf -i ; git submodule update --init --recursive ; ./configure ; make ; make install

# Samtools
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/samtools/samtools.git && \
cd $SETUPDIR/samtools && \
autoheader ; autoconf ; ./configure ; make ; make install

# Bcftools
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/samtools/bcftools.git && \
cd $SETUPDIR/bcftools && \
autoheader ; autoconf ; ./configure ; make ; make install

# Bamtools
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/pezmaster31/bamtools.git && \
cd $SETUPDIR/bamtools && \
mkdir build && \
cd $SETUPDIR/bamtools/build && \
cmake .. ; make ; make install

# VCFtools
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/vcftools/vcftools.git && \
cd $SETUPDIR/vcftools && \
./autogen.sh ; ./configure ; make ; make install

# Bedtools
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/arq5x/bedtools2.git && \
cd $SETUPDIR/bedtools2 && \
make ; make install

# deepTools
###########
RUN cd $SETUPDIR/ && \
git clone https://github.com/deeptools/deepTools && \
cd $SETUPDIR/deepTools && \
python setup.py install

# BEDOPS
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/bedops/bedops.git && \
cd $SETUPDIR/bedops && \
make ; make install ; mv ./bin/* /usr/local/bin/

# SAMBAMBA
##########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/biod/sambamba/releases/download/v0.8.2/sambamba-0.8.2-linux-amd64-static.gz && \
gzip -d sambamba-0.8.2-linux-amd64-static.gz && \
mv sambamba-0.8.2-linux-amd64-static /usr/local/bin/sambamba && \
chmod +x /usr/local/bin/sambamba

##########################################################################################
##########################################################################################

# Misc
######
######

# Docker
########

# RUN cd $SETUPDIR/ && \
# wget -t 0 https://get.docker.com/ -O docker.sh && \
# sh docker.sh

# Miniconda
###########
RUN cd $SETUPDIR/ && \
wget -t 0 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda3

# Nexflow
#########
RUN cd $SETUPDIR/ && \
wget -qO- https://get.nextflow.io | bash && \
chmod +x nextflow && \
mv nextflow /usr/local/bin/

##########################################################################################
##########################################################################################

# Variant Calling
#################
#################

# GATK
######
RUN mkdir -p /apps/ && \
cd /apps/ && \
git clone https://github.com/broadinstitute/gatk.git && \
cd /apps/gatk && \
./gradlew

# IGV
#####
RUN cd /apps/ && \
wget -t 0 https://data.broadinstitute.org/igv/projects/downloads/snapshot/IGV_Linux_snapshot_WithJava.zip && \
unzip IGV_Linux_snapshot_WithJava.zip && \
mv IGV_Linux_snapshot IGV

# VEP
#####
RUN cd /apps/ && \
apt-get -y install cpanminus libtry-tiny-perl libperl4-corelibs-perl && \
cpanm autodie && \
cpanm Module::Build && \
cpanm Bio::DB::HTS::Tabix && \
git clone https://github.com/Ensembl/ensembl-vep.git && \
cd ensembl-vep
# RUN cd /apps/ensembl-vep && perl INSTALL.pl --NO_HTSLIB --AUTO alcfp --SPECIES homo_sapiens --ASSEMBLY GRCh38 --PLUGINS all


##########################################################################################
##########################################################################################

# Population Genetics
#####################
#####################

# PLINK 1.90
############
RUN cd $SETUPDIR/ && \
wget -t 0 https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip && \
unzip plink_linux_x86_64_20230116.zip && \
mv plink /usr/local/bin/ && \
mv prettify /usr/local/bin/ && \


##########################################################################################
##########################################################################################


# Finishing
###########
###########
# Removing /usr/local/lib/libgomp.so.1 (seems to be broken) and use /usr/lib/x86_64-linux-gnu/libgomp.so.1 instead
# RUN rm -fr /usr/local/lib/libgomp.so.1

# Creating a soft link because mothur looks for libreadline.so.7
RUN ln -s /usr/lib/x86_64-linux-gnu/libreadline.so.8 /usr/lib/x86_64-linux-gnu/libreadline.so.7

RUN cd $SETUPDIR/
RUN echo "#!/usr/bin/bash" > $SETUPDIR/init.sh
RUN echo "export PATH=$PATH:/usr/local/ncbi/sra-tools/bin/:/usr/local/ncbi/ngs-tools/bin/:/usr/local/ncbi/ncbi-vdb/bin:/usr/local/miniconda3/bin:/apps/gatk:/apps/IGV:/apps/ensembl-vep:" >> $SETUPDIR/init.sh
RUN echo "source /etc/profile.d/*" >> $SETUPDIR/init.sh
RUN echo "echo '****************************************'" >> $SETUPDIR/init.sh
RUN echo "echo 'Welcome to ECRRM Genomics Commons (v1.0)'" >> $SETUPDIR/init.sh
RUN echo "echo '****************************************'" >> $SETUPDIR/init.sh
RUN echo "echo 'CRRM Genomics Commons is a docker container for bioinformatics'" >> $SETUPDIR/init.sh
RUN echo "echo " >> $SETUPDIR/init.sh
RUN echo "echo 'For a list of installed tools, please visit: '" >> $SETUPDIR/init.sh
RUN echo "echo 'https://github.com/ECRRM/genomics-commons/blob/master/Tools.md'" >> $SETUPDIR/init.sh
RUN echo "echo " >> $SETUPDIR/init.sh
RUN echo "echo 'If you would like to request adding certain tools or report a problem,'" >> $SETUPDIR/init.sh
RUN echo "echo 'please submit an issue https://github.com/ECRRM/genomics-commons/issues'" >> $SETUPDIR/init.sh
RUN echo "echo " >> $SETUPDIR/init.sh
RUN echo "echo 'Have fun!'" >> $SETUPDIR/init.sh
RUN echo "echo ''" >> $SETUPDIR/init.sh
RUN echo "echo ''" >> $SETUPDIR/init.sh
RUN echo "" >> $SETUPDIR/init.sh
RUN mv $SETUPDIR/init.sh /etc/genomics-commons.sh
RUN chmod a+x /etc/genomics-commons.sh

WORKDIR /root/
CMD ["/etc/genomics-commons.sh"]
RUN rm -fr $SETUPDIR

# Versions
##########
RUN python --version ; \
java -version ; \
R --version ; \
blastn -version ; \
diamond --version ; \
vsearch --version ; \
muscle -version ; \
mafft --version ; \
# tophat --version ; \
hisat2 --version ; \
bowtie2 --version ; \
STAR --version ; \
salmon --version ; \
bbmap.sh --version ; \
hts_Stats --version ; \
samtools  --version ; \
bcftools  --version ; \
bamtools --version ; \
vcftools --version ; \
bedtools --version ; \
deeptools --version ; \
bedops --version ; \
seqkit version ; \
fastp --version ; \
fqtrim -V ; \
seqmagick --version ; \
/apps/gatk/gatk --list ; \
/apps/IGV/igv.sh --version ; \
/usr/local/miniconda3/bin/conda --version ; \
nextflow -version ;\
plink --version

##########################################################################################
##########################################################################################
