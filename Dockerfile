FROM ubuntu:24.04

# Metadata
LABEL description="Bioinformatics Docker Image"
LABEL maintainer="amoustafa@aucegypt.edu"

# Suppress interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

##########################################################################################
##########################################################################################

RUN apt-get update --fix-missing && \
apt-get -y upgrade && \
apt-get -y install apt-utils dialog software-properties-common && \
add-apt-repository universe && \
add-apt-repository multiverse && \
add-apt-repository restricted

##########################################################################################
##########################################################################################

ARG SETUPDIR=/tmp/bioinformatics-toolbox-setup/
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
libreadline* \
libharfbuzz-dev libfribidi-dev

##########################################################################################
##########################################################################################

# RUN ln -s /usr/bin/python3 /usr/bin/python && \
RUN update-alternatives --set java /usr/lib/jvm/java-17-openjdk-amd64/bin/java && \
update-alternatives --set javac /usr/lib/jvm/java-17-openjdk-amd64/bin/javac && \
update-alternatives --set jar /usr/lib/jvm/java-17-openjdk-amd64/bin/jar

##########################################################################################
##########################################################################################

# Virtual Environments
######################
######################

# Miniconda
###########
RUN cd $SETUPDIR/ && wget -t 0 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda3

# Docker
########
RUN cd $SETUPDIR/ && curl -fsSL https://get.docker.com -o get-docker.sh && sh get-docker.sh

# Nexflow
#########
RUN cd $SETUPDIR/ && wget -qO- https://get.nextflow.io | bash && chmod +x nextflow && mv nextflow /usr/local/bin/

##########################################################################################
##########################################################################################

# Bio Progamming
################
################

# BioPerl
#########
RUN apt-get -y install bioperl

# Biopython & Other Python packages
###################################
RUN pip install --no-cache-dir --break-system-packages -U biopython numpy pandas matplotlib scipy seaborn statsmodels plotly bokeh scikit-learn tensorflow keras torch jupyterlab

# R
###
RUN apt update -qq && apt install -y --no-install-recommends software-properties-common dirmngr && \
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
apt-get -y install --no-install-recommends r-base r-base-dev && \
R -e "install.packages (c('tidyverse', 'tidylog', 'readr', 'dplyr', 'reticulate', \
'doParallel', 'devtools', 'tictoc',  \
'knitr', 'printr', 'rmarkdown', 'shiny', \
'gplots', 'cowplot', 'patchwork', 'plotly', 'ggpubr', 'rbokeh', 'circlize', 'RColorBrewer', \
'formattable', 'reshape2', 'data.table', 'readxl',  \
'gganimate', 'gifski', 'av', 'magick', 'ggvis', 'googleVis', \
'randomForest', 'randomForestExplainer', 'randomForestSRC', 'ggRandomForests', \
'rpart', 'rpart.plot', 'xgboost', 'gbm', 'iml', 'boot', 'pROC', 'pdp', 'vip', 'caret', \
'Seurat', 'pheatmap', 'Rtsne', 'vsn', 'vegan', \
'remotes', 'BiocManager'), ask = FALSE)" && \
R -e "BiocManager::install(c('DESeq2', 'edgeR', 'dada2', 'phyloseq', 'metagenomeSeq', 'biomaRt'), ask = FALSE, update = TRUE)"

##########################################################################################
##########################################################################################

# NCBI Tools
############
############

# Had to reinstall camke for some reason
RUN apt-get -y install --reinstall cmake

RUN mkdir -p $SETUPDIR/ncbi && cd $SETUPDIR/ncbi && \
git clone https://github.com/ncbi/ncbi-vdb.git && \
git clone https://github.com/ncbi/ngs.git && \
git clone https://github.com/ncbi/ngs-tools.git && \
git clone https://github.com/ncbi/sra-tools.git

RUN cd $SETUPDIR/ncbi/ncbi-vdb && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-java && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-sdk && ./configure && make && make install && \
cd $SETUPDIR/ncbi/ngs/ngs-python && ./configure && make && make install && \
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
wget -t 0 https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz && \
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
RUN pip install --no-cache-dir --break-system-packages -U seqmagick

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
wget -t 0 https://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz && \
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
wget -t 0 https://github.com/rcedgar/muscle/archive/refs/tags/v5.2.tar.gz && \
tar xvf v5.2.tar.gz && \
cd $SETUPDIR/muscle-5.2/src && \
chmod a+x ./build_linux.bash && ./build_linux.bash && mv ../bin/muscle /usr/local/bin/

# MAFFT
#######
RUN cd $SETUPDIR/ && \
wget -t 0 https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz && \
tar zxvf mafft-7.525-with-extensions-src.tgz && \
cd $SETUPDIR/mafft-7.525-with-extensions/core && \
make clean && make && make install && \
cd $SETUPDIR/mafft-7.525-with-extensions/extensions/ && \
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
wget -t 0 https://downloads.sourceforge.net/project/bbmap/BBMap_39.09.tar.gz && \
tar zxvf BBMap_39.09.tar.gz && \
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
RUN pip install --break-system-packages deeptools

# BEDOPS
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/bedops/bedops.git && \
cd $SETUPDIR/bedops && \
make ; make install ; mv ./bin/* /usr/local/bin/

# SAMBAMBA
##########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz && \
gzip -d sambamba-1.0.1-linux-amd64-static.gz && \
mv sambamba-1.0.1-linux-amd64-static /usr/local/bin/sambamba && \
chmod +x /usr/local/bin/sambamba

##########################################################################################
##########################################################################################

# Assemblers
############
############

# SPAdes
########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz  && \
tar zxvf SPAdes-4.0.0-Linux.tar.gz  && \
mv SPAdes-4.0.0-Linux/bin/* /usr/local/bin/  && \
mv SPAdes-4.0.0-Linux/share/* /usr/local/share/

# ABySS
#######
# RUN cd $SETUPDIR/ && git clone https://github.com/sparsehash/sparsehash.git && cd $SETUPDIR/sparsehash && ./autogen.sh && ./configure && make && make install
RUN /usr/local/miniconda3/bin/conda install -c bioconda -c conda-forge btllib
RUN /usr/local/miniconda3/bin/conda install -c bioconda -c conda-forge abyss

# Velvet
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/dzerbino/velvet.git && \
cd $SETUPDIR/velvet/ && \
make && mv velvet* /usr/local/bin/

# MEGAHIT
#########
# RUN cd $SETUPDIR/ && \
# git clone https://github.com/voutcn/megahit.git && \
# cd $SETUPDIR/megahit && \
# git submodule update --init && \
# mkdir build && \
# cd $SETUPDIR/megahit/build && \
# cmake .. -DCMAKE_BUILD_TYPE=Release && make -j4 && make simple_test  && make install

RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz && \
tar zxvf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz && \
cd $SETUPDIR/MEGAHIT-1.2.9-Linux-x86_64-static/bin/ && \
mv * /usr/local/bin/ && \
cd $SETUPDIR/MEGAHIT-1.2.9-Linux-x86_64-static/share/ && \
mv * /usr/local/share/

# MetaVelvet
############
RUN cd $SETUPDIR/ && \
git clone https://github.com/hacchy/MetaVelvet.git && \
cd $SETUPDIR/MetaVelvet && \
make && mv meta-velvetg /usr/local/bin/

##########################################################################################
##########################################################################################

# Phylogenetics
###############
###############

# TreeTime
##########
RUN pip install --break-system-packages phylo-treetime

# FastTree
##########
RUN cd $SETUPDIR/ && \
wget -t 0 http://www.microbesonline.org/fasttree/FastTree.c && \
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm && \
gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm && \
mv FastTree /usr/local/bin && \
mv FastTreeMP /usr/local/bin

# RAxML
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/stamatak/standard-RAxML.git && \
cd $SETUPDIR/standard-RAxML  && \
rm -fr *.o && make -f Makefile.gcc && cp raxmlHPC /usr/local/bin/  && \
rm -fr *.o && make -f Makefile.SSE3.gcc && cp raxmlHPC-SSE3 /usr/local/bin/  && \
rm -fr *.o && make -f Makefile.PTHREADS.gcc && cp raxmlHPC-PTHREADS /usr/local/bin/  && \
rm -fr *.o && make -f Makefile.SSE3.PTHREADS.gcc && cp raxmlHPC-PTHREADS-SSE3 /usr/local/bin/  && \
rm -fr *.o && make -f Makefile.MPI.gcc && cp raxmlHPC-MPI /usr/local/bin/  && \
rm -fr *.o && make -f Makefile.SSE3.MPI.gcc && cp raxmlHPC-MPI-SSE3 /usr/local/bin/

# RAxML NG
##########
RUN cd $SETUPDIR/ && \
git clone --recursive https://github.com/amkozlov/raxml-ng && \
cd $SETUPDIR/raxml-ng && \
mkdir build && \
cd $SETUPDIR/raxml-ng/build && \
cmake .. && make && mv ../bin/raxml-ng /usr/local/bin/  && \
cmake -DSTATIC_BUILD=ON -DENABLE_RAXML_SIMD=OFF -DENABLE_PLLMOD_SIMD=OFF .. && make && mv ../bin/raxml-ng-static /usr/local/bin/  && \
cmake -DUSE_MPI=ON .. && make && mv ../bin/raxml-ng-mpi /usr/local/bin/

# PhyML
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/stephaneguindon/phyml.git && \
cd $SETUPDIR/phyml/ && \
sh ./autogen.sh && ./configure && make && make install

# Pplacer
#########
# RUN cd $SETUPDIR/ && \
# wget -t 0 https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && \
# unzip pplacer-linux-v1.1.alpha19.zip && \
# cd $SETUPDIR/pplacer-Linux-v1.1.alpha19/ && \
# mv guppy /usr/local/bin/ && \
# mv pplacer /usr/local/bin/ && \
# mv rppr /usr/local/bin/ && \
# cd $SETUPDIR/pplacer-Linux-v1.1.alpha19/scripts/ && \
# python setup.py install

##########################################################################################
##########################################################################################

# Gene Prediction
#################
#################

RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
mv prodigal.linux /usr/local/bin/prodigal && \
chmod +x /usr/local/bin/prodigal

# GlimmerHMM
############
RUN cd $SETUPDIR/ && \
wget -t 0 http://ccb.jhu.edu/software/glimmerhmm/dl/GlimmerHMM-3.0.4.tar.gz && \
tar zxvf GlimmerHMM-3.0.4.tar.gz && \
mv $SETUPDIR/GlimmerHMM/bin/glimmerhmm_linux_x86_64 /usr/local/bin/glimmerhmm && \
chmod +x /usr/local/bin/glimmerhmm

# Infernal
##########
RUN cd $SETUPDIR/ && \
wget -t 0 http://eddylab.org/infernal/infernal-1.1.5.tar.gz && \
tar zxvf infernal-1.1.5.tar.gz && \
cd $SETUPDIR/infernal-1.1.5/ && \
./configure && make && make install

# GECCO
#######
RUN pip install --break-system-packages gecco-tool

# DeepBGC
#########
# Broken because of BioPython removed Bio.Alphabet
# RUN pip install --break-system-packages kiwisolver --force
# RUN pip install --break-system-packages deepbgc
# RUN pip install --break-system-packages deepbgc[hmm]
# RUN deepbgc download

# antiSMASH
###########
# RUN apt-get update && apt-get -y install hmmer2 hmmer fasttree
RUN wget http://dl.secondarymetabolites.org/antismash-stretch.list -O /etc/apt/sources.list.d/antismash.list && \
wget -q -O- http://dl.secondarymetabolites.org/antismash.asc | apt-key add -
RUN cd $SETUPDIR/ && \
wget https://dl.secondarymetabolites.org/releases/7.1.0/antismash-7.1.0.tar.gz && tar -zxf antismash-7.1.0.tar.gz && \
pip install --break-system-packages  ./antismash-7.1.0 --no-deps
# RUN download-antismash-databases

##########################################################################################
##########################################################################################

# Variant Calling
#################
#################

RUN mkdir -p /apps/ 

# GATK
######
RUN cd /apps/ && \
git clone https://github.com/broadinstitute/gatk.git && \
cd /apps/gatk && \
./gradlew

# IGV
#####
# RUN cd /apps/ && \
# wget -t 0 https://data.broadinstitute.org/igv/projects/downloads/snapshot/IGV_Linux_snapshot_WithJava.zip && \
# unzip IGV_Linux_snapshot_WithJava.zip && \
# mv IGV_Linux_snapshot IGV
# RUN rm /apps/IGV_Linux_snapshot_WithJava.zip

# VEP
#####
# RUN cd /apps/ && \
# apt-get -y install cpanminus libtry-tiny-perl libperl4-corelibs-perl && \
# cpanm autodie && \
# cpanm Module::Build && \
# cpanm Test::Warnings && \
# cpanm Bio::DB::HTS::Tabix && \
# git clone https://github.com/Ensembl/ensembl-vep.git && \
# cd ensembl-vep && \
# perl INSTALL.pl --NO_HTSLIB --AUTO alcfp --SPECIES homo_sapiens --ASSEMBLY GRCh38 --PLUGINS all


##########################################################################################
##########################################################################################

# Microbiome
############
############

# Centrifuge
############

RUN cd $SETUPDIR/ && \
git clone https://github.com/DaehwanKimLab/centrifuge.git && \
cd $SETUPDIR/centrifuge/ && \
make && make install

# Pavian
########
RUN R -e  "remotes::install_github('fbreitwieser/pavian')"

# Kraken2
#########
RUN cd $SETUPDIR/ && \
git clone https://github.com/DerrickWood/kraken2.git && \
cd $SETUPDIR/kraken2/src/ && \
make && \
cd $SETUPDIR/kraken2/ && \
bash install_kraken2.sh /usr/local/bin/

RUN cd $SETUPDIR/ && \
git clone https://github.com/jenniferlu717/Bracken.git && \
cd $SETUPDIR/Bracken && \
bash install_bracken.sh && \
mv bracken /usr/local/bin/ && \
mv bracken-build /usr/local/bin/

# MetaPhlAn
###########
RUN pip install --break-system-packages metaphlan panphlan humann
# RUN pip install --break-system-packages graphlan
# RUN metaphlan --install

# mothur
########
# RUN apt-get -y install libreadline*
RUN cd $SETUPDIR/ && \
wget -t 0 http://drive5.com/uchime/uchime4.2.40_i86linux32 && \
mv uchime4.2.40_i86linux32 /usr/local/bin/uchime && \
chmod a+x /usr/local/bin/uchime

RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/mothur/mothur/releases/download/v1.48.1/Mothur.linux_8.zip && \
unzip Mothur.linux_8.zip && \
mv $SETUPDIR/mothur/mothur /usr/local/bin/

# QIIME2
########

## QIIME 2 Amplicon Distribution
RUN /usr/local/miniconda3/bin/conda env create -n qiime2-amplicon-2024.5 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.5-py39-linux-conda.yml

## QIIME 2 Metagenome Distribution¶
RUN /usr/local/miniconda3/bin/conda env create -n qiime2-metagenome-2024.5 --file https://data.qiime2.org/distro/metagenome/qiime2-metagenome-2024.5-py39-linux-conda.yml

##########################################################################################
##########################################################################################

# Finishing
###########
###########
# Removing /usr/local/lib/libgomp.so.1 (seems to be broken) and use /usr/lib/x86_64-linux-gnu/libgomp.so.1 instead
# RUN rm -fr /usr/local/lib/libgomp.so.1

# Creating a soft link because mothur looks for libreadline.so.7
RUN ln -s /usr/lib/x86_64-linux-gnu/libreadline.so.8 /usr/lib/x86_64-linux-gnu/libreadline.so.7

RUN cd $SETUPDIR/ && \
    echo "#!/usr/bin/bash" > $SETUPDIR/init.sh && \
    echo "export PATH=$PATH:/usr/local/ncbi/sra-tools/bin/:/usr/local/ncbi/ngs-tools/bin/:/usr/local/ncbi/ncbi-vdb/bin:/usr/local/miniconda3/bin:/apps/gatk:/apps/IGV:/apps/ensembl-vep:" >> $SETUPDIR/init.sh && \
    echo "source /etc/profile.d/*" >> $SETUPDIR/init.sh && \
    echo "echo '****************************************'" >> $SETUPDIR/init.sh && \
    echo "echo 'Welcome to Bioinformatics Toolbox (v1.4)'" >> $SETUPDIR/init.sh && \
    echo "echo '****************************************'" >> $SETUPDIR/init.sh && \
    echo "echo 'Bioinformatics Toolbox is a docker container for bioinformatics'" >> $SETUPDIR/init.sh && \
    echo "echo " >> $SETUPDIR/init.sh && \
    echo "echo 'For a list of installed tools, please visit: '" >> $SETUPDIR/init.sh && \
    echo "echo 'https://github.com/ahmedmoustafa/bioinformatics-toolbox/blob/master/Tools.md'" >> $SETUPDIR/init.sh && \
    echo "echo " >> $SETUPDIR/init.sh && \
    echo "echo 'If you would like to request adding certain tools or report a problem,'" >> $SETUPDIR/init.sh && \
    echo "echo 'please submit an issue https://github.com/ahmedmoustafa/bioinformatics-toolbox/issues'" >> $SETUPDIR/init.sh && \
    echo "echo " >> $SETUPDIR/init.sh && \
    echo "echo 'If you use Bioinformatics Toolbox in your work, please cite: '" >> $SETUPDIR/init.sh && \
    echo "echo '10.5281/zenodo.8103969'"  >> $SETUPDIR/init.sh && \
    echo "echo 'Have fun!'" >> $SETUPDIR/init.sh && \
    echo "echo ''" >> $SETUPDIR/init.sh && \
    echo "echo ''" >> $SETUPDIR/init.sh && \
    echo "/usr/bin/bash" >> $SETUPDIR/init.sh && \
    mv $SETUPDIR/init.sh /etc/bioinformatics-toolbox.sh && \
    chmod a+x /etc/bioinformatics-toolbox.sh


WORKDIR /root/
ENTRYPOINT ["/etc/bioinformatics-toolbox.sh"]
RUN rm -fr $SETUPDIR
RUN rm -rf /var/lib/apt/lists/*

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
FastTree ; \
treetime --version ; \
raxmlHPC -v ; \
raxml-ng --version ; \
# pplacer --version ; \
samtools  --version ; \
bcftools  --version ; \
bamtools --version ; \
vcftools --version ; \
bedtools --version ; \
deeptools --version ; \
bedops --version ; \
spades.py --version ; \
megahit --version ; \
ABYSS --version ; \
seqkit version ; \
fastp --version ; \
fqtrim -V ; \
seqmagick --version ; \
gecco --version ; \
# deepbgc info ; \
/apps/gatk/gatk --list ; \
# /apps/IGV/igv.sh --version ; \
# vep ; \
centrifuge --version ; \
kraken2 --version ; \
bracken -v ; \
metaphlan --version ; \
humann --version ; \
uchime --version ; \
mothur --version ; \
/usr/local/miniconda3/bin/conda --version ; \
nextflow -version ; \
docker --version