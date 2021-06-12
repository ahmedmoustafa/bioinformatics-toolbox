FROM ubuntu:20.04

LABEL description="Bioinformatics Docker Container"
LABEL maintainer="amoustafa@aucegypt.edu"

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

##########################################################################################
##########################################################################################

RUN apt-get update apt-get --fix-missing && \
apt-get -y upgrade && \
apt-get -y install apt-utils dialog software-properties-common

RUN add-apt-repository universe && \
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
default-jre default-jdk ant \
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
caffe-cpu \
cargo \
ffmpeg \
libmagick++-dev \
libavfilter-dev \
dos2unix

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
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
apt-get update && \
apt-get -y install r-base r-base-dev && \
R -e "install.packages (c('tidyverse', 'tidylog', 'readr', 'dplyr', 'knitr', 'printr', 'rmarkdown', 'shiny', \
'ggplot2', 'gplots', 'plotly', 'rbokeh', 'circlize', 'RColorBrewer', 'formattable', \
'reshape2', 'data.table', 'readxl', 'devtools', 'cowplot', 'tictoc', 'ggpubr', 'patchwork', 'reticulate', \
'rpart', 'rpart.plot', 'randomForest', 'randomForestExplainer', 'randomForestSRC', 'ggRandomForests', 'xgboost', 'gbm', 'iml', \
'gganimate', 'gifski', 'av', 'magick', 'ggvis', 'googleVis', \
'pheatmap', 'Rtsne', 'vsn', 'vegan', 'BiocManager'))" && \
R -e "BiocManager::install(c('DESeq2', 'edgeR', 'dada2', 'phyloseq', 'metagenomeSeq', 'biomaRt'), ask = FALSE, update = TRUE)"

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
wget -t 0 https://github.com/shenwei356/seqkit/releases/download/v0.16.1/seqkit_linux_amd64.tar.gz && \
tar zxvf seqkit_linux_amd64.tar.gz && \
mv seqkit /usr/local/bin/

# fastp
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/OpenGene/fastp.git && \
cd $SETUPDIR/fastp && \
make && make install

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
RUN apt-get -y install ncbi-blast+ hmmer
RUN cd $SETUPDIR/ && wget -t 0 http://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz && tar zxvf diamond-linux64.tar.gz && mv diamond /usr/local/bin/

# CD-HIT
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/weizhongli/cdhit.git && \
cd $SETUPDIR/cdhit && \
make && \
make install

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
wget -t 0 https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz && \
tar zxvf muscle3.8.31_src.tar.gz && \
cd $SETUPDIR/muscle3.8.31/src && \
make && mv muscle /usr/local/bin/

# MAFFT
#######
RUN cd $SETUPDIR/ && \
wget -t 0 https://mafft.cbrc.jp/alignment/software/mafft-7.481-with-extensions-src.tgz && \
tar zxvf mafft-7.481-with-extensions-src.tgz && \
cd $SETUPDIR/mafft-7.481-with-extensions/core && \
make clean && make && make install && \
cd $SETUPDIR/mafft-7.481-with-extensions/extensions/ && \
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
wget -t 0 https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz && \
tar zxvf salmon-1.4.0_linux_x86_64.tar.gz && \
mv $SETUPDIR/salmon-latest_linux_x86_64/bin/* /usr/local/bin/ && \
mv $SETUPDIR/salmon-latest_linux_x86_64/lib/* /usr/local/lib/

# kallisto
##########
RUN cd $SETUPDIR/ && \
git clone https://github.com/pachterlab/kallisto.git && \
cd $SETUPDIR/kallisto/ext/htslib && \
autoheader && autoconf && \
cd $SETUPDIR/kallisto/ && \
mkdir build && \
cd $SETUPDIR/kallisto/build && \
cmake .. && make && make install && \
R -e "BiocManager::install('pachterlab/sleuth', ask = FALSE, update = TRUE)"

# BBMap
#######
RUN cd $SETUPDIR/ && \
wget -t 0 https://downloads.sourceforge.net/project/bbmap/BBMap_38.90.tar.gz && \
tar zxvf BBMap_38.90.tar.gz && \
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
git clone git://github.com/samtools/samtools.git && \
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
git clone git://github.com/pezmaster31/bamtools.git && \
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
wget -t 0 https://github.com/biod/sambamba/releases/download/v0.8.0/sambamba-0.8.0-linux-amd64-static.gz && \
gzip -d sambamba-0.8.0-linux-amd64-static.gz && \
mv sambamba-0.8.0-linux-amd64-static /usr/local/bin/sambamba && \
chmod +x /usr/local/bin/sambamba


##########################################################################################
##########################################################################################

# Assemblers
############
############

# SPAdes
########
RUN cd $SETUPDIR/ && \
wget -t 0 http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Linux.tar.gz  && \
tar zxvf SPAdes-3.15.2-Linux.tar.gz  && \
mv SPAdes-3.15.2-Linux/bin/* /usr/local/bin/  && \
mv SPAdes-3.15.2-Linux/share/* /usr/local/share/


# ABySS
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/sparsehash/sparsehash.git && \
cd $SETUPDIR/sparsehash && \
./autogen.sh && ./configure && make && make install && \
cd $SETUPDIR/ && \
git clone https://github.com/bcgsc/abyss.git && \
cd $SETUPDIR/abyss && \
./autogen.sh && ./configure && make && make install


# Velvet
########
RUN cd $SETUPDIR/ && \
git clone https://github.com/dzerbino/velvet.git && \
cd $SETUPDIR/velvet/ && \
make && mv velvet* /usr/local/bin/


# MEGAHIT
#########
RUN cd $SETUPDIR/ && \
git clone https://github.com/voutcn/megahit.git && \
cd $SETUPDIR/megahit && \
git submodule update --init && \
mkdir build && \
cd $SETUPDIR/megahit/build && \
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j4 && make simple_test  && make install


# MetaVelvet
############
RUN cd $SETUPDIR/ && \
git clone git://github.com/hacchy/MetaVelvet.git && \
cd $SETUPDIR/MetaVelvet && \
make && mv meta-velvetg /usr/local/bin/

##########################################################################################
##########################################################################################

# Phylogenetics
###############
###############

# TreeTime
##########
RUN pip3 install phylo-treetime


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
cmake -DUSE_MPI=ON .. && make && mv ../bin/raxml-ng /usr/local/bin/raxml-ng-mpi


# PhyML
#######
RUN cd $SETUPDIR/ && \
git clone https://github.com/stephaneguindon/phyml.git && \
cd $SETUPDIR/phyml/ && \
sh ./autogen.sh && ./configure && make && make install


# Pplacer
#########
RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && \
unzip pplacer-linux-v1.1.alpha19.zip && \
cd $SETUPDIR/pplacer-Linux-v1.1.alpha19/ && \
mv guppy /usr/local/bin/ && \
mv pplacer /usr/local/bin/ && \
mv rppr /usr/local/bin/ && \
cd $SETUPDIR/pplacer-Linux-v1.1.alpha19/scripts/ && \
python setup.py install

##########################################################################################
##########################################################################################

# Gene Prediction
#################
#################

RUN cd $SETUPDIR/ && \
wget -t 0 https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
mv prodigal.linux /usr/local/bin/prodigal && \
chmod +x /usr/local/bin/prodigal

# Infernal
##########
RUN cd $SETUPDIR/ && \
wget -t 0 http://eddylab.org/infernal/infernal-1.1.4.tar.gz && \
tar zxvf infernal-1.1.4.tar.gz && \
cd $SETUPDIR/infernal-1.1.4/ && \
./configure && make && make install

##########################################################################################
##########################################################################################

# Misc
######
######

# Docker
########
RUN cd $SETUPDIR/ && \
wget -t 0 https://get.docker.com/ -O docker.sh && \
sh docker.sh

# Miniconda
###########
RUN cd $SETUPDIR/ && \
wget -t 0 https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
sh Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda3

##########################################################################################
##########################################################################################

# Finishing
###########
###########
# Removing /usr/local/lib/libgomp.so.1 (seems to be broken) and use /usr/lib/x86_64-linux-gnu/libgomp.so.1 instead
RUN rm -fr /usr/local/lib/libgomp.so.1

RUN cd $SETUPDIR/
RUN echo "#!/usr/bin/bash" > $SETUPDIR/init.sh
RUN echo "export PATH=$PATH:/usr/local/ncbi/sra-tools/bin/:/usr/local/ncbi/ngs-tools/bin/:/usr/local/ncbi/ncbi-vdb/bin:/usr/local/miniconda3/bin" >> $SETUPDIR/init.sh
RUN echo "source /etc/profile.d/*" >> $SETUPDIR/init.sh
RUN echo "echo '-----------------'" >> $SETUPDIR/init.sh
RUN echo "echo 'Welcome to Bioinformatics Toolbox (v1.0)'" >> $SETUPDIR/init.sh
RUN echo "echo '----------------------------------------'" >> $SETUPDIR/init.sh
RUN echo "echo 'Bioinformatics Toolbox is a docker container for bioinformatics'" >> $SETUPDIR/init.sh
RUN echo "echo 'Maintained by Ahmed Moustafa (amoustafa@aucegypt.edu)'" >> $SETUPDIR/init.sh
RUN echo "echo 'For a list of installed tools, please visit: '" >> $SETUPDIR/init.sh
RUN echo "echo 'https://github.com/ahmedmoustafa/bioinformatics-toolbox/blob/master/Tools.md'" >> $SETUPDIR/init.sh
RUN echo "echo 'If you use Bioinformatics Toolbox in your work, please cite: '" >> $SETUPDIR/init.sh
RUN echo "echo '10.5281/zenodo.3723585'"  >> $SETUPDIR/init.sh
RUN echo "echo 'Have fun!'" >> $SETUPDIR/init.sh
RUN echo "echo ''" >> $SETUPDIR/init.sh
RUN echo "echo ''" >> $SETUPDIR/init.sh
RUN echo "/usr/bin/bash" >> $SETUPDIR/init.sh
RUN echo "" >> $SETUPDIR/init.sh
RUN mv $SETUPDIR/init.sh /etc/bioinformatics-toolbox.sh
RUN chmod a+x /etc/bioinformatics-toolbox.sh

WORKDIR /root/
ENTRYPOINT ["/etc/bioinformatics-toolbox.sh"]
RUN rm -fr $SETUPDIR

# Versions
##########
RUN python --version ; \
R --version ; \
blastn -version ; \
diamond --version ; \
muscle -version ; \
mafft --version ; \
# tophat --version ; \
hisat2 --version ; \
bowtie2 --version ; \
STAR --version ; \
salmon --version ; \
bbmap.sh --version ; \
hts_Stats --version ; \
treetime --version ; \
# RUN FastTree
# RUN phyml --version
raxmlHPC -v ; \
raxml-ng --version ; \
pplacer --version ; \
samtools  --version ; \
bcftools  --version ; \
bamtools --version ; \
vcftools --version ; \
bedtools --version ; \
deeptools --version ; \
bedops --version ; \
spades.py --version ; \
megahit --version ; \
spades.py --version ; \
seqkit version ; \
fastp --version ; \
fqtrim -V ; \
seqmagick --version ; \
docker --version ; \
/usr/local/miniconda3/bin/conda --version

##########################################################################################
##########################################################################################
