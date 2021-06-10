FROM ubuntu:20.04

LABEL description="Bioinformatics Docker Container"
LABEL maintainer="amoustafa@aucegypt.edu"

RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

##########################################################################################
##########################################################################################

RUN apt-get update && \
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
