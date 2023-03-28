FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
MAINTAINER Gabriel Barrett <eaglesnatcher123@gmail.com>

# Build applications
RUN apt-get update && apt-get install --yes --no-install-recommends \
    nano \
    apt-utils \
    liblapack-dev \
    libblas3 \
    apt-transport-https \
    gnupg2 \
    wget \
    build-essential \
    perl \
    python2-dev \
    python-numpy \
    Cython \
    bcftools \
    gfortran \
    plink1.9 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src

RUN wget http://www1.montpellier.inra.fr/CBGP/software/baypass/files/baypass_2.3.tar.gz && tar xvzf baypass_2.3.tar.gz && cd baypass_2.3/sources && make clean all FC=gfortran && make clean && chmod +x g_baypass && rm -R ../../baypass_2.3.tar.gz

ENV PATH=${PATH}:/usr/src/:/usr/src/baypass_2.3/sources/

# Install R
RUN apt-get update -qq && apt-get install --yes --no-install-recommends software-properties-common dirmngr &&\
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc &&\
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" &&\
    apt-get install --yes --no-install-recommends r-base  

# Packages are heavy
RUN R -e 'install.packages("dplyr",repos="http://cloud.r-project.org/",version="1.0.9");install.packages("vcfR",repos="http://cloud.r-project.org/",version="1.12.0");install.packages("ggplot2",repos="http://cloud.r-project.org/",version="3.3.6");install.packages("stringr",repos="http://cloud.r-project.org/",version="1.4.0");install.packages("tidyr",repos="http://cloud.r-project.org/",version="1.2.0")'

#install.packages("BiocManager",repos="http://cloud.r-project.org/",version="3.15");BiocManager::install("LEA",version="3.2.0")
#RUN R -e 'library("BiocManager"); BiocManager::install("LEA",version="3.2.0")'

RUN wget http://bioconductor.org/packages/release/bioc/src/contrib/LEA_3.8.0.tar.gz && R CMD INSTALL LEA_3.8.0.tar.gz

CMD ["bash"]
