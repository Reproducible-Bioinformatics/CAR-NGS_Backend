FROM library/ubuntu:20.04 as UBUNTU_BASE
MAINTAINER alessandri.luca1991@gmail.com
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

#COPY ./R-4.1.1 /tmp
RUN apt-get update
RUN apt-get -y install gfortran
RUN apt-get -y install build-essential
RUN apt-get -y install fort77
RUN apt-get -y install xorg-dev
RUN apt-get -y install liblzma-dev  libblas-dev gfortran
RUN apt-get -y install gobjc++
RUN apt-get -y install aptitude
RUN apt-get -y install libbz2-dev
RUN apt-get -y install libpcre3-dev
RUN aptitude -y install libreadline-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get -y install libcairo2-dev
RUN apt-get -y install libxt-dev
#COPY pcre2-10.37.tar.gz /tmp/prec/
#RUN tar -zxvf /tmp/prec/pcre2-10.37.tar.gz
#RUN /pcre2-10.37/configure
#RUN make -j 24
#RUN make install
#RUN chmod +x /tmp/configure
#RUN /tmp/configure --with-x=yes --enable-R-shlib
#RUN make
#RUN make install
#COPY p7zip_16.02 /tmp/

RUN apt-get update \
    && apt-get -y install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev wget \
    && apt-get -y install libblas-dev liblapack-dev \
    && apt-get -y install libatlas-base-dev \
    && apt-get -y install gfortran \
    && apt-get -y install libxml2-dev libxmlsec1-dev \
    && apt-get install -y libbz2-dev \
    && apt-get -y install pkg-config\
    && apt-get update && apt-get install -y --no-install-recommends apt-utils \
    && apt-get -y install libhdf5-dev
RUN apt-get update
RUN apt-get install -y liblzma-dev
RUN apt-get update
RUN apt-get install -y sqlite3
RUN apt-get install -y libsqlite3-dev

ENV USER=user
RUN apt-get update && apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
RUN echo \
  "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
RUN apt-get update && apt-get install -y docker-ce docker-ce-cli containerd.io
#COPY ss.sh /ss.sh

RUN apt install sudo
RUN curl -L http://xrl.us/installperlnix | bash
#RUN apt-get install libxslt
RUN apt-get install libxml2
RUN apt-get install libxml2-dev
RUN apt-get install libxslt-dev
RUN yes Y | apt-get update
RUN yes Y | apt-get upgrade
RUN apt-get install build-essential
RUN apt-get install gcc
RUN yes Y | apt install imagemagick
RUN apt-get install bzip2
#RUN apt-get install bzip2-dev
RUN apt install zlib1g
RUN apt install zlib1g-dev
RUN yes Y | apt install --fix-missing build-essential libxml2-dev libssl-dev
RUN yes Y | apt install libharfbuzz-dev libfribidi-dev 
RUN sudo apt update
RUN apt install -y libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN yes Y | apt-get install gawk
WORKDIR /home
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
RUN tar xvzf trim_galore.tar.gz
RUN yes Y | apt install cutadapt
WORKDIR /home
#install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2
RUN tar xvjf samtools-1.3.tar.bz2 --no-same-owner
WORKDIR samtools-1.3
RUN ./configure
RUN make
RUN sudo make install
#ENV PATH="$PATH:/home/samtools-1.3"
WORKDIR /home
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
RUN tar xvjf samtools-1.7.tar.bz2
WORKDIR samtools-1.7
RUN ./configure
RUN make
RUN  make install
#ENV PATH="$PATH:/home/samtools-1.7"
WORKDIR /home

#STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz
RUN tar -xzf 2.5.1b.tar.gz --no-same-owner
WORKDIR STAR-2.5.1b
WORKDIR ./source
RUN make STAR
WORKDIR /home

#install bedtools
RUN apt install python-is-python3
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
RUN tar -zxvf bedtools-2.27.1.tar.gz
WORKDIR bedtools2
RUN make
WORKDIR /home
# to check # RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
# to check # RUN tar -xzf 2.7.10a.tar.gz
# to ckeck# RUN apt-get install -y rsem
#to check # RUN apt install python-is-python3
#to check # RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
#to check # RUN tar -zxvf bedtools-2.27.1.tar.gz
RUN apt-get install -y bedops
#to check # RUN yes Y | apt-get install autoconf
#ENV R_VERSION=4.2.3
RUN sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" | sudo tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | sudo apt-key add -
RUN sudo apt-get update
RUN sudo apt-get install -y r-base r-base-dev r-cran-xml r-cran-rjava libcurl4-openssl-dev
RUN sudo apt-get install -y libssl-dev libxml2-dev openjdk-7-* libgdal-dev libproj-dev libgsl-dev
RUN sudo apt-get install -y xml2 default-jre default-jdk mesa-common-dev libglu1-mesa-dev freeglut3-dev
RUN sudo apt-get install -y mesa-common-dev libx11-dev r-cran-rgl r-cran-rglpk r-cran-rsymphony r-cran-plyr
RUN sudo apt-get install -y  r-cran-reshape  r-cran-reshape2 r-cran-rmysql
RUN sudo R CMD javareconf
RUN Rscript -e 'install.packages("dplyr", repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'
COPY ./GENOMES /home/GENOMES
COPY ./jspacker /home/jspacker
COPY ./scripts /home

