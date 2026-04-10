# Image with all dependencies for AGAT
FROM ubuntu:latest

# Needed to avoid to be stuck when installing r-base
ENV DEBIAN_FRONTEND noninteractive

RUN  apt-get update

RUN apt-get -y install libdb-dev

RUN apt-get -y install r-base

# Perl dependencies
RUN apt-get -y install libbio-perl-perl libclone-perl libgraph-perl liblwp-useragent-determined-perl libstatistics-r-perl libcarp-clan-perl libsort-naturally-perl libfile-share-perl libfile-sharedir-install-perl libyaml-perl liblwp-protocol-https-perl libfile-sharedir-perl libmoose-perl libterm-progressbar-perl libdevel-cover-perl
RUN apt-get -y install libparallel-forkmanager-perl libipc-sharelite-perl libfile-chdir-perl

RUN apt-get install locales && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

RUN apt-get -y install curl
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Créer un utilisateur non-root
#RUN useradd -ms /bin/bash vscode
#USER vscode
WORKDIR /workspace

# Set environment variables
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

RUN apt-get -y install git wget nano