FROM ubuntu:20.10

LABEL author="Thilo Schulze, Clemens Mauksch"

RUN apt-get update --fix-missing && apt-get install sudo
RUN sudo apt-get update --fix-missing

RUN sudo apt-get install -y apt-utils openjdk-14-jdk default-jre curl python3 cutadapt fastqc wget git ncbi-blast+ trim-galore --fix-missing
RUN sudo apt-get install -y pigz

RUN adduser --disabled-password --gecos '' symbiontdivider
RUN adduser symbiontdivider sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER symbiontdivider
WORKDIR /home/symbiontdivider/Docker
RUN sudo  chmod a+rwx /home/symbiontdivider/Docker

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
RUN rm Miniconda3-latest-Linux-x86_64.sh 

ENV PATH /home/symbiontdivider/miniconda3/bin:$PATH

RUN conda update conda
RUN conda update --all

RUN conda install -c anaconda pip
RUN pip install better_fasta_grep