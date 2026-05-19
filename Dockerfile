FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
    wget \
 && wget https://repo.anaconda.com/miniconda/Miniconda3-py312_25.9.1-1-Linux-x86_64.sh \
 && bash Miniconda3-py312_25.9.1-1-Linux-x86_64.sh -b -p /root/miniconda3 \
 && rm Miniconda3-py312_25.9.1-1-Linux-x86_64.sh \
 && apt-get clean

ENV PATH=/root/miniconda3/bin:$PATH

RUN conda config --add channels bioconda \
 && conda config --add channels conda-forge \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
 && conda install -y \
    star==2.7.10b \
    fuzzysearch==0.8.1 \
    samtools==1.21 \
    pandas==2.3.3 \
    matplotlib==3.10 \
    umi_tools==1.1.6 \
    seaborn==0.13.2 \
    cutadapt==5.2 \
    scanpy==1.11.5 \
    python-igraph==0.11.9 \
    leidenalg==0.10.2 \
    fastp==1.0.1 \
    bio==1.8.1 \
    seqtk==1.5 \
 && conda create -n stardist python=3.10 -y \
 && /root/miniconda3/envs/stardist/bin/pip install \
    tensorflow-cpu==2.21.0 \
    stardist==0.9.2 \
    imagecodecs==2025.3.30 \
    pandas==2.3.3 \
    opencv-python-headless==4.13.0.92 \
 && conda clean -afy