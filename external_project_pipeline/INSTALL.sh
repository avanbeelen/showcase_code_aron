#!/bin/bash

sudo apt-get update
sudo apt-get upgrade 

# Install Python3.7
echo "Installing Python3.7..."
sudo add-apt-repository ppa:deadsnakes/ppa   
sudo apt-get update   
sudo apt install python3.7  

sudo apt-get update

# Install pip
echo "Installing pip..."
sudo apt install python3-pip
sudo apt-get update
pip3 install --upgrade pip

# Install biopython3.6
sudo apt-get update
echo "Installing biopython3.6..."
python3 -m pip install biopython

# Install bcbio-gff
sudo apt-get update
echo "Installing bcbio-gff..."
python3 -m pip install bcbio-gff

# Install Snakemake
sudo apt-get update
echo "Installing Snakemake..."
sudo apt-get install snakemake

# Install blast
echo "Installing NCBI-BLAST+..."
sudo apt-get install ncbi-blast+

# Install R
echo "Installing R..."
echo "deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/" | sudo tee -a /etc/apt/sources.list
echo "deb https://cloud.r-project.org/bin/linux/ubuntu cosmic-cran35/" | sudo tee -a /etc/apt/sources.list
echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" | sudo tee -a /etc/apt/sources.list
echo "deb https://cloud.r-project.org/bin/linux/ubuntu xenial-cran35/" | sudo tee -a /etc/apt/sources.list
echo "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" | sudo tee -a /etc/apt/sources.list
sudo apt-get update
sudo apt-get install r-base
sudo apt-get install build-essential


sudo apt-get update
sudo apt-get upgrade 


