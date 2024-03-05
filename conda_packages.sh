# Installing Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh

# Installing sequence manipulation programs at <base>
conda install -c bioconda seqkit
conda install -c bioconda seqtk

# Creating an environment, activating it, and installing packeages at <fastq>
## Packages for retriving fastq files and access to SRA
conda create -n fastq
conda activate fastq
conda install -c bioconda sra-tools
conda install -c bioconda parallel-fastq-dump

# Deactivating previous environment, returning to <base>
conda deactivate

# Creating an environment, activating it, and installing packeages at <trimmomatic>
## Packages for trimming and cleanup of raw sequence runs
conda create -n trimmomatic
conda activate trimmomatic
conda install -c bioconda trimmomatic

# Deactivating previous environment, returning to <base>
conda deactivate

# Creating an environment, activating it, and installing packaages at <megahit>
conda create -n megahit
conda activate megahit

## Packages for assembly of (mapped) reads
conda install -c bioconda megahit

# Deactivating previous environment, returning to <base>
conda deactivate

# Creating an environment, activating it, and installing packeages at <mafft>
## Packages for building and manupulating alignment of sequences
conda create -n mafft
conda activate mafft
conda install -c bioconda mafft
conda install -c bioconda gblocks
conda install -c bioconda amas

# Deactivating previous environment, returning to <base>
conda deactivate

# Creating an environment, activating it, and installing packeages at <repeatfinder>
## Packages for finding repeats (TDR, TIR, LINEs, etc) within sequences
conda create -n repeatfinder
conda activate repeatfinder
conda install -c bioconda genericrepeatfinder

# Deactivating previous environment, returning to <base>
conda deactivate

# Downloading other packages for NGS analysis
## Create and move to a working directory (for ease of use) in the root of driver C:
mkdir /mnt/c/Ubuntu
cd /mnt/c/Ubuntu

## Mapping of reads 
## Downloading and extracting BBMap to the working directory
wget https://github.com/BioInfoTools/BBMap/releases/download/v35.85/BBMap_35.85.tar.gz
tar -zxvf BBMap_35.85.tar.gz -C /mnt/c/Ubuntu/
rm BBMap_35.85.tar.gz

## Installing Java JDK to run BBMap in a environment of choice, such as <megahit>
conda activate megahit
conda install -c bioconda java-jdk


