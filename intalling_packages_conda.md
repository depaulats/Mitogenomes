# Installation of software and packages in different environments

## Intalling Miniconda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
rm Miniconda3-latest-Linux-x86_64.sh
```

##  Sequence manipulation 
Installing [***Seqkit***](https://github.com/shenwei356/seqkit) and [***Seqtk***](https://github.com/lh3/seqtk) packages at `<base>`
```
conda install -c bioconda seqkit
conda install -c bioconda seqtk
```

## Retriving FASTQ files from SRA
Creating and activating an environment, and installing [***SRA Tools***](https://github.com/ncbi/sra-tools) and [***parallel-fastq-dump***](https://github.com/rvalieris/parallel-fastq-dump) packages at `<fastq>`
```
conda create --name fastq
conda activate fastq
conda install -c bioconda sra-tools
conda install -c bioconda parallel-fastq-dump
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

## Trimming and cleaning up raw reads
Creating and activating an environment, and installing [***Trimmomatic***](https://github.com/usadellab/Trimmomatic) packages at `<trimmomatic>`
```
conda create --name trimmomatic
conda activate trimmomatic
conda install -c bioconda trimmomatic
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

## Assembly of reads (mapped or *de novo*)
Creating and activating an environment, and installing [***Megahit***](https://github.com/voutcn/megahit) packaages at `<megahit>`
```
conda create --name megahit
conda activate megahit
conda install -c bioconda megahit
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

## Mapping reads to a reference

Downloading and extracting [***BBMap***](https://sourceforge.net/projects/bbmap/) to a working directory (in Windows `C:\`, in Linux `/mnt/c/`). Mind the version of the software in the code.
```
mkdir /mnt/c/Ubuntu
cd /mnt/c/Ubuntu
wget https://sourceforge.net/projects/bbmap/files/BBMap_39.06.tar.gz
tar -zxvf BBMap_39.06.tar.gz -C /mnt/c/Ubuntu/
rm BBMap_39.06.tar.gz
```

Installing [***Open JDK***](https://github.com/openjdk/) to run ***BBMap*** in a environment of choice, such as `<megahit>`.
```
conda activate megahit
conda install -c bioconda java-jdk
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

Creating and activating an environment, and installing [***Bowtie2***](https://github.com/BenLangmead/bowtie2) packaages at `<bowtie2>`
```
conda create --name bowtie2
conda activate bowtie2
conda install -c bioconda bowtie2
```


## Handling files with mapped reads and variant calls
Creating and activating an environment, and installing [***Samtools***](https://github.com/samtools/samtools) packaages at `<samtools>`
```
conda create --name samtools
conda activate samtools
conda install -c bioconda samtools
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

## Building and handling multiple sequence alignments (MSA)
Creating and activating an environment, and installing [***MAFFT***](https://mafft.cbrc.jp/alignment/software/), [***GBlocks***](https://github.com/atmaivancevic/Gblocks), and [***AMAS***](https://github.com/marekborowiec/AMAS) packages at `<mafft>`
```
conda create --name mafft
conda activate mafft
conda install -c bioconda mafft
conda install -c bioconda gblocks
conda install -c bioconda amas
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

## Finding repeats (TDR, TIR, LINEs, etc) within sequences
Creating and activating an environment, and installing [***Generic Repeat***](https://github.com/bioinfolabmu/GenericRepeatFinder) Finder packeages at `<repeatfinder>`
```
conda create --name repeatfinder
conda activate repeatfinder
conda install -c bioconda genericrepeatfinder
```

Deactivating the environment, returning to `<base>`
```
conda deactivate
```

