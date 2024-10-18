# 

## Download reads from the Sequence Read Archive (SRA)
If needed, deactivate current environment and activate the environment to use the [***SRA Tools***](https://github.com/ncbi/sra-tools) package .

```
conda deactivate
conda activate sratools
```

Run *fasterq-dump* using the settings:
- Retrieve FASTQ files from the run `SRA0000000`;
- Save output (`-O`) in the folder `/mnt/c/Ubuntu/Sample_folder/`;
- Keep R1 and R2 reads in different files (`--split-files`);
- Use `8` logical processors (`--threads`);
- Show progress bar (`--threads`).

```
fasterq-dump SRA0000000 -O /mnt/c/Ubuntu/Sample_folder/ --split-files --threads 8 --progress
```

Compress your FASTQ files to save space.
```
gzip /mnt/c/Ubuntu/Sample_folder/SRA0000000_1.fastq
gzip /mnt/c/Ubuntu/Sample_folder/SRA0000000_1.fastq
```

## Cleaning and trimming raw reads

If needed, deactivate current environment and activate the environment to use the [***Trimmomatic***](https://github.com/usadellab/Trimmomatic) software.

```
conda deactivate
conda activate trimmomatic
```

\
Run *trimmomatic* using the settings:

- Using pair-end sequences (`PE`), 8 threads (`-threads 8`), quality scores Phred+33 (`-phred33`);
- Input R1 and R2 sequence files (with /path/) and corresponding paired and unpaired output files;
- Trimming TruSeq3 PE adpaters (`TruSeq3-PE-2.fa`; see manual for others, such as Nextera adapters, `NexteraPE-PE.fa`);
- Removing 3 bases at start (`LEADING:3`) and end (`TRAILING:3`);
- Removing low quality bases using 4 base average of sliding window with qualiity 20 (`SLIDINGWINDOW:4:20`);
- Discarding reads shorter than 50 bp (`MINLEN:50`).

**IMPORTANT**: Edit the code into a single line, removing paragraph breaks and backslahes (\).
  
```
trimmomatic PE -threads 8 -phred33 \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R1.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R2.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R1_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R1_unpaired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R2_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R2_unpaired.fastq.gz \
  ILLUMINACLIP:/mnt/c/Ubuntu/Trimmomatic/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 \
  SLIDINGWINDOW:4:20 \
  MINLEN:50
```


## Mapping reads
If needed, deactivate current environment and activate the environment to use the [***BBMap***](https://github.com/BioInfoTools/BBMap) package.

```
conda deactivate
conda activate megahit
```

\
Run *bbmap.sh* with to following settings:

- Using a specific FASTA file as reference (ref=file) without index (`nodisk`);
- Using as input the paired reads R! (`in1=file`) and R2 (`in2=file`) files;
- Recording only the mapped reads (`mappedonly=t`);
- Saving an output file with interleaved mapped reads (`out=file`).

**IMPORTANT**: Edit the code into a single line, removing paragraph breaks and backslahes (\).

```
/mnt/c/Ubuntu/bbmap/bbmap.sh \
  ref=/mnt/c/Ubuntu/reference.fasta nodisk \
  in1=/mnt/c/Ubuntu/Sample_folder/Sample_file_R1_paired.fastq.gz \
  in2=/mnt/c/Ubuntu/Sample_folder/Sample_file_R2_paired.fastq.gz \
  mappedonly=t \
  out=/mnt/c/Ubuntu/Sample_folder/Sample_file_mapped_data.fastq
```

## Assembling reads
If needed, deactivate current environment and activate the environment to use the [***Megahit***](https://github.com/voutcn/megahit) software.

```
conda deactivate
conda activate megahit
```

\
Run *megahit* with the following settings:

- Using an interleaved file as input (`--12 file`);
- Saving output files in a specific directory (`--out-dir directory`).

**IMPORTANT**: Edit the code into a single line, removing paragraph breaks and backslahes (\).
  
```
megahit --12 /mnt/c/Ubuntu/Sample_folder/Sample_file_mapped_data.fastq \
  --out-dir /mnt/c/Sample_folder/Sample_file_mapped_data
```

## Reads statistics
If needed, deactivate current environment to use the [***SeqKit***](https://bioinf.shenwei.me/seqkit/) package at <base>.

```
conda deactivate
```

\
Generate statiscts for raw, filtered, and mapped reads. 

**IMPORTANT**: Edit the code into a single line, removing paragraph breaks and backslahes (\).

```
seqkit stats -a \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R1.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R2.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R1_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_R2_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_file_mapped_data.fastq
```
