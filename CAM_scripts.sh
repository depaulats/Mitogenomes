# Cleaning and trimming raw reads
## If needed, deactivate current environment
conda deactivate

## Activating environment to use the Trimmomatic package
conda activate trimmomatic

## Running Trimmomatic
## Using pair-end sequences (PE), 8 threads (-threads 8), quality scores Phred+33 (-phred33)
## Replace sequence file names and locations on the R1 and R2 sequence files and corresponding paired and unpaired output files
## Trimming TruSeq3 PE adpaters (see manual for others, such as Nextera adapters (NexteraPE-PE.fa)
## Removing 3 bases at start (LEADING:3) and end (TRAILING:3)
## Removing low quality bases using 4 base average of sliding window with qualiity 20 (SLIDINGWINDOW:4:20)
## Discarding reads shorter than 50 bp (MINLEN:50)
## IMPORTANT: Edit it into a single line, removing paragraph breaks and backslahes
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

# Mapping reads
## If needed, deactivate current environment
conda deactivate

## Activating environment to use the mapping and assembly package
conda activate megahit

## Running BBMap
## Using a specific FASTA file as reference (ref=file) without index (nodisk)
## Using as input the paired reads R! (in1=file) and R2 (in2=file) files
## Record only the mapped reads (mappedonly=t)
## Saving an output file with interleaved mapped reads (out=file)
## IMPORTANT: Edit it into a single line, removing paragraph breaks and backslahes
/mnt/c/Ubuntu/bbmap/bbmap.sh \
  ref=/mnt/c/Ubuntu/reference.fasta nodisk \
  in1=/mnt/c/Ubuntu/Sample_folder/Sample_file_R1_paired.fastq.gz \
  in2=/mnt/c/Ubuntu/Sample_folder/Sample_file_R2_paired.fastq.gz \
  mappedonly=t \
  out=/mnt/c/Ubuntu/Sample_folder/Sample_file_mapped_data.fastq

# Assembling reads
## If needed, deactivate current environment
conda deactivate

## Activating environment to use the mapping and assembly package
conda activate megahit

## Running Megahit
## Using an interleaved file as input (--12 file)
## Saving output files in a specific directory (--out-dir directory)
## IMPORTANT: Edit it into a single line, removing paragraph breaks and backslahes
megahit --12 /mnt/c/Ubuntu/Sample_folder/Sample_file_mapped_data.fastq \
  --out-dir /mnt/c/Sample_folder/Sample_file_mapped_data

