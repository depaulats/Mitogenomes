# Cleaning raw reads
## Using pair-end sequences (PE), 8 threads (-threads 8), quality scores Phred+33 (-phred33)
## Replace sequence file names and locations on the R1 and R2 sequence files and corresponding paired and unpaired output files
## Trimming TruSeq3 PE adpaters (see manual for others, such as Nextera adapters (NexteraPE-PE.fa)
## Removing 3 bases at start (LEADING:3) and end (TRAILING:3)
## Removing low quality bases using 4 base average of sliding window with qualiity 20 (SLIDINGWINDOW:4:20)
## Discarding reads shorter than 50 bp (MINLEN:50)
## IMPORTANT: MIght need to edit into a single line, removing paragraph breaks and backslahes
trimmomatic PE -threads 8 -phred33 \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R1.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R2.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R1_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R1_unpaired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R2_paired.fastq.gz \
  /mnt/c/Ubuntu/Sample_folder/Sample_sequence_R2_unpaired.fastq.gz \
  ILLUMINACLIP:/mnt/c/Ubuntu/Trimmomatic/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# 
