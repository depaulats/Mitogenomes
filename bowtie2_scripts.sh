## Installing bowtie2 via conda
conda create --name bowtie2
conda activate bowtie2
conda install bioconda::bowtie2

## Installing samtools via conda (remember to return to <basee> environment using $ conda deactivate)
conda create --name samtools
conda activate samtools
conda install bioconda::samtools

## Creating an indexed reference genome with bowtie2 for local alignment (remember to move to <bowtie2> environment using $ conda activate bowtie2)
bowtie2-build /path-to-file/fasta-file.fas /path-to-index/index-name

## Local aligning of pair-end reads to an indexed reference genome (remember to move to <bowtie2> environment using $ conda activate bowtie2)
bowtie2 --fr --local -x /path-to-index/index-name -1 /path-to-input/R1-file.fastq -2 /path-to-input/R2-file.fastq -S /path-to-output/output-file.sam

## Converting a SAM file to a BAM file using samtools (remember to move to <samtools> environment using $ conda activate samtools)
samtools view -b -S /path-to-input/input-file.sam -o /path-to-output/output-file.bam

## Extracting mapped reads from a BAM file using samtools (remember to move to <samtools> environment using $ conda activate samtools)
samtools view -b -F 4 /path-to-input/input-file.sam > /path-to-output/output-file_mapped.bam

## Sorting mapped reads in a BAM file using samtools (remember to move to <samtools> environment using $ conda activate samtools)
samtools sort /path-to-input/input-file_mapped.bam > /path-to-output/output-file_mapped_sorted.bam

## Creating an index for the sorted reads in a BAM file using samtools (remember to move to <samtools> environment using $ conda activate samtools)
samtools index /path-to-input/input-file_mapped_sorted.bam > /path-to-output/output-file_mapped_sorted.bai
