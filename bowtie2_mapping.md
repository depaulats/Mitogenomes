# Mapping reads using Bowtie2

## Installing essential softwares
If you haven't yet, install [***Bowtie2***](https://github.com/BenLangmead/bowtie2) and [***Samtools***](https://github.com/samtools/samtools) via Conda from `<base>`. 
```
conda create --name bowtie2
conda activate bowtie2
conda install bioconda::bowtie2
conda deactivate
conda create --name samtools
conda activate samtools
conda install bioconda::samtools
conda deactivate
```

## Mapping reads to a reference
Activate the `<bowtie2>` environment.
```
conda activate bowtie2
```

Run ***bowtie2-build*** to creating an indexed reference genome for the local alignment.
```
bowtie2-build /path-to-file/fasta-file.fas /path-to-index/index-name
```

Run ***bowtie2*** with the following settings:
- Run a local alignment of reads (`--local`) where R1 appears upstream of the reverse complement of R2 (`--fr`);
- Use the reference genome (`-x`) in the index `/path-to-index/index-name`;
- Use as input the paired reads R1 (`-1`) and R2 (`-2`) in the FASTQ files `/path-to-input/R1-file.fastq` and R2 `/path-to-input/R2-file.fastq`, respectively;
- Save the output file (`-S`) in the SAM file `/path-to-output/output-file.sam`.
```
bowtie2 --local --fr \
  -x /path-to-index/index-name \
  -1 /path-to-input/R1-file.fastq \
  -2 /path-to-input/R2-file.fastq \
  -S /path-to-output/output-file.sam
```

Deactivate the environment, returning to `<base>`.
```
conda deactivate
```


## Recovering mapped reads
Activate the `<samtools>` environment.
```
conda activate samtools
```

Convert the SAM file to a BAM file.
```
samtools view -b -S /path-to-input/input-file.sam -o /path-to-output/output-file.bam
```

Extract mapped reads from the BAM file.
```
samtools view -b -F 4 /path-to-input/input-file.sam > /path-to-output/output-file_mapped.bam
```

Sort the mapped reads in the BAM file.
```
samtools sort /path-to-input/input-file_mapped.bam > /path-to-output/output-file_mapped_sorted.bam
```

Creat an index for the sorted reads in the BAM file.
```
samtools index /path-to-input/input-file_mapped_sorted.bam > /path-to-output/output-file_mapped_sorted.bai
```

Deactivate the environment, returning to `<base>`.
```
conda deactivate
```
