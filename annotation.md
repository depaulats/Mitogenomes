# Annotating mitogenomes

In order to annotate you mitogenomes into GenBank files (gb/gbk) you will need you sequence file (FAS/FASTA) and a tab-delimited, feature file (GFF/GTF).

## Sequence files

For compatibility issues, keep your final mtDNA assemblies in a FASTA file wrapped in strings of 60 bp.

```
seqkit seq /path-to-input/input-file.fas -w 60 > /path-to-output/output-file.fas
```

Also, keep their headers simple, with an unique ID number, such as their accession or voucher numbers, followed by organism and molecule, for instance:

```
XX000000 Genus epitet voucher MUSEUM 00000 mitochondrion, complete genome
```

## Feature files


## GenBank files

You can now combine your sequence and feature files to create a GenBank (GB/GBK) file with [***Seqret***](https://galaxy-iuc.github.io/emboss-5.0-docs/seqret.html). 

```
seqret -sequence /path-to-sequence/input-sequence.fas -feature -fformat gff -fopenfile /path-to-feature/input-feature.gff -osformat genbank -osname_outseq /path-to-output/output-file.gb -ofdirectory_outseq gbk_file -auto
```
