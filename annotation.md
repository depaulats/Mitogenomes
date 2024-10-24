# Annotating mitogenomes

In order to annotate you mitogenomes into GenBank files (gb/gbk) you will need you sequence file (FAS/FASTA) and a tab-delimited, feature file (GFF/GTF).

## FASTA File

For compatibility issues, keep your final mtDNA assemblies in a FASTA file wrapped in strings of 60 bp.

```
seqkit seq path-to-input/input-file.fasta -w 60 > path-to-output/output-file.fasta
```

Also, keep their headers simple, with an unique ID number, such as their accession or voucher numbers, followed by organism and molecule, for instance:

```
XX000000 Genus epitet voucher MUSEUM 00000 mitochondrion, complete genome
```

## GFF file
