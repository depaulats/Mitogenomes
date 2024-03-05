# Extracting biological data from Genbank files

\
**IMPORTANT**: If you have the softwares in different Conda environments, activate them accordingly.

\
Move into the working folder, for instance:
```
cd /mnt/c/Ubuntu/MT_GB
```

## If you have already individual Genbank files...

\
Compile a list of all Genbank (gb/gbk) files in the folder.
```
ls *.gb > gb.txt
```

\
Create a bash file to use with the python scripts (`$ python script.py infile.gb`)
```
sed '/^/ s/^/python gbex_tranlation.py /g' gb.txt > gb00_extract_aa.sh
sed '/^/ s/^/python gbex_rrnaseq.py /g' gb.txt > gb00_extract_rrna.sh
```

\
Extract tranlation sequences from CDS entries in Genbank (gb/gbk) files.
```
bash gb00_extract_aa.sh
bash gb00_extract_rrna.sh
```

\
Concatenate all protein and rRNA sequences in the folder into a single file.
```
cat *.faa > GB_all.faa
cat *.fas > GB_all.fas
```

## If you don't have the Genbank files yet...

\
Download all Genbank entries at once using [Batch Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez)

\
**IMPORTANT**: In order to retrieve the sequences using the codes below, you might need to manually rename the name/product of genes within genbank entries for congruence 
(*e.g.*, replacing "**nad1**" to "**ND1**" or "**COI**" to "**COX1**")

\
Concatenate your Genbank files (*e.g.*, MT08.gb) with new ones you downloaded into a single file.
```
cat MT08.gb sequence.gb > gb_all.gb
```

\
Extract amino acid sequences from translation features in Genbank (gb/gbk) files using the Python script [gbex_tranlation.py](https://github.com/depaulats/Mitogenomes/blob/main/gbex_tranlation.py).
```
python gbex_tranlation.py gb_all.gb
```

\
Extract DNA sequences from rRNA features in Genbank (gb/gbk) files using the Python script [gbex_rrnaseq.py](https://github.com/depaulats/Mitogenomes/blob/main/gbex_rrnaseq.py).
python gbex_rrnaseq.py gb_all.gb


# Splitting data for individual gene aligment

\
Create a folder to hold temporary files (*e.g.*, genes) but dont move to it, staying in the working folder set above. 
```
mkdir genes
```

\
Retrieve sequences from specific genes into individual files. 
```
grep "@ATP6_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP6.fa
grep "@ATP8_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP8.fa
grep "@ATP9_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP9.fa
grep "@COX1_" -A 1 --no-group-separator gb_all.faa > genes/gb_COX1.fa
grep "@COX2_" -A 1 --no-group-separator gb_all.faa > genes/gb_COX2.fa
grep "@COX3_" -A 1 --no-group-separator gb_all.faa > genes/gb_COX3.fa
grep "@CYTB_" -A 1 --no-group-separator gb_all.faa > genes/gb_CYTB.fa
grep "@ND1_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND1.fa
grep "@ND2_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND2.fa
grep "@ND3_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND3.fa
grep "@ND4_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND4.fa
grep "@ND4L_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND4L.fa
grep "@ND5_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND5.fa
grep "@ND6_" -A 1 --no-group-separator gb_all.faa > genes/gb_ND6.fa
grep "@rnl_" -A 1 --no-group-separator gb_all.fas > genes/gb_rnl.fa
grep "@rns_" -A 1 --no-group-separator gb_all.fas > genes/gb_rns.fa
```

\
If you prefer, create a bash shell (*e.g.*, gb01_grep.sh) from the code above using the ***echo*** command, as follow:
```
echo 'grep "@ATP6_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP6.fa' > gb01_grep.sh
echo 'grep "@ATP8_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP8.fa' >> gb01_grep.sh
...
```

And run it using the ***bash*** command.
```
bash gb01_grep.sh
```

# Multiple sequence alignment (MSA)

If you dont have the [MAFFT](https://mafft.cbrc.jp/alignment/software/) sofware installed, install it using Conda on its own enviroment.
```
conda create --name mafft
conda activate mafft
conda install -c bioconda mafft
```

\
Run the MAFFT software to align protein (global) and rRNA (motif) sequences.
```
mafft --globalpair --maxiterate 1000 genes/gb_ATP6.fa > genes/gb_ATP6_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ATP8.fa > genes/gb_ATP8_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ATP9.fa > genes/gb_ATP9_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_COX1.fa > genes/gb_COX1_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_COX2.fa > genes/gb_COX2_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_COX3.fa > genes/gb_COX3_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_CYTB.fa > genes/gb_CYTB_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND1.fa > genes/gb_ND1_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND2.fa > genes/gb_ND2_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND3.fa > genes/gb_ND3_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND4.fa > genes/gb_ND4_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND4L.fa > genes/gb_ND4L_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND5.fa > genes/gb_ND5_align.fa
mafft --globalpair --maxiterate 1000 genes/gb_ND6.fa > genes/gb_ND6_align.fa
mafft --genafpair --maxiterate 1000 genes/gb_rnl.fa > genes/gb_rnl_align.fa
mafft --genafpair --maxiterate 1000 genes/gb_rns.fa > genes/gb_rns_align.fa
```

If you prefer, create a bash shell (*e.g.*, gb02_mafft.sh) from the code above using the ***echo*** command, as follow:
```
echo 'mafft --globalpair --maxiterate 1000 genes/gb_ATP6.fa > genes/gb_ATP6_align.fa' > gb02_mafft.sh
echo 'mafft --globalpair --maxiterate 1000 genes/gb_ATP8.fa > genes/gb_ATP8_align.fa' >> gb02_mafft.sh
...
```

And run it using the ***bash*** command.
```
bash gb02_mafft.sh
```

# Removing ambiguously aligned positions from MSA (optional)

\
If you dont have the [Gblocks](https://github.com/atmaivancevic/Gblocks) sofware installed, install it using Conda on the <mafft> enviroment.
```
conda install -c bioconda gblocks
```

\
Run the [Gblocks](https://github.com/atmaivancevic/Gblocks) software to remove ambiguiously aligned positions.
```
Gblocks genes/gb_ATP6_align.fa -t=p -b5=a
Gblocks genes/gb_ATP8_align.fa -t=p -b5=a
Gblocks genes/gb_ATP9_align.fa -t=p -b5=a
Gblocks genes/gb_COX1_align.fa -t=p -b5=a
Gblocks genes/gb_COX2_align.fa -t=p -b5=a
Gblocks genes/gb_COX3_align.fa -t=p -b5=a
Gblocks genes/gb_CYTB_align.fa -t=p -b5=a
Gblocks genes/gb_ND1_align.fa -t=p -b5=a
Gblocks genes/gb_ND2_align.fa -t=p -b5=a
Gblocks genes/gb_ND3_align.fa -t=p -b5=a
Gblocks genes/gb_ND4_align.fa -t=p -b5=a
Gblocks genes/gb_ND4L_align.fa -t=p -b5=a
Gblocks genes/gb_ND5_align.fa -t=p -b5=a
Gblocks genes/gb_ND6_align.fa -t=p -b5=a
Gblocks genes/gb_rnl_align.fa -t=d -b5=a
Gblocks genes/gb_rns_align.fa -t=d -b5=a
```

 ($ bash gb03_gblocks.sh)

If you prefer, create a bash shell (*e.g.*, gb03_gblocks.sh) from the code above using the ***echo*** command, as follow:
```
echo 'Gblocks genes/gb_ATP6_align.fa -t=p -b5=a' > gb03_gblocks.sh
echo 'Gblocks genes/gb_ATP8_align.fa -t=p -b5=a' >> gb03_gblocks.sh
...
```

And run it using the ***bash*** command.
```
bash gb03_gblocks.sh
```

# Concatenating individual aligments into a single matrix

\
Fix the sequences by removing empty spaces and renaming the headers; headers must be ">ID gene" for concatenation later on.
```
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ATP6_align.fa-gb > genes/gb_ATP6_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ATP8_align.fa-gb > genes/gb_ATP8_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ATP9_align.fa-gb > genes/gb_ATP9_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_COX1_align.fa-gb > genes/gb_COX1_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_COX2_align.fa-gb > genes/gb_COX2_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_COX3_align.fa-gb > genes/gb_COX3_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_CYTB_align.fa-gb > genes/gb_CYTB_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND1_align.fa-gb > genes/gb_ND1_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND2_align.fa-gb > genes/gb_ND2_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND3_align.fa-gb > genes/gb_ND3_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND4_align.fa-gb > genes/gb_ND4_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND4L_align.fa-gb > genes/gb_ND4L_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND5_align.fa-gb > genes/gb_ND5_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ND6_align.fa-gb > genes/gb_ND6_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_rnl_align.fa-gb > genes/gb_rnl_fix.fa
sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_rns_align.fa-gb > genes/gb_rns_fix.fa
```

\
If you prefer, create a bash shell (*e.g.*, gb04_fix.sh) from the code above using the ***echo*** command, as follow:
```
echo 'sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ATP6_align.fa-gb > genes/gb_ATP6_fix.fa' > gb04_fix.sh
echo 'sed -e 's/ //g; /^>/ s/@/ /g; /^>/ s/_$//g' genes/gb_ATP8_align.fa-gb > genes/gb_ATP8_fix.fa' >> gb04_fix.sh
...
```

And run it using the ***bash*** command.
```
bash gb04_fix.sh
```

\
Unwrap sequences in FASTA files in a single line
```
seqkit seq genes/gb_ATP6_fix.fa -w 0 > genes/gb_ATP6_set.fa
seqkit seq genes/gb_ATP8_fix.fa -w 0 > genes/gb_ATP8_set.fa
seqkit seq genes/gb_ATP9_fix.fa -w 0 > genes/gb_ATP9_set.fa
seqkit seq genes/gb_COX1_fix.fa -w 0 > genes/gb_COX1_set.fa
seqkit seq genes/gb_COX2_fix.fa -w 0 > genes/gb_COX2_set.fa
seqkit seq genes/gb_COX3_fix.fa -w 0 > genes/gb_COX3_set.fa
seqkit seq genes/gb_CYTB_fix.fa -w 0 > genes/gb_CYTB_set.fa
seqkit seq genes/gb_ND1_fix.fa -w 0 > genes/gb_ND1_set.fa
seqkit seq genes/gb_ND2_fix.fa -w 0 > genes/gb_ND2_set.fa
seqkit seq genes/gb_ND3_fix.fa -w 0 > genes/gb_ND3_set.fa
seqkit seq genes/gb_ND4_fix.fa -w 0 > genes/gb_ND4_set.fa
seqkit seq genes/gb_ND4L_fix.fa -w 0 > genes/gb_ND4L_set.fa
seqkit seq genes/gb_ND5_fix.fa -w 0 > genes/gb_ND5_set.fa
seqkit seq genes/gb_ND6_fix.fa -w 0 > genes/gb_ND6_set.fa
seqkit seq genes/gb_rnl_fix.fa -w 0 > genes/gb_rnl_set.fa
seqkit seq genes/gb_rns_fix.fa -w 0 > genes/gb_rns_set.fa
```

\
If you prefer, create a bash shell (*e.g.*, gb05_unwrap.sh) from the code above using the ***echo*** command, as follow:
```
echo 'seqkit seq genes/gb_ATP6_fix.fa -w 0 > genes/gb_ATP6_set.fa' > gb05_unwrap.sh
echo 'seqkit seq genes/gb_ATP8_fix.fa -w 0 > genes/gb_ATP8_set.fa' >> gb05_unwrap.sh
...
```

And run it using the ***bash*** command.
```
bash gb05_unwrap.sh
```

\
Run the [seqkit](https://github.com/shenwei356/seqkit) software to concatenate the sequences by their IDs. 
```
seqkit concat genes/*_set.fa > gb_all_genes.fa
```

If you choose, you can concatenate first the AA and then the DNA sequences to know the range position of each partition, 
specifying the order of each input files, one by one, as for instance:
```
seqkit concat in1 in2 > out
```

\
The final matrix will be obtained after renaming the headers, removing the list of genes.
```
sed '/^>/ s/ .*$//g' gb_all_genes.fa > gb_all_accession.fas
```

\
**VERY IMPORTANT**: The AA and DNA partitions must be defined manually.

\
For such, retrieve the lenght of the individual aligments in order to define the partitions later on.
```
for i in genes/*_set.fa; do seqkit fx2tab -l -n $i > $i'-len'; done
for i in genes/*_set.fa-len; do head -n 1 $i >> gb_genes_lenght.txt ; done
```

\
Rename the sequences in the final fasta files using a list of names in another file (*e.g.*, gb_names.txt). The list must have the accession number and the name of the entry 
separated by a tab, for instance `XX000000  XX000000_Genus_epitheton`
```
seqkit replace -p "(.+)" -r '{kv}' -k gb_names.txt gb_all_accession.fas > gb_all_final.fas
```

\
Now you retrieved the final matrix, it is safe to remove the temporary files in the 'genes' folder.

**IMPORTANT**: You can always generate them again starting from [here](#splitting-data-for-individual-gene-aligment).
```
rm -rv genes/
```

# Reconstructing a phylogenetic tree


If you dont have the [RAxML-NG](https://github.com/amkozlov/raxml-ng/) sofware installed, install it using Conda on its own enviroment.
```
conda create --name raxml
conda activate raxml
conda install bioconda::raxml-ng
```

Or download the pre-compiled binary from the GiHub page. 

**IMPORTANT**: Check the version of the file before running the code.

```
wget https://github.com/amkozlov/raxml-ng/releases/download/1.2.1/raxml-ng_v1.2.1_linux_x86_64.zip
unzip raxml-ng_v1.2.1_linux_x86_64.zip -d /mnt/c/Ubuntu/
rm raxml-ng_v1.2.1_linux_x86_64.zip
```

Create a text file to hold model and partition data (*e.g.*, cat_partitions.txt). The file present two partitions:
- Protein coding genes (PCGs), using the model `mtZOA+G10+FO`, ranging from `1` to `3830` bp; and
- Ribosomal RNAs (rRNAs), using the model `GTR+G10+FO`, ranging from `3831` to `6877` bp.

- Using pair-end sequences (PE), 8 threads (-threads 8), quality scores Phred+33 (-phred33);
- Replace sequence file names and locations on the R1 and R2 sequence files and corresponding paired and unpaired output files;The code below, t

**IMPORTANT**: Refer to the [RAxML-NG Wiki](https://github.com/amkozlov/raxml-ng/wiki) to specify the evolutionary models.

**IMPORTANT**: Refer to the file holding length data of the individual alignments (e.g., gb_genes_lenght.txt) to set partition sizes.

```
echo 'mtZOA+G10+FO, PCGs = 1-3830' > cat_partitions.txt
echo 'GTR+G10+FO, rRNAs = 3831-6877' >> cat_partitions.txt
```

\
Run RAxML-NG to recover the best-scoring ML tree with Bootstrap support. 
**IMPORTANT**: The code below is running the binary from the downloaded file. In order to run the binary installed via Conda, remove the path `/mnt/c/Ubuntu/` from the code.
```
/mnt/c/Ubuntu/raxml-ng/raxml-ng --msa gb_all_final.fas --msa-format FASTA --seed 12345 --model cat_partitions.txt --all --bs-trees autoMRE{500} --bs-metric TBE
```


