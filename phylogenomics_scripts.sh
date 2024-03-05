# VERY IMPORTANT!! Do not run this as a bash shell script! 
# If you have the programs in different Conda environments, activate them accordingly

# Move into the working folder
cd /mnt/f/Ubuntu/MT_GB

### If you have individual genbank files already donwloaded...

# Compile a list of all Genbank (gb/gbk) files in the folder
ls *.gb > gb.txt

# Create a bash file to use with the pythin script ($ python script.py infile.gb)
sed '/^/ s/^/python gbex_tranlation.py /g' gb.txt > gb00_extract_aa.sh
sed '/^/ s/^/python gbex_rrnaseq.py /g' gb.txt > gb00_extract_rrna.sh

# Extract tranlation sequences from CDS entries in Genbank (GB/GBK) files
bash gb00_extract_aa.sh
bash gb00_extract_rrna.sh

# Concatenate all protein and rRNA sequences in the folder into a single file
cat *.faa > GB_all.faa
cat *.fas > GB_all.fas

### If you don't have the genbank files yet...

# Download all genbank entries at once using Batch Entrez (https://www.ncbi.nlm.nih.gov/sites/batchentrez)

# VERY IMPORTANT!! You might need to manually rename the genes for congruence among entries.

# Concatenate your Genbank files with the donwloaded ones into a single file
cat MT08.gb sequence.gb > gb_all.gb

# Extract amino acid sequences from translation features in Genbank (gb/gbk) files (custom Python script provided here)
python gbex_tranlation.py gb_all.gb

# Extract DNA sequences from rRNA features in Genbank (gb/gbk) files (custom Python script provided here)
python gbex_rrnaseq.py gb_all.gb

# Create a folder to hold temporary files (example /genes) but dont move to it (stay in the working folder set above)
mkdir genes

# Retrieve sequences from specific genes into individual files ($ bash gb01_grep.sh)
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

# To create the bash shell (BS) file...
echo 'grep "@ATP6_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP6.fa' > gb01_grep.sh
echo 'grep "@ATP8_" -A 1 --no-group-separator gb_all.faa > genes/gb_ATP8.fa' >> gb01_grep.sh
# ... and so on

# Run MAFFT to align protein (global) and rrna (motif) sequences and ($ bash gb02_mafft.sh)
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

# If you dont have a MAFFT installed...
conda install -c bioconda mafft

# If you dont have MAFFT environment activated...
conda activate mafft

# Run Gblocks to remove ambiguiously aligned positions ($ bash gb03_gblocks.sh)
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

# If you dont have Gblocks installed...
conda install -c bioconda gblocks

# If you dont have a Gblocks environment activated... (in my case, MAFFT and Gblocks are in the same environment)
conda activate gblocks

# Fix the sequences by removing empty spaces and renaming the headers; headers must be ">ID gene" for concatenation later on ($ bash gb04_fix.sh)
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

# Unwrap sequences in the FASTA files in a single-line ($ bash gb05_unwrap.sh)
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

# Concatenate the sequences by IDs
# If you choose, you can concatenate first the AA and then the DNA sequences to know the range position of each. For such specify the input files one by one (example, 'seqkit concat in1 in2 > out'.
seqkit concat genes/*_set.fa > gb_all_genes.fa

# The final alignment will be retrieved after renaming the headers, removing the list of genes
sed '/^>/ s/ .*$//g' gb_all_genes.fa > gb_all_accession.fas

# VERY IMPORTANT! The AA and DNA partitions must be defined manually.
# Retrieve lenght of the individual genes for partition definition
for i in genes/*_set.fa; do seqkit fx2tab -l -n $i > $i'-len'; done
for i in genes/*_set.fa-len; do head -n 1 $i >> gb_genes_lenght.txt ; done

# Rename the sequences in the final fasta files using a list in another file
seqkit replace -p "(.+)" -r '{kv}' -k gb_names.txt gb_all_accession.fas > gb_all_final.fas

# Run RAxML-NG to recover the best-scoring ML tree with Bootstrap support
/mnt/f/Ubuntu/raxml-ng/raxml-ng --msa gb_all_final.fas --msa-format FASTA --seed 12345 --model cat_partitions.txt --all --bs-trees autoMRE{500} --bs-metric TBE

# You can remove the temporary files in the 'genes' folder now
rm -rv genes/


