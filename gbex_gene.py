#Usage: python gbex_gene.py sequence.gb
#Do not work with 'locus_tag'
#Sequence headers will start, after '>', with 'locus' then 'gene'; 'gene' between '@' and "_' symbols to ease GREP later on, example '>NC_010171@COX1_'

import sys, os
from Bio import GenBank, SeqIO

gbk_filename = sys.argv[1]
root_name = os.path.splitext(gbk_filename)[0]
faa_filename = root_name + ".fas"

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="gene" :
            output_handle.write(">%s@%s_\n%s\n" % (
                   seq_record.name,
                   seq_feature.qualifiers['gene'][0],
                   seq_feature.extract(seq_record.seq)))

output_handle.close()
input_handle.close()
print("Done")
