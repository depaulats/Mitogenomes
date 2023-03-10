#Usage: python gbex_tranlation.py sequence.gb
#Sequence headers will start, after '>', with 'locus' then 'product'; 'gene' between '@' and "_' symbols to ease GREP later on, example '>NC_010171@large subunit ribosomal RNA_'
#VERY IMPORTANT! You will need to manually rename the headers for congruence among sequences.

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
        if seq_feature.type=="rRNA" :
            output_handle.write(">%s@%s_\n%s\n" % (
                   seq_record.name,
                   seq_feature.qualifiers['product'][0],
                   seq_feature.extract(seq_record.seq)))

output_handle.close()
input_handle.close()
print("Done")
