#Usage: python gbex_tranlation.py sequence.gb

import sys, os
from Bio import GenBank, SeqIO

gbk_filename = sys.argv[1]
root_name = os.path.splitext(gbk_filename)[0]
faa_filename = root_name + ".faa"

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s@%s_\n%s\n" % (
                   seq_record.name,
                   seq_feature.qualifiers['gene'][0],
                   seq_feature.qualifiers['translation'][0]))

#Headers will start with '>' then 'locus' then '@' then 'gene' then '_', example '>NC_010171@COX1_'

output_handle.close()
input_handle.close()
print("Done")
