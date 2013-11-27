#!/usr/bin/env python

import glob
import os
import string
from Bio import SeqIO
import shutil


for i in glob.glob("output/gene*fasta"):
    gene_id = os.path.basename(i.strip()).split(":")[0]
    gene_id = gene_id.replace("gene_", "")
    for seq_record in SeqIO.parse("grefs/Danaus_exons.fasta", "fasta"):
        danaus_id = seq_record.id
        danaus_id = danaus_id.split(":")[0]
        if gene_id == danaus_id:
            # joing sequences into one file
            new_file = i.replace("gene", "merged_gene")
            shutil.copyfile(i, new_file)
            f = open(new_file, "a")
            f.write("\n" + ">" + str(seq_record.id) + "\n" + str(seq_record.seq))
            f.close()
            print "Copying sequence %s into file %s" % (str(seq_record.id),
                            new_file)
