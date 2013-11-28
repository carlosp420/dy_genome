#!/usr/bin/env python

from pyphylogenomics import BLAST
import sys
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



if len(sys.argv) < 2:
    print """This script takes as input a FASTA file and will put all sequences
            in same direction using blast results."""
    sys.exit()

infile = sys.argv[1].strip()
BLAST.blastn(infile, infile)

# remove BLAST files
for file in glob.glob(infile + ".*"):
    os.remove(file)
if os.path.isfile(infile + "_dust.asnb"):
    os.remove(infile + "_dust.asnb")

# parse BLAST output file
blast_file = infile.replace(".fasta", "_blastn_out.csv")
f = open(blast_file, "r")
lines = f.readlines()
f.close()

def reverse(id, infile):
    # reverse complement sequence of ID in FASTA file infile
    for seq_record in SeqIO.parse(infile, "fasta"):
        if str(seq_record.id) == id:
            seq = seq_record.seq.reverse_complement()
            this_seq_record = SeqRecord(seq)
            this_seq_record.id = id

            return this_seq_record

def get_this_sequence(id, infile):
    # extract sequence from file and return as seqrecord object
    for seq_record in SeqIO.parse(infile, "fasta"):
        if str(seq_record.id) == id:
            return seq_record



sequences = []
ids = []
for line in lines:
    line = line.strip().split(",")
    if line[0] != line[1]:
        # direction of query and subject
        if int(line[7]) > int(line[8]):
            # reverse
            record = reverse(line[0], infile)
            if record.id not in ids:
                print "Appending reverse of %s" % str(line[0])
                sequences.append(record)
                ids.append(record.id)
        else:
            record = get_this_sequence(line[0], infile)
            if record.id not in ids:
                print "Appendig unchanged %s" % str(line[0])
                sequences.append(record)
                ids.append(record.id)


# write sequences to infile
print "We have %i sequences to process" % len(sequences)
f = open(infile + "sd.fasta", "w")
for i in sequences:
    if i != None:
        f.write(">" + str(i.id) + "\n")
        f.write(str(i.seq) + "\n")
        print "Writing sequence %s to file %s" % (str(i.id), str(infile))
f.close()

