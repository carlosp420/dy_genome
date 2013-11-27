#!/usr/bin/env python

from pyphylogenomics import BLAST

BLAST.blastn("grefs/Bombyx_exons.fas", "grefs/Dp_genome_v2.fasta")
BLAST.blastParser("grefs/Bombyx_exons_blastn_out.csv",
                    "grefs/Dp_genome_v2.fasta",     
                    "grefs/Danaus_exons.fasta",
                    sp_name = "Danaus")
