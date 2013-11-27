#!/usr/bin/env python

from pyphylogenomics import BLAST

BLAST.blastn("grefs/Bombyx_exons.fas", "grefs/Dp_genome_v2.fasta")
