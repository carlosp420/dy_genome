from pyphylogenomics import BLAST;
import sys

query_seqs = sys.argv[1].strip()
genome = sys.argv[2].strip()
BLAST.blastn(query_seqs, genome);

