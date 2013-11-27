from pyphylogenomics import OrthoDB
from pyphylogenomics import BLAST

in_file = 'grefs/OrthoDB7_Arthropoda_tabtext'
genes = OrthoDB.single_copy_genes(in_file, 'Bombyx mori')
cds_file = 'grefs/silkcds.fa'
BLAST.get_cds(genes, cds_file)
BLAST.blastn('pulled_seqs.fasta', 'grefs/silkgenome.fa')
exons = BLAST.getLargestExon("pulled_seqs_blastn_out.csv", E_value=0.001, ident=98, exon_len=300)
exons = BLAST.eraseFalsePosi(exons)
BLAST.storeExonsInFrame(exons, "pulled_seqs.fasta", "grefs/Bombyx_exons.fas")
