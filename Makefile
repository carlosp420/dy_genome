SRC = $(wildcard *.md)

DOCS = $(SRC:.md=.docx)
PDFS = $(SRC:.md=.pdf)

doc: clean $(DOCS)

pdf: clean $(PDFS)

%.docx: %.md refs.bib
	pandoc -f markdown -V geometry:margin=1in -t docx $< --bibliography=refs.bib -o $@

%.pdf: %.md header.latex refs.bib
	pandoc --latex-engine=xelatex -s -S --template header.latex -f markdown -V geometry:margin=1in $< --bibliography=refs.bib -o $@



analysis: grefs/Bombyx_exons.fas data/DpleKU_DAS5_blastn_out.csv fig_blast_identity.png data/new_blastn_out.csv data/DpleKU_DAS5.fa output/gene*fasta grefs/Danaus_exons.fasta output/merged*fasta output/*aligned.fasta output/*trimmed output/recovered_genes/*fasta output/recovered_genes/mismatches.csv

grefs/Bombyx_exons.fas: grefs/silkgenome.fa grefs/silkcds.fa grefs/OrthoDB7_Arthropoda_tabtext code/search_genes_from_Bmori.py
	python code/search_genes_from_Bmori.py

# prepare data
data/DpleKU_DAS5.fa: data/DpleKU_DAS5.fa.orig
	rm -rf tmp.fa
	fasta_formatter -w 0 -i data/DpleKU_DAS5.fa.orig -o data/tmp.fa
	cat data/tmp.fa | sed -r 's/^>.+_contig_/>/g' > data/DpleKU_DAS5.fa

data/DpleKU_DAS5_blastn_out.csv: code/blast_against_Bmori_genes.py data/DpleKU_DAS5.fa grefs/Bombyx_exons.fas
	python code/blast_against_Bmori_genes.py data/DpleKU_DAS5.fa grefs/Bombyx_exons.fas

# Make a pretty histogram
fig_blast_identity.png: data/DpleKU_DAS5_blastn_out.csv code/plot_blast_output.py 
	python code/plot_blast_output.py data/DpleKU_DAS5_blastn_out.csv
	
data/new_blastn_out.csv: data/DpleKU_DAS5_blastn_out.csv
	# We want matches longer than 300bp
	cat data/DpleKU_DAS5_blastn_out.csv | awk -F ',' '{ if( $$4 > 300 ) { print $$0 } }' > data/tmp.csv
	# Let us select matches with NO gaps
	cat data/tmp.csv |  awk -F ',' '{ if( $$6 < 1 ) {print $$0}}' > data/new_blastn_out.csv

output/gene%fasta: code/parse_blast_results.py data/new_blastn_out.csv data/DpleKU_DAS5.fa
	python code/parse_blast_results.py data/new_blastn_out.csv data/DpleKU_DAS5.fa
	# Let us take only those bins with one sequence inside
	python code/remove_gene_bins_with_more_than_1_sequence.py

# get the sequences for Danaus from the set of Bmori genes
grefs/Danaus_exons.fasta: grefs/Bombyx_exons.fas grefs/Dp_genome_v2.fasta code/get_danaus_seqs.py
	python code/get_danaus_seqs.py

output/merged%fasta: output/gene*fasta code/merge_files_with_homologous_seqs.py grefs/Danaus_exons.fasta
	python code/merge_files_with_homologous_seqs.py

# Use BLAST to get sequences in the same direction
output/%sd.fasta: code/get_seqs_same_direction.py output/merged*fasta
	ls output/merged_gene_BG* | parallel -I {} python code/get_seqs_same_direction.py {}

# align sequences
output/%aligned.fasta: output/*sd.fasta
	ls output/*sd.fasta | parallel -I {} sed -i 's/N/-/g' {}
	ls output/*sd.fasta | parallel -I {} clustalw -TYPE=DNA -INFILE={} -OUTPUT=FASTA -OUTFILE={}_aligned.fasta	

# trim sequences using trimal
output/%trimmed: output/*aligned.fasta
	ls output/*aligned.fasta | parallel -I {} trimal -nogaps -in {} -out {}_trimmed

output/recovered_genes/%fasta: output/*trimmed
	mkdir -p output/recovered_genes
	ls output/*trimmed | parallel -I {} cp {} output/recovered_genes/{/}
	rename 's/merged_gene_(.+):.+/$$1.fasta/g' output/recovered_genes/*

# get mismatches between sequences
output/recovered_genes/mismatches.csv: output/recovered_genes/*fasta code/compare_seqs.py
	ls output/recovered_genes/*fasta | parallel -I {} python code/compare_seqs.py {} > output/recovered_genes/mismatches.csv

clean:
	rm -rf *pdf *docx
	rm -rf output/*
	rm -rf data/*csv
	rm -rf data/*fa
	rm -rf grefs/silkgenome.fa.*
	rm -rf grefs/silkgenome.fa_*
	rm -rf grefs/Bombyx_exons.fas.*
	rm -rf grefs/Bombyx_exons.fas_*
