from pyphylogenomics import NGS;
import sys

blast_table = sys.argv[1].strip()
ion_file    = sys.argv[2].strip()
NGS.parse_blast_results(blast_table, ion_file);
