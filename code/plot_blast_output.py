import sys
import csv
import prettyplotlib as ppl
from prettyplotlib import plt
import numpy as np

blast_data = sys.argv[1].strip()

ident = []
with open(blast_data, "r") as handle:
    reader = csv.reader(handle)
    for row in reader:
        ident.append(float(row[2]))


fig, ax = plt.subplots(1)

ppl.hist(ax, ident)
plt.xlabel(u"Percentage of identity", fontdict={'fontsize':14})
plt.ylabel(u"Frecuency", fontdict={'fontsize':14})
ax.set_title("Plotted Blast results: % of identity")

fig.savefig('fig_blast_identity.png')
print "Made a plot of the percentags of identity from the blast data."
print "It is in the file `fig_blast_identity.png`"


