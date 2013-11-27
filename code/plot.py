import sys
import fileinput
import csv
import prettyplotlib as ppl
from prettyplotlib import plt
import numpy as np

"""
Takes numers as input and plots a histrogram from them. Very useful to plot
columns of blast output data
"""

data = []
for line in fileinput.input():
    data.append(float(line.strip()))

fig, ax = plt.subplots(1)

ppl.hist(ax, data)
plt.xlabel(u"Percentage of identity", fontdict={'fontsize':14})
plt.ylabel(u"Frecuency", fontdict={'fontsize':14})
ax.set_title("Plotted Blast results: % of identity")

fig.savefig('plot.png')
print "Made a plot from your input data."
print "It is in the file `plot.png`"


