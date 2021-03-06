import sys
import fileinput
import csv
import prettyplotlib as ppl
from prettyplotlib import plt
import numpy as np

"""
Takes numers as input and plots a barbplo from them. Very useful to plot
columns of blast output data
"""

data = []
for line in fileinput.input():
    data.append(float(line.strip()))

fig, ax = plt.subplots(1)

ppl.bar(ax, np.arange(len(data)), data)
plt.xlabel(u"gene number", fontdict={'fontsize':14})
plt.ylabel(u"Number of mismatches", fontdict={'fontsize':14})
ax.set_title("Number of mismatches between \nfound and expected gene sequences")

fig.savefig('plot.png')
print "Made a plot from your input data."
print "It is in the file `plot.png`"


