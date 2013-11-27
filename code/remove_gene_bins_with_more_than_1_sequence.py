#!/usr/bin/env python

import glob
import subprocess
import os


for i in glob.glob("output/gene*fasta"):
    cmd = "grep -c '>' " + i
    p = subprocess.check_output(cmd, shell=True)
    if int(p.strip()) > 1:
        os.remove(i)
    
