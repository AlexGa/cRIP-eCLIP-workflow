#!/usr/bin python
_author_ = "Alexander Gabel"

import sys
import os
import subprocess
import re
from pathlib import Path


query_bed = snakemake.input.query
target_beds = snakemake.input.target

if not os.path.exists(os.path.dirname(target_beds[0])):
	Path(os.mkdir(os.path.dirname(target_beds[0]))).mkdir(parents = True, exist_ok = True)

cmd_line = "bedtools intersect -a " + query_bed + " -b " + " ".join(target_beds)+ " > " + snakemake.output[0]
# cmd_line = "bedtools intersect -a " + query_bed + " -b " + " ".join(target_beds)+ " -u > " + snakemake.output[0]

print(cmd_line)
child = subprocess.call(cmd_line,shell=True)
print("Exit status: " + str(child))