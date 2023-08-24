__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.sample) # just work if list is input -> for singleend we have to check if it is a string by isinstance(snakemake.input.sample, str)
print(snakemake.input.sample)
print(n)
assert isinstance(snakemake.input.sample, str) or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

if snakemake.input.sample[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
elif isinstance(snakemake.input.sample, str) and snakemake.input.sample.endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

outprefix = os.path.dirname(snakemake.output[0]) + "/"

print(readcmd)
print(outprefix)

# command = "STAR " + snakemake.params.star_usr  +" --limitGenomeGenerateRAM " + str(snakemake.params.ram)  +"--runThreadN " + str(snakemake.threads)  +" --genomeDir " + snakemake.input.index + " --readFilesIn " + snakemake.input.sample +" "+ readcmd + " --outFileNamePrefix " + outprefix + " " + snakemake.output[0] +" "+ log

# print(command)

shell(
    "STAR "
    "{snakemake.params.star_usr} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.input.index} "      #have to set to input! so snakemake can build it.
    "--readFilesIn {snakemake.input.sample} "
    "--outFileNamePrefix {snakemake.params.outprefix} "
    "{readcmd} "
    "{log} ")

