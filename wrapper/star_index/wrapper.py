__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"

import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

gtf_command = ""
if len(snakemake.input) == 2:
	gtf_command = "--sjdbGTFfile " + snakemake.input.gtf

shell("rm -fr {snakemake.params.dir}")
		
shell("mkdir {snakemake.params.dir}")    

shell(
	"STAR "
	"--runMode genomeGenerate "
    "--limitGenomeGenerateRAM {snakemake.params.ram} "
	"--genomeDir {snakemake.params.dir} "
	"--genomeFastaFiles {snakemake.input.fasta} "
	"{gtf_command} "
	"--runThreadN {snakemake.threads} "
	"{snakemake.params.star_usr} "
	"{log}")