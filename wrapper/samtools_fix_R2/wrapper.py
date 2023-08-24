__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

shell("samtools view -b -S <(cat <(samtools view -H {snakemake.input[0]}) <(samtools view {snakemake.input[0]} | awk '$2 = ($2 == \"163\" ? \"0\" : \"16\")' | tr \" \" \"\\t\")) > {snakemake.output[0]}")