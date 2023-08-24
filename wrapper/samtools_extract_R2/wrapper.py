__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

shell("samtools view -hb -f 130 -@ {snakemake.threads} {snakemake.input[0]} -o {snakemake.output[0]}")
shell("samtools index {snakemake.output[0]}")
