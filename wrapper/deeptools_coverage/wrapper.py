__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)


extra = snakemake.params.extra

strand = snakemake.params.get("strand", "")

if isinstance(strand, str):
    strand = strand.replace("_", "")
    if strand in ["forward", "reverse"]:
        extra = extra + " --filterRNAstrand " + strand


shell("bamCoverage "
      "{extra} "
      "-b {snakemake.input.bam} "
      "-o {snakemake.output.bigwig} "
      "-of bigwig {log}")