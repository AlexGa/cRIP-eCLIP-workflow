__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
import os

samtools_opts = get_samtools_opts(snakemake)
read_type = snakemake.params.get("read_type", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


bai_file = snakemake.output[1].replace(".bam", ".bai")

shell("samtools sort --threads {snakemake.threads} -o {snakemake.output[0]} {snakemake.input[0]}")

if isinstance(region, str):

	sorted_file_name = snakemake.output[0].replace(".bam", ".bai")
	header_file = snakemake.output[0].replace(".bam", "_header.sam")
	temp_sam = snakemake.output[1].replace(".bam", ".sam")

	# Index prep for region extraction
	shell("samtools index {snakemake.output[0]} {sorted_file_name}")

	# Extract region of interest
	shell("samtools view --threads {snakemake.threads} -o {snakemake.output[1]}_temp.sam --output-fmt SAM  {snakemake.output[0]} {region}")
	
	# Extract header and only region of interest as SQ  
	shell("samtools view -H {snakemake.output[0]} | grep -E \"@(HD|PG|CO)|{region}\" > {header_file}")

	# Concat region specific header and region specific sam
	shell("cat {header_file} {snakemake.output[1]}_temp.sam > {temp_sam}")

	# Create bam file out of region sam
	shell("samtools view --threads {snakemake.threads} --write-index -o {snakemake.output[1]}##idx##{bai_file} --output-fmt BAM  {temp_sam}")

	# Clean up
	os.unlink(sorted_file_name)
	os.unlink(header_file)
	os.unlink(temp_sam)
	os.unlink(snakemake.output[1]+"_temp.sam")
	
else:
	shell("samtools view --threads {snakemake.threads} --write-index -o  {snakemake.output[1]}##idx##{bai_file} --output-fmt BAM {snakemake.output[0]}")


def is_paired_sequencing(bamfile):
    # TODO: this is scary. Should check for unpaired being 0, and number paired == total number
    r = pysam.flagstat(bamfile)
    paired = int(r[5].split()[0])
    if paired != 0:
        return True
    else:
        return False

