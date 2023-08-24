__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts
import os 
import shutil

samtools_opts = get_samtools_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

fwd_file = snakemake.output.get("fwd")
rev_file = snakemake.output.get("rev")
paired_file = snakemake.output.get("paired")

tmp_dir = snakemake.params.get("tmpdir")

if not os.path.exists(tmp_dir):
	os.makedirs(tmp_dir)

def extract_strands(bam_file, suffix, flag_array, tmp_dir, out_file):
	bam_filename = os.path.basename(bam_file)
	header_file = os.path.join(tmp_dir, bam_filename.replace(".bam", "_header.sam"))
	strand_tmp_file = os.path.join(tmp_dir, bam_filename.replace(".bam", "_" + suffix))
	''' Extract header from bam file '''
	shell("samtools view -H {bam_file} > {header_file}")
	''' Split bam file according to array of flags'''
	strand_files = dict([[i, strand_tmp_file + "_" + str(i) +".sam"] for i in flag_array])
	for flag in strand_files:
		out_strand_file = strand_files[flag]
		shell("samtools view -f {flag} {bam_file} > {out_strand_file}")
	''' Merge, sort and index flag specific file'''
	joined_files = " ".join([header_file] + list(strand_files.values()))
	shell("cat {joined_files} > {strand_tmp_file}_tmp.sam")
	shell("samtools sort -O BAM -o {out_file} {strand_tmp_file}_tmp.sam")
	shell("samtools index {out_file}")

extract_strands(snakemake.input[0], "fwd", [163, 83], tmp_dir, fwd_file)
extract_strands(snakemake.input[0], "rev", [147, 99], tmp_dir, rev_file)
extract_strands(snakemake.input[0], "paired", [163, 83, 147, 99], tmp_dir, paired_file)

try:
	shutil.rmtree(tmp_dir)
except OSError as error:
	print("Error: %s : %s:" %(tmp_dir, error.strerror))