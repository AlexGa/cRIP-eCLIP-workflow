__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
import tempfile
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

tmpdir = snakemake.params.get("tmpdir")

fq1 = snakemake.input.get("fq1")
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = (
    [snakemake.input.fq1]
    if isinstance(snakemake.input.fq1, str)
    else snakemake.input.fq1
)
fq2 = snakemake.input.get("fq2")
if fq2:
    fq2 = (
        [snakemake.input.fq2]
        if isinstance(snakemake.input.fq2, str)
        else snakemake.input.fq2
    )
    assert len(fq1) == len(
        fq2
    ), "input-> equal number of files required for fq1 and fq2"
input_str_fq1 = ",".join(fq1)
input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
input_str = " ".join([input_str_fq1, input_str_fq2])

if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

index = snakemake.params.get("idx")
# if not index:
#     index = snakemake.params.get("idx", "")

# with tempfile.TemporaryDirectory() as tmpdir:

shell(
    "STAR "
    " --runThreadN {snakemake.threads}"
    " --genomeDir {index}"
    " --readFilesIn {input_str}"
    " {readcmd}"
    " {extra}"
    " --outFileNamePrefix {tmpdir}/"
    " --outStd Log "
    " {log}"
)

if "SortedByCoordinate" in extra:
    bamprefix = "Aligned.sortedByCoord.out"
else:
    bamprefix = "Aligned.out"

if snakemake.output.get("bam"):
    if tmpdir + "/" + bamprefix + ".bam" != snakemake.output.get("bam"):  
        shell("mv {tmpdir}/{bamprefix}.bam {snakemake.output.bam:q}")

if snakemake.output.get("sam"):
    if tmpdir + "/" + bamprefix + ".sam" != snakemake.output.get("sam"):  
        shell("mv {tmpdir}/{bamprefix}.sam {snakemake.output.sam:q}")

if snakemake.output.get("reads_per_gene"):
    if tmpdir + "/ReadsPerGene.out.tab" != snakemake.output.get("reads_per_gene"):  
        shell("mv {tmpdir}/ReadsPerGene.out.tab {snakemake.output.reads_per_gene:q}")

if snakemake.output.get("chim_junc"):
    if tmpdir + "/Chimeric.out.junction" != snakemake.output.get("chim_junc"):  
        shell("mv {tmpdir}/Chimeric.out.junction {snakemake.output.chim_junc:q}")

if snakemake.output.get("sj"):
    if tmpdir + "/SJ.out.tab" != snakemake.output.get("sj"):  
        shell("mv {tmpdir}/SJ.out.tab {snakemake.output.sj:q}")

if snakemake.output.get("log"):
    if tmpdir + "/Log.out" != snakemake.output.get("log"):  
        shell("mv {tmpdir}/Log.out {snakemake.output.log:q}")

if snakemake.output.get("log_progress"):
    if tmpdir + "/Log.progress.out" != snakemake.output.get("log_progress"): 
        shell("mv {tmpdir}/Log.progress.out {snakemake.output.log_progress:q}")

if snakemake.output.get("log_final"):
    if tmpdir + "/Log.final.out" != snakemake.output.get("log_final"):
        shell("mv {tmpdir}/Log.final.out {snakemake.output.log_final:q}")
