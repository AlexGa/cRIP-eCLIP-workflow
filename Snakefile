from snakemake.utils import min_version
import yaml
import pandas as pd
from yaml.loader import SafeLoader
##### set minimum snakemake version #####
min_version("7.8.0")


##### setup report #####
configfile: "config/config.yaml"


report: "eclip-workflow/report/workflow.rst"


units = pd.read_csv(config["units"], sep="\t", dtype={"sample": str, "unit": str}).set_index(["sample", "unit"], drop=False).sort_index()
peaks = pd.read_csv(config["peak_comparison"], sep="\t", dtype={"unit": str}).set_index(["unit"], drop=False).sort_index()

with open(config['sample_comparison']) as fh:
  sample_comp = yaml.load(fh, SafeLoader)

##### load rules #####

include: "eclip-workflow/rules/trim.smk"
include: "eclip-workflow/rules/alignment.smk"
include: "eclip-workflow/rules/peak_calling.smk"
include: "eclip-workflow/rules/crosslinking_sites.smk"

rule all:
  input:
    expand("results/clipTrim/{file.sample}-{file.unit}/{file.sample}_{file.unit}_R1.fq.gz", file = units.itertuples())
    ,expand("results/trimmed/round1/{file.sample}-{file.unit}/{file.sample}_{file.unit}_R1.fq.gz", file = units.itertuples())
    ,expand("results/DBA_strand/sars_cov2/{peak.unit}/{peak.unit}_fwd_peaks_raw.tsv", peak = peaks.itertuples())
    ,expand("results/DBA_strand/sars_cov2/{peak.unit}/{peak.unit}_rev_peaks_raw.tsv", peak = peaks.itertuples())
    ,expand("results/DBA/sars_cov2/{unit}/{unit}_fwd_raw.tsv", unit = sample_comp['DBA'].keys())
    ,expand("results/DBA/sars_cov2/{unit}/{unit}_rev_raw.tsv", unit = sample_comp['DBA'].keys())
    ,expand("results/normed_tracks_peaks/sars_cov2/KL/{peak.unit}/{peak.unit}_KL_fwd.bw", peak = peaks.itertuples())
    ,expand("results/normed_tracks_peaks/sars_cov2/KL/{peak.unit}/{peak.unit}_KL_rev.bw", peak = peaks.itertuples())
    ,expand("results/DBA_XL_sites/sars_cov2/{unit}/{unit}_fwd_raw.tsv", unit = sample_comp['DBA'].keys())
    ,expand("results/DBA_XL_sites/sars_cov2/{unit}/{unit}_rev_raw.tsv", unit = sample_comp['DBA'].keys())
    ,expand("results/normed_tracks_XL_sites/single/sars_cov2/KL/{peak.unit}/{peak.unit}_KL_fwd.bw", peak = peaks.itertuples())
    ,expand("results/normed_tracks_XL_sites/single/sars_cov2/KL/{peak.unit}/{peak.unit}_KL_rev.bw", peak = peaks.itertuples())
