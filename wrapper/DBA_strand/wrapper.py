__author__ = "Alexander Gabel"
__copyright__ = "Copyright 2022, Alexander Gabel"
__email__ = "alexander.gabel@helmholtz-hzi.de."
__license__ = "MIT"

import os
import numpy as np
import pysam
import pandas as pd
import csv
import re
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests

def paired_read_number(bam_file):
    read_stats = pysam.flagstat(bam_file)
    [line.split("+") for line in read_stats.split("\n")]# 12.Element
    proper_pairs = read_stats.split("\n")
    number_prob_paired = int([line.split("+") for line in read_stats.split("\n")][11][0])
    return(number_prob_paired * 2)

def count_peak_reads(bamFile, peak_ints):
    # initialize vars
    count_vec = np.zeros((len(peak_ints)))
    bam_content = pysam.Samfile(bamFile, "rb")
    peak_list = [""] * len(peak_ints)
    # loop over intervals
    for i in range(0,len(peak_ints)):
        rList = []
        s_int = int(peak_ints.loc[i,1])
        e_int = int(peak_ints.loc[i,2])
        peak_list[i] = 'peak_' + peak_ints.loc[i,0] + '_' + str(s_int) + '_' + str(e_int)
        for p2_rds in bam_content.fetch(peak_ints.loc[i,0], max(0,s_int-2000), e_int+2000):
            if p2_rds.qname in rList:
                continue
            if p2_rds.is_reverse:
                insPos = p2_rds.pos + p2_rds.alen
            else:
                insPos = p2_rds.pos
            # if within window
            if insPos > s_int and insPos < e_int:
                count_vec[i] += 1
                rList.append(p2_rds.qname)
    # eof for
    df_count_mat = pd.DataFrame(count_vec,
                                index = peak_list)
    return(df_count_mat)

# peak_file = "results/confirmed_peaks/only_sars_cov2_woFilter/eCLIP_cnbp-IP-SM/eCLIP_cnbp_IP_SM_rev_signif.bed"
# fwd_bam = "results/dedup/only_sars_cov2_woFilter/IP-eCLIP_cnbp/IP-eCLIP_cnbp_fwd.bam"
# rev_bam = "results/dedup/only_sars_cov2_woFilter/IP-eCLIP_cnbp/IP-eCLIP_cnbp_rev.bam"

# strand = "rev"
# count_threshold = 10

# peak_file = "../results/confirmed_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_12h_CTRL2-IP-SM/cRIP_fourth_NSP9_12h_CTRL2_IP_SM_fwd_signif.bed" 
# fwd_bam = "../results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_12h_CTRL2/IP-cRIP_fourth_NSP9_12h_CTRL2_fwd.bam" 
# rev_bam = "../results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_12h_CTRL2/IP-cRIP_fourth_NSP9_12h_CTRL2_rev.bam" 
# strand = "fwd"
# count_threshold = 10

peak_file = snakemake.input.peaks
fwd_bam = snakemake.input.fwd
rev_bam = snakemake.input.rev

strand = snakemake.params.strand
count_threshold = snakemake.params.count_threshold

output_file_raw = snakemake.output.raw
output_file_filtered = snakemake.output.filtered

try:
    peak_ints = pd.read_csv(peak_file, delimiter = "\t", header = None, index_col = None)
except:
    print("Peak file is empty.")
    with open(output_file_raw, 'w') as out:
        print("Could not find peaks for comparison!")
    with open(output_file_filtered, 'w') as out:
        print("Could not find peaks for comparison!")
    exit(0)

libsize_fwd = paired_read_number(fwd_bam)
libsize_rev = paired_read_number(rev_bam)

peaks_fwd = count_peak_reads(fwd_bam, peak_ints)
peaks_rev = count_peak_reads(rev_bam, peak_ints)

res_df = pd.DataFrame(data = {"peak_name": list(peaks_fwd.index), 
                              "count_fwd": list(peaks_fwd[0]), 
                              "count_rev": list(peaks_rev[0]),
                              "libSize.fwd": libsize_fwd, 
                              "libSize.rev": libsize_rev})


res_df = res_df.loc[(res_df['count_fwd'] > count_threshold) & (res_df['count_rev'] > count_threshold)]

res_df['Peak_prob'] = (res_df.count_fwd + res_df.count_rev)/(libsize_fwd + libsize_rev)

if strand == 'fwd':
    res_df['lambda'] = res_df.count_rev

if strand == 'rev':
    res_df['lambda'] = res_df.count_fwd

res_df['CPM_fwd'] = (res_df.count_fwd)/(libsize_fwd + libsize_rev) * 1e6
res_df['CPM_rev'] = (res_df.count_rev)/(libsize_fwd + libsize_rev) * 1e6

if strand == 'fwd':
    res_df['log2FC(fwd/rev)'] = np.log2(res_df['CPM_fwd'] + 1)  - np.log2(res_df['CPM_rev'] + 1) 

if strand == 'rev':
    res_df['log2FC(rev/fwd)'] =  np.log2(res_df['CPM_rev'] + 1)  - np.log2(res_df['CPM_fwd'] + 1)

res_df['fraction_fwd'] = (res_df.count_fwd)/(res_df.count_fwd + res_df.count_rev)
res_df['fraction_rev'] = (res_df.count_rev)/(res_df.count_fwd + res_df.count_rev)

if strand == 'fwd':
    res_df['p.value'] = poisson.pmf(k = res_df.count_fwd, mu = res_df.count_rev)

if strand == 'rev':
    res_df['p.value'] = poisson.pmf(k = res_df.count_rev, mu = res_df.count_fwd)

res_df['p.adj(BY)'] = multipletests(list(res_df['p.value']), method = 'fdr_by')[1]

res_df_filtered = res_df[res_df['p.adj(BY)'] < 0.05]

res_df.to_csv(output_file_raw, sep = "\t")
res_df_filtered.to_csv(output_file_filtered, sep = "\t")
