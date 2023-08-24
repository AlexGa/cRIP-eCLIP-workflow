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
from scipy.stats import fisher_exact
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

# peak_file = "results/peak_intersect/only_sars_cov2_woFilter/cRIP_fourth_NSP9_12h-IP-SM_fwd_peaks.bed" 
# ctrl_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_12h_CTRL2/IP-cRIP_fourth_NSP9_12h_CTRL2_fwd.bam" 
# ko_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_12h_KO6_plus_eV/IP-cRIP_fourth_NSP9_12h_KO6_plus_eV_fwd.bam"

peak_file = snakemake.input.peaks
ctrl_bam = snakemake.input.ctrl
ko_bam = snakemake.input.ko

output_file_raw = snakemake.output.raw
output_file_filtered = snakemake.output.filtered
output_file_narrow = snakemake.output.bed

peak_ints = pd.read_csv(peak_file, delimiter = "\t", header = None)
 #np.loadtxt(peak_file, delimiter="\t", dtype=bytes).astype(str)
libsize_ctrl = paired_read_number(ctrl_bam)
libsize_ko = paired_read_number(ko_bam)

peaks_ctrl = count_peak_reads(ctrl_bam, peak_ints)
peaks_ko = count_peak_reads(ko_bam, peak_ints)

odds_list = np.zeros((len(peak_ints)))
pval_list = np.zeros((len(peak_ints)))
no_peak_ctrl = libsize_ctrl - peaks_ctrl
no_peak_ko = libsize_ko - peaks_ko
fc_list = np.zeros((len(peak_ints)))

for i in range(0, len(peaks_ctrl)):   
    odds_list[i], pval_list[i] = fisher_exact([[peaks_ko[0][i], no_peak_ko[0][i]],[peaks_ctrl[0][i], no_peak_ctrl[0][i]]])


norm_ko = (peaks_ko + 1)/libsize_ko * 1e6
norm_ctrl = (peaks_ctrl + 1)/libsize_ctrl * 1e6

fc_list = np.log2(norm_ko) - np.log2(norm_ctrl)

padj_list = multipletests(list(pval_list), method = 'fdr_by')

d = {'peak_name': list(peaks_ctrl.index), 'count.ctrl': list(peaks_ctrl[0]), 'count.ko': list(peaks_ko[0]),
     'noPeak_ctrl':  list(no_peak_ctrl[0]), 'noPeak_ko':  list(no_peak_ko[0]), 
     'CPM.ctrl': list(norm_ctrl[0]), 'CPM.ko': list(norm_ko[0]), 'log2FC(ko/ctrl)': list(fc_list[0]), 
     'oddsRatio': list(odds_list), 'p.value': list(pval_list), 'p.adj': list(padj_list[1])}

df_ctrl_ko = pd.DataFrame(data=d)

df_ctrl_ko_filtered = df_ctrl_ko.loc[padj_list[1] < 0.05]
ints_df = pd.DataFrame(peak_ints).loc[padj_list[1] < 0.05]

df_ctrl_ko_filtered = df_ctrl_ko_filtered.assign(score = -1 * np.log10(df_ctrl_ko_filtered.loc[:,'p.adj']))
df_ctrl_ko_filtered.loc[np.isinf(df_ctrl_ko_filtered.loc[:,'score']), 'score'] = 400
df_ctrl_ko_filtered = df_ctrl_ko_filtered.assign(logpVal = -1 * np.log10(df_ctrl_ko_filtered.loc[:,'p.value']))
df_ctrl_ko_filtered.loc[np.isinf(df_ctrl_ko_filtered.loc[:,'logpVal']), 'logpVal'] = 400


narrow_dict = {'chr': list(ints_df[0]), 'start': list(ints_df[1]), 'end': list(ints_df[2]), 'name': list(df_ctrl_ko_filtered.loc[:,'peak_name']), 'score': list(df_ctrl_ko_filtered.loc[:,'score'] * 10), 'strand': ['.'] * len(ints_df), 'pValue': list(df_ctrl_ko_filtered.loc[:,'logpVal']), 'qValue': list(df_ctrl_ko_filtered.loc[:,'score']), 'peak': [1] * len(ints_df)}

narrow_peak_df = pd.DataFrame(narrow_dict)
narrow_peak_df.to_csv(output_file_narrow, sep = "\t", header = False, index = False)


df_ctrl_ko.to_csv(output_file_raw, sep = "\t")
df_ctrl_ko_filtered.to_csv(output_file_filtered, sep = "\t")
