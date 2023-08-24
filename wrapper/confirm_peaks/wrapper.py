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
    # loop through intervals
    for i in range(0,len(peak_ints)):
        rList = []
        s_int = int(peak_ints[i][1])
        e_int = int(peak_ints[i][2])
        peak_list[i] = peak_ints[i][0] + '_' + str(s_int) + '_' + str(e_int)
        for p2_rds in bam_content.fetch(peak_ints[i][0].tolist(), max(0,s_int-2000), e_int+2000):
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
    df_count_mat = pd.DataFrame({'name': peak_list, 'count': count_vec})
    df_count_mat[['chr', 'start', 'end']] = df_count_mat.name.str.split("_", expand = True)
    return(df_count_mat[['name', 'chr', 'start', 'end', 'count']])

def write_empty_files_and_exit(list_files_to_write, libsize_ip, libsize_sm):
    for filename in list_files_to_write:   
        with open(filename, 'w') as out:
            print("Unsufficient read support!")
            print("Reads in IP: " + str(libsize_ip))
            print("Reads in SM: " + str(libsize_sm))
    sys.exit(0)


    # with open(output_file_raw, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))        
    # with open(output_file_filtered, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))
    # with open(output_file_filtered_fc1, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))
    # with open(output_file_filtered_fc2, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))
    # with open(output_file_narrow, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))
    # with open(output_file_narrow_fc1, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))
    # with open(output_file_narrow_fc2, 'w') as out:
    #     print("Unsufficient read support!")
    #     print("Reads in IP: " + str(libsize_ip))
    #     print("Reads in SM: " + str(libsize_sm))

# ip_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_CTRL2/IP-cRIP_fourth_NSP9_8h_CTRL2_fwd.bam"
# sm_bam = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_fourth_NSP9_8h_CTRL2/SM-cRIP_fourth_NSP9_8h_CTRL2_fwd.bam"
# ip_bed = "results/peaks_coverage/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h_CTRL2-IP-SM/IP-cRIP_fourth_NSP9_8h_CTRL2_fwd.bedgraph"
# sm_bed = "results/peaks_coverage/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h_CTRL2-IP-SM/SM-cRIP_fourth_NSP9_8h_CTRL2_fwd.bedgraph"
# peak_file = "results/callpeaks_stranded/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h_CTRL2-IP-SM/fwd_peaks.narrowPeak"

# ip_bam = "results/dedup/only_sars_cov2_woFilter/IP-eCLIP_N_rep1/IP-eCLIP_N_rep1_fwd.bam"
# sm_bam = "results/dedup/only_sars_cov2_woFilter/SM-eCLIP_N_rep1/SM-eCLIP_N_rep1_fwd.bam"
# ip_bed = "results/peaks_coverage/only_sars_cov2_woFilter/eCLIP_N_rep1-IP-SM/IP-eCLIP_N_rep1_fwd.bedgraph"
# sm_bed = "results/peaks_coverage/only_sars_cov2_woFilter/eCLIP_N_rep1-IP-SM/SM-eCLIP_N_rep1_fwd.bedgraph"
# peak_file = "results/callpeaks_stranded/only_sars_cov2_woFilter/eCLIP_N_rep1-IP-SM/fwd_peaks.narrowPeak"


# ip_bam = "results/dedup/only_h1n1_woFilter/IP-eCLIP_NP_SZ/IP-eCLIP_NP_SZ_rev.bam"
# sm_bam = "results/dedup/only_h1n1_woFilter/SM-eCLIP_NP_SZ/SM-eCLIP_NP_SZ_rev.bam"
# ip_bed = "results/peaks_coverage/only_h1n1_woFilter/eCLIP_NP_SZ-IP-SM/IP-eCLIP_NP_SZ_rev.bedgraph"
# sm_bed = "results/peaks_coverage/only_h1n1_woFilter/eCLIP_NP_SZ-IP-SM/SM-eCLIP_NP_SZ_rev.bedgraph"
# peak_file = "results/callpeaks_stranded/only_h1n1_woFilter/eCLIP_NP_SZ-IP-SM/rev_peaks.narrowPeak"

# peaks_ip_bed = np.loadtxt(ip_bed, delimiter="\t", dtype=bytes).astype(str)
# peaks_ip = pd.DataFrame(peaks_ip_bed)
# peaks_ip.columns = ['chr', 'start', 'end', 'count']
# peaks_ip = pd.concat([pd.DataFrame({'name': peaks_ip['chr'] + "_" + peaks_ip['start'] + "_" + peaks_ip['end']}), peaks_ip], axis = 1)
    
# peaks_sm_bed = np.loadtxt(sm_bed, delimiter="\t", dtype=bytes).astype(str)
# peaks_sm = pd.DataFrame(peaks_sm_bed)
# peaks_sm.columns = ['chr', 'start', 'end', 'count']
# peaks_sm = pd.concat([pd.DataFrame({'name': peaks_sm['chr'] + "_" + peaks_sm['start'] + "_" + peaks_sm['end']}), peaks_sm], axis = 1)


peak_file = snakemake.input.peaks
ip_bam = snakemake.input.ip
sm_bam = snakemake.input.sm

use_read_cutoff_and_odds_ratio = snakemake.params.count_threshold
count_threshold = snakemake.params.count_threshold


peak_ints = np.loadtxt(peak_file, delimiter="\t", dtype=bytes).astype(str)

if snakemake.input.get("ip_bed", ""):
    peaks_ip_bed = np.loadtxt(snakemake.input.ip_bed, delimiter="\t", dtype=bytes).astype(str)
    # peak_names = ["_".join(list(peak_coord[:3])) for peak_coord in peaks_ip_bed]
    # peak_covs = [int(peak_coord[3]) for peak_coord in peaks_ip_bed]
    # peaks_ip = pd.DataFrame(peak_covs, index = peak_names)
    peaks_ip = pd.DataFrame(peaks_ip_bed)
    peaks_ip.columns = ['chr', 'start', 'end', 'count']
    peaks_ip = pd.concat([pd.DataFrame({'name': peaks_ip['chr'] + "_" + peaks_ip['start'] + "_" + peaks_ip['end']}), peaks_ip], axis = 1)
    
else:
    peaks_ip = count_peak_reads(ip_bam, peak_ints)


if snakemake.input.get("sm_bed", ""):  
    peaks_sm_bed = np.loadtxt(snakemake.input.sm_bed, delimiter="\t", dtype=bytes).astype(str)
    # peak_names = ["_".join(list(peak_coord[:3])) for peak_coord in peaks_sm_bed]
    # peak_covs = [int(peak_coord[3]) for peak_coord in peaks_sm_bed]
    # peaks_sm = pd.DataFrame(peak_covs, index = peak_names)
    peaks_sm = pd.DataFrame(peaks_sm_bed)
    peaks_sm.columns = ['chr', 'start', 'end', 'count']
    peaks_sm = pd.concat([pd.DataFrame({'name': peaks_sm['chr'] + "_" + peaks_sm['start'] + "_" + peaks_sm['end']}), peaks_sm], axis = 1)
else:
    peaks_sm = count_peak_reads(sm_bam, peak_ints)

peaks_ip = peaks_ip.astype({'count':'int'})
peaks_sm = peaks_sm.astype({'count':'int'})

output_file_raw = snakemake.output.raw
output_file_filtered = snakemake.output.filtered
output_file_filtered_fc1 = snakemake.output.filtered_fc1
output_file_filtered_fc2 = snakemake.output.filtered_fc2
output_file_narrow = snakemake.output.narrow
output_file_narrow_fc1 = snakemake.output.narrow_fc1
output_file_narrow_fc2 = snakemake.output.narrow_fc2

libsize_ip = paired_read_number(ip_bam)
libsize_sm = paired_read_number(sm_bam)

df_ip_sm = peaks_ip.merge(peaks_sm, on = ['name', 'chr', 'start', 'end'], suffixes = ['.IP', '.SM'])

if len(df_ip_sm) == 0:
    write_empty_files_and_exit([output_file_raw, output_file_filtered, output_file_filtered_fc1, 
                                output_file_filtered_fc2, output_file_narrow, output_file_narrow_fc1, 
                                output_file_narrow_fc2], libsize_ip, libsize_sm)

no_peak_ip = libsize_ip - df_ip_sm['count.IP']
no_peak_sm = libsize_sm - df_ip_sm['count.SM']

no_peak_df = pd.DataFrame({'noPeak.IP': no_peak_ip, 
                           'noPeak.SM': no_peak_sm})

df_ip_sm = pd.concat([df_ip_sm, no_peak_df], axis = 1)
df_ip_sm['log2FC(IP/SM)'] = np.log2((df_ip_sm['count.IP'] + 1)/(libsize_ip + len(df_ip_sm))) - np.log2((df_ip_sm['count.SM'] + 1)/(libsize_sm + len(df_ip_sm)))


odds_list = np.zeros((len(df_ip_sm)))
pval_list = np.zeros((len(df_ip_sm)))
for i in range(0, len(df_ip_sm)):   
    odds_list[i], pval_list[i] = fisher_exact([[df_ip_sm['count.IP'][i], df_ip_sm['noPeak.IP'][i]],
                                               [df_ip_sm['count.SM'][i], df_ip_sm['noPeak.SM'][i]]], alternative = 'greater')

fisher_df = pd.DataFrame({'oddsRatio': odds_list, 'p.value': pval_list})

df_ip_sm = pd.concat([df_ip_sm, fisher_df], axis = 1)

if len(df_ip_sm) == 0:
    df_ip_sm.to_csv(output_file_raw, sep = "\t")
    write_empty_files_and_exit([output_file_filtered, output_file_filtered_fc1, 
                                output_file_filtered_fc2, output_file_narrow, output_file_narrow_fc1, 
                                output_file_narrow_fc2], libsize_ip, libsize_sm)

df_ip_sm.to_csv(output_file_raw, sep = "\t")



if use_read_cutoff_and_odds_ratio:
    df_ip_sm = df_ip_sm[df_ip_sm['oddsRatio'] > 1].reset_index(drop = True)

if len(df_ip_sm) == 0:
    df_ip_sm.to_csv(output_file_raw, sep = "\t")
    write_empty_files_and_exit([output_file_filtered, output_file_filtered_fc1, 
                                output_file_filtered_fc2, output_file_narrow, output_file_narrow_fc1, 
                                output_file_narrow_fc2], libsize_ip, libsize_sm)


df_ip_sm['p.adj'] = multipletests(list(df_ip_sm['p.value']), method = 'fdr_by')[1]
df_ip_sm.to_csv(output_file_raw, sep = "\t")

df_ip_sm_filtered = df_ip_sm[df_ip_sm['p.adj'] < 0.05]
df_ip_sm_filtered = df_ip_sm_filtered.reset_index(drop = True)


df_ip_sm_filtered = df_ip_sm_filtered.assign(score = -1 * np.log10(df_ip_sm_filtered.loc[:,'p.adj']))
df_ip_sm_filtered.loc[np.isinf(df_ip_sm_filtered.loc[:,'score']), 'score'] = 400
df_ip_sm_filtered = df_ip_sm_filtered.assign(logpVal = -1 * np.log10(df_ip_sm_filtered.loc[:,'p.value']))
df_ip_sm_filtered.loc[np.isinf(df_ip_sm_filtered.loc[:,'logpVal']), 'logpVal'] = 400

df_ip_sm_filtered.to_csv(output_file_filtered, sep = "\t")

# narrow_dict = {'chr': list(ints_df[0]), 
#                'start': list(ints_df[1]), 
#                'end': list(ints_df[2]), 
#                'name': list(df_ip_sm_filtered.loc[:,'peak_name']), 
#                'score': list(df_ip_sm_filtered.loc[:,'score'] * 10), 
#                'strand': ['.'] * len(ints_df), 
#                'pValue': list(df_ip_sm_filtered.loc[:,'logpVal']), 
#                'qValue': list(df_ip_sm_filtered.loc[:,'score']), 
#                'peak': [1] * len(ints_df)}

# narrow_peak_df = pd.DataFrame(narrow_dict)

df_narrow = df_ip_sm_filtered[['chr', 'start', 'end', 'name', 'score', 'p.value', 'p.adj']]
# df_narrow['peak_score'] = df_narrow['score'] * 10
df_narrow_spec = pd.DataFrame({'peak_score': df_narrow['score'] * 10,
                               'strand': pd.Series(["." for x in range(len(df_narrow))]),
                               'peak': pd.Series([1 for x in range(len(df_narrow))])})
df_narrow = pd.concat([df_narrow, df_narrow_spec], axis = 1)

df_narrow = df_narrow[['chr', 'start', 'end', 'name', 'peak_score', 'strand', 'p.value', 'p.adj', 'peak']]
df_narrow.to_csv(output_file_narrow, sep = "\t", header = False, index = False)



df_ip_sm_filtered_fc1 = df_ip_sm_filtered.loc[np.abs(df_ip_sm_filtered['log2FC(IP/SM)']) > 1]
df_ip_sm_filtered_fc1 = df_ip_sm_filtered_fc1.reset_index(drop = True)
df_ip_sm_filtered_fc2 = df_ip_sm_filtered.loc[np.abs(df_ip_sm_filtered['log2FC(IP/SM)']) > 2]
df_ip_sm_filtered_fc2 = df_ip_sm_filtered_fc2.reset_index(drop = True)

df_ip_sm_filtered_fc1.to_csv(output_file_filtered_fc1, sep = "\t")
df_ip_sm_filtered_fc2.to_csv(output_file_filtered_fc2, sep = "\t")

df_narrow_fc1 = df_ip_sm_filtered_fc1[['chr', 'start', 'end', 'name', 'score', 'p.value', 'p.adj']]
# df_narrow['peak_score'] = df_narrow['score'] * 10
df_narrow_spec_fc1 = pd.DataFrame({'peak_score': df_narrow_fc1['score'] * 10,
                               'strand': pd.Series(["." for x in range(len(df_narrow_fc1))]),
                               'peak': pd.Series([1 for x in range(len(df_narrow_fc1))])})
df_narrow_fc1 = pd.concat([df_narrow_fc1, df_narrow_spec_fc1], axis = 1)

df_narrow_fc1 = df_narrow_fc1[['chr', 'start', 'end', 'name', 'peak_score', 'strand', 'p.value', 'p.adj', 'peak']]
df_narrow_fc1.to_csv(output_file_narrow_fc1, sep = "\t", header = False, index = False)

df_narrow_fc2 = df_ip_sm_filtered_fc1[['chr', 'start', 'end', 'name', 'score', 'p.value', 'p.adj']]
# df_narrow['peak_score'] = df_narrow['score'] * 10
df_narrow_spec_fc2 = pd.DataFrame({'peak_score': df_narrow_fc2['score'] * 10,
                               'strand': pd.Series(["." for x in range(len(df_narrow_fc2))]),
                               'peak': pd.Series([1 for x in range(len(df_narrow_fc2))])})
df_narrow_fc2 = pd.concat([df_narrow_fc2, df_narrow_spec_fc2], axis = 1)

df_narrow_fc2 = df_narrow_fc2[['chr', 'start', 'end', 'name', 'peak_score', 'strand', 'p.value', 'p.adj', 'peak']]
df_narrow_fc2.to_csv(output_file_narrow_fc2, sep = "\t", header = False, index = False)

