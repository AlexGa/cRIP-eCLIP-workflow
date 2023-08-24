import pandas as pd
import pysam
import numpy as np
import re
import os
import sys
import pyBigWig
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def total_xl_reads(bam_file):
    read_stats = pysam.flagstat(bam_file)
    [line.split("+") for line in read_stats.split("\n")]# 12.Element
    proper_pairs = read_stats.split("\n")
    number_prob_paired = int([line.split("+") for line in read_stats.split("\n")][11][0])
    return(number_prob_paired/2)

def count_x_link_sites(bam_file, orientation = 'fwd', peak_ints_file = None):
    pairedreads = pysam.AlignmentFile(bam_file, "rb")
    chr_length_dict = dict(zip(*[list(pairedreads.references), list(pairedreads.lengths)]))
    pwm = {}
    count_matrix = {}
    for key in chr_length_dict.keys():
        pwm[key] = np.zeros((chr_length_dict[key] + 1, 5)) # A -> 0, C -> 1, G -> 2, T -> 3
        count_matrix[key] = np.zeros((1, chr_length_dict[key]))
    bp2int = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
    chr_length = 0
    chr_name = ""
    if peak_ints_file != None:
        peak_ints = pd.read_csv(peak_ints_file, delimiter = "\t", header = None)
    for read in pairedreads.fetch():
        chr_name = read.reference_name
        chr_length = chr_length_dict[chr_name]
        print(chr_name)
        print(count_matrix[chr_name])
        if peak_ints_file != None:
            if np.sum((read.reference_start <= peak_ints[2]) & (read.reference_end >= peak_ints[1])) == 0:
                continue
        if orientation == 'fwd' and read.flag == 163:
            if (read.reference_start - 1) >= 0:
                count_matrix[chr_name][0][read.reference_start-1] += 1 # start is inclusive
                pwm[chr_name][read.reference_start][bp2int[read.query_sequence[0]]] += 1
            else:
                count_matrix[chr_name][0][0] += 1
                pwm[chr_name][bp2int[read.query_sequence[0]]] += 1
        if orientation == 'rev' and read.flag == 147:
            if (read.reference_end) < chr_length:
                count_matrix[chr_name][0][read.reference_end] += 1 # start in count_matrix at 0 but later at 1 and bed file end is not inclusive => read.reference_end - 1 points at
                pwm[chr_name][read.reference_end - 1][bp2int[read.query_sequence[len(read.query_sequence)-1]]] += 1
            else:
                count_matrix[chr_name][0][chr_length-1] += 1
                pwm[chr_length][bp2int[read.query_sequence[len(read.query_sequence)-1]]] += 1
    # close file
    pairedreads.close()
    # create dataframe
    print("Create DataFrame...")
    bed_df = pd.DataFrame()
    for key in count_matrix.keys():
        bed_sub_df = pd.DataFrame({'chr': key,
                                   'start': range(0, chr_length_dict[key]),
                                   'end': range(1, chr_length_dict[key] + 1),
                                   'strand': '*',
                                   'score': pd.Series(list(count_matrix[key][0]))})
        bed_df = pd.concat([bed_df, bed_sub_df])
    bed_df = bed_df.set_index(['chr', 'start'])
    # bed_df.index = bed_df['start']
    xl_data = {"bed": bed_df, "pwm": pwm}
    return(xl_data)

def write_bw_from_bed(bw_df, output_file, chr_length_set, chr_length_dict):
    bw = pyBigWig.open(output_file, "w")
    bw.addHeader(chr_length_set, maxZooms=0)
    for chrom in chr_length_dict.keys():
        chroms = np.array(list(bw_df[bw_df['chr'] == chrom]['chr']))
        starts = np.array(list(bw_df[bw_df['chr'] == chrom]['start']), dtype=np.int64)
        ends = np.array(list(bw_df[bw_df['chr'] == chrom]['end']), dtype=np.int64)
        values0 = np.array(list(bw_df[bw_df['chr'] == chrom]['score']), dtype=np.float64)
        if len(values0) > 0:
            bw.addEntries(chroms, starts, ends=ends, values=values0)
    bw.close()

# bamFile_ip = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bam"
# bamFile_sm = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bam"

# bedFile = "results/confirmed_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h_CTRL2-IP-SM/cRIP_fifth_rep1_NSP9_12h_CTRL2_IP_SM_fwd_signif.narrowPeak"
# bamFile_ip = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV_fwd.bam"
# bamFile_sm = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_fourth_NSP9_8h_KO6_plus_eV/SM-cRIP_fourth_NSP9_8h_KO6_plus_eV_fwd.bam"

# bamFile_ip = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_8h_CTRL2/IP-cRIP_fifth_rep1_NSP9_8h_KO6_plus_SND1_rev.bam"
# bamFile_sm = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_8h_CTRL2/SM-cRIP_fifth_rep1_NSP9_8h_KO6_plus_SND1_rev.bam"

# bamFile_ip = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep2_NSP9_12h_KO6_plus_SND1/IP-cRIP_fifth_rep2_NSP9_12h_KO6_plus_SND1_rev.bam"
# bamFile_sm = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep2_NSP9_12h_KO6_plus_SND1/SM-cRIP_fifth_rep2_NSP9_12h_KO6_plus_SND1_rev.bam"

# orientation = "fwd"
# chr_length =  29903
# count_threshold = 10

# bamFile_ip = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd.bam"
# bamFile_sm = "/vol/projects/agabel/NSP9_new/results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV/SM-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd.bam"
# bed_ip = "/vol/projects/agabel/NSP9_new/results/XL_sites_bedgraph/only_sars_cov2_woFilter/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd.bedgraph"
# bed_sm = "/vol/projects/agabel/NSP9_new/results/XL_sites_bedgraph/only_sars_cov2_woFilter/SM-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV/SM-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd.bedgraph"

# orientation = "fwd"
# count_threshold = 10

# bamFile_ip = "/vol/projects/agabel/eCLIP_analysis/Influenza/results/dedup/only_h1n1_woFilter/IP-eCLIP_CNBP_SZ/IP-eCLIP_CNBP_SZ_fwd.bam"
# bamFile_sm = "/vol/projects/agabel/eCLIP_analysis/Influenza/results/dedup/only_h1n1_woFilter/SM-eCLIP_CNBP_SZ/SM-eCLIP_CNBP_SZ_fwd.bam"
# bed_ip = "/vol/projects/agabel/eCLIP_analysis/Influenza/results/XL_sites_bedgraph/only_h1n1_woFilter/IP-eCLIP_CNBP_SZ/IP-eCLIP_CNBP_SZ_fwd.bedgraph"
# bed_sm = "/vol/projects/agabel/eCLIP_analysis/Influenza/results/XL_sites_bedgraph/only_h1n1_woFilter/SM-eCLIP_CNBP_SZ/SM-eCLIP_CNBP_SZ_fwd.bedgraph"


# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_IP_SM_fwd_raw.tsv
# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_IP_SM_fwd_signif.tsv
# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_IP_SM_fwd_signif.bed
# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd_raw.bw
# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/SM-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd_raw.bw
# results/XL_sites_all/only_sars_cov2_woFilter/cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV-IP-SM/IP-cRIP_fifth_rep12_NSP9_12h_KO6_plus_eV_fwd_signif.bw

# bamFile_ip = "results/dedup/only_sars_cov2_woFilter/IP-eCLIP_N_rep1/IP-eCLIP_N_rep1_fwd.bam"
# bamFile_sm = "results/dedup/only_sars_cov2_woFilter/SM-eCLIP_N_rep1/SM-eCLIP_N_rep1_fwd.bam"
# bed_ip = "results/XL_sites_bedgraph/only_sars_cov2_woFilter/IP-eCLIP_N_rep1/IP-eCLIP_N_rep1_fwd.bedgraph"
# bed_sm = "results/XL_sites_bedgraph/only_sars_cov2_woFilter/SM-eCLIP_N_rep1/SM-eCLIP_N_rep1_fwd.bedgraph"
# bedFile = "results/confirmed_peaks/only_sars_cov2_woFilter/eCLIP_N_rep1-IP-SM/eCLIP_N_rep1_IP_SM_rev_signif.bed"

bamFile_ip = snakemake.input.ip
bamFile_sm = snakemake.input.sm

bed_ip = snakemake.input.ip_bed
bed_sm = snakemake.input.sm_bed

xl_results_ip = pd.read_csv(bed_ip, delimiter = "\t", header = None)
xl_results_ip.columns = ['chr', 'start', 'end', 'score']
xl_results_ip = xl_results_ip.set_index(['chr', 'start']) #.index = xl_results_ip['start']

try:
    xl_results_sm = pd.read_csv(bed_sm, delimiter = "\t", header = None)
    xl_results_sm.columns = ['chr', 'start', 'end', 'score']
    xl_results_sm = xl_results_sm.set_index(['chr', 'start']) #= xl_results_sm['start']
except:
    xl_results_sm = xl_results_ip.copy()
    xl_results_sm['score'] = 0

# check if bed file with significant peaks is given
# if so than only consider reads/XL-sites overlapping with significant 
# peaks
# bedFile = snakemake.input.get("bed", "")


orientation = snakemake.params.orientation
# chr_length =  snakemake.params.chr_length
count_threshold = snakemake.params.count_threshold

output_file_raw = snakemake.output.tsv_raw
output_file_filtered = snakemake.output.tsv_filtered
output_file_narrow = snakemake.output.bed

output_file_bw_ip = snakemake.output.bw_ip
output_file_bw_sm = snakemake.output.bw_sm

output_file_bw_filtered_ip = snakemake.output.bw_filtered

xl_reads_ip = total_xl_reads(bamFile_ip)
xl_reads_sm = total_xl_reads(bamFile_sm)


bamFile = pysam.AlignmentFile(bamFile_ip, "rb")
chr_length_set = list(zip(*[list(bamFile.references), list(bamFile.lengths)]))
chr_length_dict = dict(zip(*[list(bamFile.references), list(bamFile.lengths)]))

# if bedFile:
#     xl_results_ip = count_x_link_sites(bamFile_ip, chr_length = chr_length, orientation = orientation, peak_ints_file = bedFile)
#     xl_results_sm = count_x_link_sites(bamFile_sm, chr_length = chr_length, orientation = orientation, peak_ints_file = bedFile)
# else:
#     xl_results_ip = count_x_link_sites(bamFile_ip, chr_length = chr_length, orientation = orientation)
#     xl_results_sm = count_x_link_sites(bamFile_sm, chr_length = chr_length, orientation = orientation)

ip_cov = xl_results_ip['score']
sm_cov = xl_results_sm['score']

write_bw_from_bed(xl_results_ip.reset_index(), output_file_bw_ip, chr_length_set, chr_length_dict)
write_bw_from_bed(xl_results_sm.reset_index(), output_file_bw_sm, chr_length_set, chr_length_dict)

idx_suff_coverage = (ip_cov > count_threshold)

# If no reads left after filtering, exit the script
if len(ip_cov[idx_suff_coverage]) == 0:
    with open(output_file_raw, 'w') as out:
        print("Unsufficient read support!")
        print("Reads in IP: " + str(xl_reads_ip))
        print("Reads in SM: " + str(xl_reads_ip))
    with open(output_file_filtered, 'w') as out:
        print("Unsufficient read support!")
        print("Reads in IP: " + str(xl_reads_ip))
        print("Reads in SM: " + str(xl_reads_ip))
    with open(output_file_narrow, 'w') as out:
        print("Unsufficient read support!")
        print("Reads in IP: " + str(xl_reads_ip))
        print("Reads in SM: " + str(xl_reads_ip))
    with open(output_file_bw_filtered_ip, 'w') as out:
        print("Unsufficient read support!")
        print("Reads in IP: " + str(xl_reads_ip))
        print("Reads in SM: " + str(xl_reads_ip))
    sys.exit(0)

xl_results_ip = xl_results_ip[idx_suff_coverage] 

xl_results_ip_sm = xl_results_ip.merge(xl_results_sm, on=['chr', 'start'], how='left')

xl_results_sm = xl_results_ip_sm[['end_x', 'score_y']]
xl_results_sm['score_y'].fillna(0, inplace=True)
xl_results_sm.columns = xl_results_ip.columns
xl_results_sm = xl_results_sm.reset_index()
xl_results_ip = xl_results_ip.reset_index()

ip_cov = xl_results_ip['score']
sm_cov = xl_results_sm['score']

df_ip_sm = pd.DataFrame()
for chrom in chr_length_dict:
    odds_dict_chrom = pval_dict_chrom = fc_dict_chrom = np.full(np.sum(xl_results_ip['chr'] == chrom), np.nan) #np.zeros(np.sum(xl_results_ip['chr'] == chrom))
    ip_cov_chrom = list(xl_results_ip[xl_results_ip['chr'] == chrom]['score'])
    sm_cov_chrom = list(xl_results_sm[xl_results_sm['chr'] == chrom]['score'])
    no_ip_cov_chrom = list(xl_reads_ip - xl_results_ip[xl_results_ip['chr'] == chrom]['score'])
    no_sm_cov_chrom = list(xl_reads_sm - xl_results_sm[xl_results_sm['chr'] == chrom]['score'])
    for i in range(0, len(ip_cov_chrom)):   
        odds_dict_chrom[i], pval_dict_chrom[i] = fisher_exact([ [ip_cov_chrom[i], no_ip_cov_chrom[i]],
                                                                  [sm_cov_chrom[i], no_sm_cov_chrom[i]]], 
                                                                  alternative = 'greater')
    fc_dict_chrom = np.log2((np.array(ip_cov_chrom) + 1)/(xl_reads_ip + 1)) - np.log2((np.array(sm_cov_chrom) + 1)/(xl_reads_sm + 1))
    d = {'chr': list(xl_results_ip[xl_results_ip['chr'] == chrom]['chr']),
         'xl_site': list(xl_results_ip[xl_results_ip['chr'] == chrom]['start']), 
         'count.ip': list(ip_cov_chrom), 
         'count.sm': list(sm_cov_chrom), 
         'no_xl_ip':  list(no_ip_cov_chrom), 
         'no_xl_sm':  list(no_sm_cov_chrom), 
         'log2FC(ip/sm)': list(fc_dict_chrom), 
         'oddsRatio': list(odds_dict_chrom), 
         'p.value': list(pval_dict_chrom)}
    df_ip_sm_chrom = pd.DataFrame(d)
    df_ip_sm = pd.concat([df_ip_sm, df_ip_sm_chrom])

df_ip_sm = df_ip_sm.reset_index(drop = True)
padj_list = multipletests(list(df_ip_sm['p.value']), method = 'fdr_by')[1]

df_ip_sm = pd.concat([df_ip_sm, pd.DataFrame({'p.adj': padj_list})], axis = 1)

df_ip_sm.to_csv(output_file_raw, sep = "\t", index = False)

df_ip_sm_filtered = df_ip_sm[df_ip_sm['p.adj'] < 0.05]

if len(df_ip_sm_filtered) == 0:
    with open(output_file_filtered, 'w') as out:
        print("No significant XL sites!")
    with open(output_file_narrow, 'w') as out:
        print("No significant XL sites!")
    with open(output_file_bw_filtered_ip, 'w') as out:
        print("No significant XL sites!")
    sys.exit(0)

df_ip_sm_filtered = df_ip_sm_filtered.assign(score = -1 * np.log10(df_ip_sm_filtered.loc[:,'p.adj']))
df_ip_sm_filtered.loc[np.isinf(df_ip_sm_filtered.loc[:,'score']), 'score'] = 400
df_ip_sm_filtered = df_ip_sm_filtered.assign(logpVal = -1 * np.log10(df_ip_sm_filtered.loc[:,'p.value']))
df_ip_sm_filtered.loc[np.isinf(df_ip_sm_filtered.loc[:,'logpVal']), 'logpVal'] = 400

narrow_dict = {'chr': list(df_ip_sm_filtered['chr']), 
               'start': list(df_ip_sm_filtered['xl_site']), 
               'end': list(df_ip_sm_filtered['xl_site']), 
               'name': list("xl_site_" + df_ip_sm_filtered['chr'] + "_" + df_ip_sm_filtered['xl_site'].astype(str)), 
               'score': list(df_ip_sm_filtered['score'] * 10), 
               'strand': ['.'] * len(df_ip_sm_filtered), 
               'pValue': list(df_ip_sm_filtered['logpVal']), 
               'qValue': list(df_ip_sm_filtered['score']), 
               'peak': [1] * len(df_ip_sm_filtered)}


narrow_peak_df = pd.DataFrame(narrow_dict)
narrow_peak_df.to_csv(output_file_narrow, sep = "\t", header = False, index = False)


df_ip_sm_filtered.to_csv(output_file_filtered, sep = "\t")

xl_results_ip = xl_results_ip.set_index(['chr', 'start'])
df_ip_sm_filtered = df_ip_sm_filtered.rename(columns = {'xl_site': 'start'})
df_ip_sm_filtered = df_ip_sm_filtered.set_index(['chr', 'start'])

bed_signif_df = xl_results_ip.merge(df_ip_sm_filtered, on = ['chr', 'start'], how='inner').reset_index()[['chr', 'start', 'end', 'score_x']]
bed_signif_df = bed_signif_df.rename(columns = {'score_x': 'score'})

write_bw_from_bed(bed_signif_df, output_file_bw_filtered_ip, chr_length_set, chr_length_dict)

