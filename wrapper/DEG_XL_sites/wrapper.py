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

# bamFile_ctrl = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_CTRL2/IP-cRIP_fourth_NSP9_8h_CTRL2_rev.bam" 
# bamFile_ko = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV_rev.bam"
# bed_ctrl = "results/XL_sites_bedgraph_from_pub_peaks/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_CTRL2/IP-cRIP_fourth_NSP9_8h_CTRL2_rev.bedgraph" 
# bed_ko = "results/XL_sites_bedgraph_from_pub_peaks/only_sars_cov2_woFilter/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV/IP-cRIP_fourth_NSP9_8h_KO6_plus_eV_rev.bedgraph"
    
# output_file_raw = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_rev_raw.tsv"
# output_file_filtered = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_rev_signif.tsv"
# output_file_narrow = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/cRIP_fourth_NSP9_8h_rev_signif.bed"
# output_file_bw_ctrl = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/CTRL-cRIP_fourth_NSP9_8h_rev_signif.bw"
# output_file_bw_ko = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fourth_NSP9_8h/KO-cRIP_fourth_NSP9_8h_rev_signif.bw"


# bamFile_ctrl = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bam"
# bamFile_ko = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_KO6_plus_eV/IP-cRIP_fifth_rep1_NSP9_12h_KO6_plus_eV_fwd.bam"
# bed_ctrl = "results/XL_sites_bedgraph_from_pub_peaks/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bedgraph"
# bed_ko = "results/XL_sites_bedgraph_from_pub_peaks/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bedgraph"

# output_file_raw = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_raw.tsv"
# output_file_filtered = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_signif.tsv"
# output_file_narrow = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/cRIP_fifth_rep1_NSP9_12h_fwd_signif.bed"
# output_file_bw_ctrl = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/CTRL-cRIP_fifth_rep1_NSP9_12h_fwd_signif.bw"
# output_file_bw_ko = "results/DEG_XL_sites_from_pub_peaks/only_sars_cov2_woFilter/cRIP_fifth_rep1_NSP9_12h/KO-cRIP_fifth_rep1_NSP9_12h_fwd_signif.bw"

# orientation = "rev"
# count_threshold = 10

bamFile_ctrl = snakemake.input.ctrl
bamFile_ko = snakemake.input.ko

bed_ctrl = snakemake.input.ctrl_bed
bed_ko = snakemake.input.ko_bed

xl_results_ctrl = pd.read_csv(bed_ctrl, delimiter = "\t", header = None)
xl_results_ctrl.columns = ['chr', 'start', 'end', 'score']
xl_results_ctrl = xl_results_ctrl.set_index(['chr', 'start']) #.index = xl_results_ctrl['start']


xl_results_ko = pd.read_csv(bed_ko, delimiter = "\t", header = None)
xl_results_ko.columns = ['chr', 'start', 'end', 'score']
xl_results_ko = xl_results_ko.set_index(['chr', 'start']) #= xl_results_ko['start']


# check if bed file with significant peaks is given
# if so than only consider reads/XL-sites overlapping with significant 
# peaks
# bedFile = snakemake.input.get("bed", "")


orientation = snakemake.params.orientation
count_threshold = snakemake.params.count_threshold

output_file_raw = snakemake.output.tsv_raw
output_file_filtered = snakemake.output.tsv_filtered
output_file_narrow = snakemake.output.bed

output_file_bw_ctrl = snakemake.output.bw_filtered_ctrl
output_file_bw_ko = snakemake.output.bw_filtered_ko


xl_reads_ctrl = total_xl_reads(bamFile_ctrl)
xl_reads_ko = total_xl_reads(bamFile_ko)


bamFile = pysam.AlignmentFile(bamFile_ctrl, "rb")
chr_length_set = list(zip(*[list(bamFile.references), list(bamFile.lengths)]))
chr_length_dict = dict(zip(*[list(bamFile.references), list(bamFile.lengths)]))


ctrl_cov = xl_results_ctrl['score']
ko_cov = xl_results_ko['score']

# idx_suff_coverage = (ctrl_cov > count_threshold)

# If no reads left after filtering, exit the script
# if len(ctrl_cov[idx_suff_coverage]) == 0:
#     with open(output_file_raw, 'w') as out:
#         print("Unsufficient read support!")
#         print("Reads in CTRL: " + str(xl_reads_ctrl))
#         print("Reads in KO: " + str(xl_reads_ko))
#     with open(output_file_filtered, 'w') as out:
#         print("Unsufficient read support!")
#         print("Reads in CTRL: " + str(xl_reads_ctrl))
#         print("Reads in KO: " + str(xl_reads_ko))
#     with open(output_file_narrow, 'w') as out:
#         print("Unsufficient read support!")
#         print("Reads in CTRL: " + str(xl_reads_ctrl))
#         print("Reads in KO: " + str(xl_reads_ko))
#     with open(output_file_bw_filtered_ctrl, 'w') as out:
#         print("Unsufficient read support!")
#         print("Reads in CTRL: " + str(xl_reads_ctrl))
#         print("Reads in KO: " + str(xl_reads_ko))
#     sys.exit(0)

# xl_results_ctrl = xl_results_ctrl[(ctrl_cov > count_threshold)]
# xl_results_ko = xl_results_ko[(ko_cov > count_threshold)]

xl_results_ctrl_ko = xl_results_ctrl.merge(xl_results_ko, on=['chr', 'start'], how='inner')

xl_results_ctrl = xl_results_ctrl_ko[['end_x', 'score_x']]
xl_results_ctrl['score_x'].fillna(0, inplace=True)
xl_results_ctrl.columns = xl_results_ko.columns

xl_results_ko = xl_results_ctrl_ko[['end_x', 'score_y']]
xl_results_ko['score_y'].fillna(0, inplace=True)
xl_results_ko.columns = xl_results_ctrl.columns

xl_results_ko = xl_results_ko.reset_index()
xl_results_ctrl = xl_results_ctrl.reset_index()

ctrl_cov = xl_results_ctrl['score']
ko_cov = xl_results_ko['score']

df_ctrl_ko = pd.DataFrame()
for chrom in chr_length_dict:
    odds_dict_chrom = pval_dict_chrom = fc_dict_chrom = np.full(np.sum(xl_results_ctrl['chr'] == chrom), np.nan) #np.zeros(np.sum(xl_results_ctrl['chr'] == chrom))
    ctrl_cov_chrom = list(xl_results_ctrl[xl_results_ctrl['chr'] == chrom]['score'])
    ko_cov_chrom = list(xl_results_ko[xl_results_ko['chr'] == chrom]['score'])
    no_ctrl_cov_chrom = list(xl_reads_ctrl - xl_results_ctrl[xl_results_ctrl['chr'] == chrom]['score'])
    no_ko_cov_chrom = list(xl_reads_ko - xl_results_ko[xl_results_ko['chr'] == chrom]['score'])
    for i in range(0, len(ctrl_cov_chrom)):   
        odds_dict_chrom[i], pval_dict_chrom[i] = fisher_exact([[ko_cov_chrom[i], no_ko_cov_chrom[i]], 
                                                               [ctrl_cov_chrom[i], no_ctrl_cov_chrom[i]]], 
                                                                  alternative = 'two-sided')
    fc_dict_chrom = np.log2((np.array(ko_cov_chrom) + 1)/(xl_reads_ko + 1)) - np.log2((np.array(ctrl_cov_chrom) + 1)/(xl_reads_ctrl + 1))
    d = {'chr': list(xl_results_ctrl[xl_results_ctrl['chr'] == chrom]['chr']),
         'xl_site': list(xl_results_ctrl[xl_results_ctrl['chr'] == chrom]['start']), 
         'count.ctrl': list(ctrl_cov_chrom), 
         'count.ko': list(ko_cov_chrom), 
         'no_xl_ctrl':  list(no_ctrl_cov_chrom), 
         'no_xl_ko':  list(no_ko_cov_chrom), 
         'log2FC(ko/ctrl)': list(fc_dict_chrom), 
         'oddsRatio': list(odds_dict_chrom), 
         'p.value': list(pval_dict_chrom)}
    df_ctrl_ko_chrom = pd.DataFrame(d)
    df_ctrl_ko = pd.concat([df_ctrl_ko, df_ctrl_ko_chrom])

df_ctrl_ko = df_ctrl_ko.reset_index()
padj_list = multipletests(list(df_ctrl_ko['p.value']), method = 'fdr_by')[1]

df_ctrl_ko = pd.concat([df_ctrl_ko, pd.DataFrame({'p.adj': padj_list})], axis = 1)

df_ctrl_ko.to_csv(output_file_raw, sep = "\t", index = False)

df_ctrl_ko_filtered = df_ctrl_ko[df_ctrl_ko['p.adj'] < 0.05]

if len(df_ctrl_ko_filtered) == 0:
    with open(output_file_filtered, 'w') as out:
        print("No significant XL sites!")
    with open(output_file_narrow, 'w') as out:
        print("No significant XL sites!")
    with open(output_file_bw_ctrl, 'w') as out:
        print("No significant XL sites!")
    with open(output_file_bw_ko, 'w') as out:
        print("No significant XL sites!")
    sys.exit(0)

df_ctrl_ko_filtered = df_ctrl_ko_filtered.assign(score = -1 * np.log10(df_ctrl_ko_filtered.loc[:,'p.adj']))
df_ctrl_ko_filtered.loc[np.isinf(df_ctrl_ko_filtered.loc[:,'score']), 'score'] = 400
df_ctrl_ko_filtered = df_ctrl_ko_filtered.assign(logpVal = -1 * np.log10(df_ctrl_ko_filtered.loc[:,'p.value']))
df_ctrl_ko_filtered.loc[np.isinf(df_ctrl_ko_filtered.loc[:,'logpVal']), 'logpVal'] = 400

narrow_dict = {'chr': list(df_ctrl_ko_filtered['chr']), 
               'start': list(df_ctrl_ko_filtered['xl_site']), 
               'end': list(df_ctrl_ko_filtered['xl_site']), 
               'name': list("xl_site_" + df_ctrl_ko_filtered['chr'] + "_" + df_ctrl_ko_filtered['xl_site'].astype(str)), 
               'score': list(df_ctrl_ko_filtered['score'] * 10), 
               'strand': ['.'] * len(df_ctrl_ko_filtered), 
               'pValue': list(df_ctrl_ko_filtered['logpVal']), 
               'qValue': list(df_ctrl_ko_filtered['score']), 
               'peak': [1] * len(df_ctrl_ko_filtered)}


narrow_peak_df = pd.DataFrame(narrow_dict)
narrow_peak_df.to_csv(output_file_narrow, sep = "\t", header = False, index = False)


df_ctrl_ko_filtered.to_csv(output_file_filtered, sep = "\t")

xl_results_ctrl = xl_results_ctrl.set_index(['chr', 'start'])
df_ctrl_ko_filtered = df_ctrl_ko_filtered.rename(columns = {'xl_site': 'start'})
df_ctrl_ko_filtered = df_ctrl_ko_filtered.set_index(['chr', 'start'])

bed_signif_df_ctrl = xl_results_ctrl.merge(df_ctrl_ko_filtered, on = ['chr', 'start'], how='inner').reset_index()[['chr', 'start', 'end', 'score_x']]
bed_signif_df_ctrl = bed_signif_df_ctrl.rename(columns = {'score_x': 'score'})

write_bw_from_bed(bed_signif_df_ctrl, output_file_bw_ctrl, chr_length_set, chr_length_dict)

bed_signif_df_ko = xl_results_ko.merge(df_ctrl_ko_filtered, on = ['chr', 'start'], how='inner').reset_index()[['chr', 'start', 'end', 'score_y']]
bed_signif_df_ko = bed_signif_df_ko.rename(columns = {'score_y': 'score'})

write_bw_from_bed(bed_signif_df_ko, output_file_bw_ko, chr_length_set, chr_length_dict)
