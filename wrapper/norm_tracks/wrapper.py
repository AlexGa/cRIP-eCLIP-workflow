import numpy as np
import pysam
import pandas as pd
import csv
import re
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import itertools
from pathlib import Path
import os
import pyBigWig

class InvalidSetException(Exception):
    "Raised when two sets are not equal"
    pass

def read_genome(fasta_file):
    """
    Reading genome fasta file and return a dict with 
    chromosome names and length of the chromosomes
    """
    chr_dict = {}
    chr_length = 0
    chr_name = ""
    with open(fasta_file, 'r') as read_fa:
        for line in read_fa.readlines():
            line = line.strip()
            if line.startswith('>'):
                if chr_length > 0:
                    chr_dict[chr_name] = chr_length
                    chr_length = 0
                chr_name = line.split(' ')[0].replace('>', '')
            else:
                chr_length = chr_length + len(line)
    chr_dict[chr_name] = chr_length
    return(chr_dict)

def paired_read_number(bam_file):
    read_stats = pysam.flagstat(bam_file)
    [line.split("+") for line in read_stats.split("\n")]# 12.Element
    proper_pairs = read_stats.split("\n")
    number_prob_paired = int([line.split("+") for line in read_stats.split("\n")][11][0])
    return(number_prob_paired * 2)

# def count_bin_reads(bamFile, chr, start, end):
#     count_vec = np.zeros(end-start)
#     bam_content = pysam.Samfile(bamFile, "rb")
#     pile_iter = bam_content.pileup(contig = chr, start = start, stop = end, stepper = "nofilter")
#     i = 0
#     for pileupcolumn in pile_iter:
#         count_vec[pileupcolumn.pos] = pileupcolumn.n
#         '''
#         print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
#         '''
#     return(count_vec)

def calculate_norm_tracks(bamFile_fwd, bamFile_rev, bwFile_fwd, bwFile_rev, chr_dict, scale = 1):#, chr, start, end):
    """
    Calculate sequencing depth
    """
    libsize_fwd = paired_read_number(bamFile_fwd)
    libsize_rev = paired_read_number(bamFile_rev)
    libsize = (libsize_fwd + libsize_rev)/scale
    """
    Count raw reads and normalize reads over lib size
    for each chromosome of the genome
    """
    tracks_dict = {}
    big_raw_fwd = pyBigWig.open(bwFile_fwd)
    big_raw_rev = pyBigWig.open(bwFile_rev)
    for chr_name in chr_dict.keys():
        """ 
        count base coverage for each chromosome in the genome
        """
        raw_fwd = big_raw_fwd.values(chr_name, 0, chr_dict[chr_name]) #count_bin_reads(bamFile_fwd, chr_name, 0, chr_dict[chr_name])
        raw_rev = big_raw_rev.values(chr_name, 0, chr_dict[chr_name]) #count_bin_reads(bamFile_rev, chr_name, 0, chr_dict[chr_name])
        raw_df = pd.DataFrame({'fwd': raw_fwd,
                               'rev': raw_rev})
        """
        Normalize bin counts
        """
        norm_fwd = pd.DataFrame(raw_fwd) / libsize
        norm_rev = pd.DataFrame(raw_rev) / libsize
        norm_df = pd.DataFrame({'fwd': norm_fwd[0].tolist(),
                                'rev': norm_rev[0].tolist()})
        """
        Normalize bin counts
        """
        rpm_fwd = norm_df['fwd'] * 1e6
        rpm_rev = norm_df['rev'] * 1e6
        rpm_df = pd.DataFrame({'fwd': rpm_fwd.tolist(),
                               'rev': rpm_rev.tolist()})
        """
        sum into data frame
        """
        tracks_dict[chr_name] = {'raw': raw_df, 'norm': norm_df, 'rpm': rpm_df}
    return(tracks_dict)

def log2fc_tracks(ip_dict, sm_dict):
    lfc2_dict = {}
    for chr_name in ip_dict.keys():
        lfc2_dict[chr_name] = np.log2(ip_dict[chr_name]['norm'] + 1e-6) - np.log2(sm_dict[chr_name]['norm'] + 1e-6)
    return(lfc2_dict)

def subtracted_tracks(ip_dict, sm_dict):
    subtracted_dict = {}
    for chr_name in ip_dict.keys():
        subtracted_dict[chr_name] = ip_dict[chr_name]['norm'] - sm_dict[chr_name]['norm']
    return(subtracted_dict)

def kl_norm_tracks(ip_dict, lfc_dict):
    kl_dict = {}
    for chr_name in ip_dict.keys():
        kl_dict[chr_name] = ip_dict[chr_name]['norm'] * lfc_dict[chr_name]
    return(kl_dict)

def write_bed(filename, gen_dict, values):
    """
    Write bed file 
    """
    chroms = list(itertools.chain(*[[key]*value for key, value in gen_dict.items()]))
    starts = list(itertools.chain(*[range(1, value+1) for value in gen_dict.values()]))
    bed_df = pd.DataFrame({'chr': chroms,
                           'start': starts,
                           'end': starts,
                           'strand': '*',
                           'score': values})
    bed_df.to_csv(filename, sep = "\t", header = False, index = False)

def write_bw(filename, gen_dict, values):
    """
    Write bigwig file
    """
    bw = pyBigWig.open(filename, "w")
    bw.addHeader(list(gen_dict.items()))
    chroms = list(itertools.chain(*[[key]*value for key, value in gen_dict.items()]))
    starts = list(itertools.chain(*[range(1, value+1) for value in gen_dict.values()]))
    bw.addEntries(chroms, starts = [value - 1 for value in starts], ends = starts, values = values)
    bw.close()

# ip_fwd_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_N_rep2/IP-cRIP_N_rep2_fwd.bam"
# ip_rev_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_N_rep2/IP-cRIP_N_rep2_rev.bam"

# sm_fwd_bam = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_N_rep2/SM-cRIP_N_rep2_fwd.bam"
# sm_rev_bam = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_N_rep2/SM-cRIP_N_rep2_rev.bam"

# ip_fwd_bw = "results/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep2-IP-SM/IP-cRIP_N_rep2_fwd_signi.bw"
# ip_rev_bw = "results/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep2-IP-SM/IP-cRIP_N_rep2_rev_signif.bw"

# sm_fwd_bw = "results/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep2-IP-SM/SM-cRIP_N_rep2_fwd_raw.bw"
# sm_rev_bw = "results/XL_sites_from_peaks/only_sars_cov2_woFilter/cRIP_N_rep2-IP-SM/SM-cRIP_N_rep2_rev_raw.bw"
# genome_fa = "references/fasta/Sars_cov_2.ASM985889v3.dna.toplevel.fa"

# output_prefix_raw = "results/normed_tracks_XL_sites_from_peaks_signif/only_sars_cov2_woFilter/Raw/cRIP_N_rep2/cRIP_N_rep2_RAW"
# output_prefix_rpm = "results/normed_tracks_XL_sites_from_peaks_signif/only_sars_cov2_woFilter/Raw/cRIP_N_rep2/cRIP_N_rep2_RPM"
# output_prefix_lfc = "results/normed_tracks_XL_sites_from_peaks_signif/only_sars_cov2_woFilter/LFC/cRIP_N_rep2/cRIP_N_rep2_LFC"
# output_prefix_kl = "results/normed_tracks_XL_sites_from_peaks_signif/only_sars_cov2_woFilter/KL/cRIP_N_rep2/cRIP_N_rep2_KL"
# output_prefix_subtract = "results/normed_tracks_XL_sites_from_peaks_signif/only_sars_cov2_woFilter/Raw/cRIP_N_rep2_SUBTRACT/cRIP_N_rep2_SUBTRACT"

# ip_fwd_bw = "results/bigwig/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bw"
# ip_rev_bw = "results/bigwig/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_rev.bw"

# sm_fwd_bw = "results/bigwig/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bw"
# sm_rev_bw = "results/bigwig/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_rev.bw"

# ip_fwd_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bam"
# ip_rev_bam = "results/dedup/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2/IP-cRIP_fifth_rep1_NSP9_12h_CTRL2_rev.bam"

# sm_fwd_bam = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_fwd.bam"
# sm_rev_bam = "results/dedup/only_sars_cov2_woFilter/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2/SM-cRIP_fifth_rep1_NSP9_12h_CTRL2_rev.bam"

# genome_fa = "references/fasta/Sars_cov_2.ASM985889v3.dna.toplevel.fa"

ip_fwd_bam = snakemake.input.ip_bam[0]
ip_rev_bam = snakemake.input.ip_bam[1]

sm_fwd_bam = snakemake.input.sm_bam[0]
sm_rev_bam = snakemake.input.sm_bam[1]

ip_fwd_bw = snakemake.input.ip_bw[0]
ip_rev_bw = snakemake.input.ip_bw[1]

sm_fwd_bw = snakemake.input.sm_bw[0]
sm_rev_bw = snakemake.input.sm_bw[1]


bamFile = pysam.AlignmentFile(ip_fwd_bam, "rb")
gen_dict = dict(zip(*[list(bamFile.references), list(bamFile.lengths)]))


output_prefix_raw = snakemake.params.raw
output_prefix_rpm = snakemake.params.rpm
output_prefix_lfc = snakemake.params.lfc
output_prefix_kl = snakemake.params.kl
output_prefix_subtract = snakemake.params.subtract

scaleLibSize = 1
if snakemake.params.get("onlyXL", ""):
    scaleLibSize = 2

# output_prefix_raw = '/vol/projects/agabel/NSP9_new/results/TrackNormalizations/Raw/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_8h_CTRL2'
# output_prefix_lfc = '/vol/projects/agabel/NSP9_new/results/TrackNormalizations/LFC/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_8h_CTRL2'
# output_prefix_kl = '/vol/projects/agabel/NSP9_new/results/TrackNormalizations/KL/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_8h_CTRL2'
# output_prefix_subtract = '/vol/projects/agabel/NSP9_new/results/TrackNormalizations/Subtract/only_sars_cov2_woFilter/IP-cRIP_fifth_rep1_NSP9_8h_CTRL2'

ip_bin_dict = calculate_norm_tracks(ip_fwd_bam, ip_rev_bam, ip_fwd_bw, ip_rev_bw, gen_dict, scale = scaleLibSize) # chr = "MN908947.3", start = 0, end = 29903)
sm_bin_dict = calculate_norm_tracks(sm_fwd_bam, sm_rev_bam, sm_fwd_bw, sm_rev_bw, gen_dict, scale = scaleLibSize) # chr = "MN908947.3", start = 0, end = 29903)

try:
    if set(ip_bin_dict.keys()) != set(sm_bin_dict.keys()):
        raise InvalidSetException
except InvalidSetException:
    print("Exception occurred: Could not find the same chromosome names in IP and SMI")
    print("IP chromosomes: ", ", ".join(list(ip_bin_dict.keys())))
    print("SMI chromosomes: ", ", ".join(list(sm_bin_dict.keys())))

lfc_dict = log2fc_tracks(ip_bin_dict, sm_bin_dict)
subtracted_dict = subtracted_tracks(ip_bin_dict, sm_bin_dict)
kl_dict = kl_norm_tracks(ip_bin_dict, lfc_dict)

# l2fc_norm = np.log2((ip_bin_array[1] / sm_bin_array[1]).fillna(1))
# subtract_norm = ip_bin_array[1] - sm_bin_array[1]
# kl_norm = ip_bin_array[1] * l2fc_norm

if not os.path.exists(os.path.dirname(output_prefix_raw)):
    Path(os.path.dirname(output_prefix_raw)).mkdir(parents=True, exist_ok=True)

if not os.path.exists(os.path.dirname(output_prefix_lfc)):
    Path(os.path.dirname(output_prefix_lfc)).mkdir(parents=True, exist_ok=True)

if not os.path.exists(os.path.dirname(output_prefix_kl)):
    Path(os.path.dirname(output_prefix_kl)).mkdir(parents=True, exist_ok=True)

if not os.path.exists(os.path.dirname(output_prefix_subtract)):
    Path(os.path.dirname(output_prefix_subtract)).mkdir(parents=True, exist_ok=True)

for strand in ('fwd', 'rev'):
    lfc_list = list(itertools.chain(*[sub_dict[strand].tolist() for sub_dict in lfc_dict.values()]))
    lfc_file = output_prefix_lfc + "_" + strand + ".bw"
    write_bw(lfc_file, gen_dict, lfc_list)
    subtracted_list = np.array(list(itertools.chain(*[sub_dict[strand].tolist() for sub_dict in subtracted_dict.values()])))
    subtracted_file = output_prefix_subtract + "_" + strand + ".bw"
    write_bw(subtracted_file, gen_dict, subtracted_list)
    kl_list = np.array(list(itertools.chain(*[sub_dict[strand].tolist() for sub_dict in kl_dict.values()])))
    kl_file = output_prefix_kl + "_" + strand + ".bw"
    write_bw(kl_file, gen_dict, kl_list)
    ip_list = np.array(list(itertools.chain(*[sub_dict['raw'][strand].tolist() for sub_dict in ip_bin_dict.values()])))
    ip_file = output_prefix_raw + "_IP_" + strand + ".bw"
    write_bw(ip_file, gen_dict, ip_list)
    sm_list = np.array(list(itertools.chain(*[sub_dict['raw'][strand].tolist() for sub_dict in sm_bin_dict.values()])))
    sm_file = output_prefix_raw + "_SM_" + strand + ".bw"
    write_bw(sm_file, gen_dict, sm_list)
    # ip_list = np.array(list(itertools.chain(*[sub_dict['norm'][strand].tolist() for sub_dict in ip_bin_dict.values()])))
    # ip_file = output_prefix_raw + "_LIBSIZE_NORM_IP_" + strand + ".bw"
    # write_bw(ip_file, gen_dict, ip_list)
    # sm_list = np.array(list(itertools.chain(*[sub_dict['norm'][strand].tolist() for sub_dict in sm_bin_dict.values()])))
    # sm_file = output_prefix_raw + "_LIBSIZE_NORM_SM_" + strand + ".bw"
    # write_bw(sm_file, gen_dict, sm_list)
    ip_list = np.array(list(itertools.chain(*[sub_dict['rpm'][strand].tolist() for sub_dict in ip_bin_dict.values()])))
    ip_file = output_prefix_rpm + "_IP_" + strand + ".bw"
    write_bw(ip_file, gen_dict, ip_list)
    sm_list = np.array(list(itertools.chain(*[sub_dict['rpm'][strand].tolist() for sub_dict in sm_bin_dict.values()])))
    sm_file = output_prefix_rpm + "_SM_" + strand + ".bw"
    write_bw(sm_file, gen_dict, sm_list)


