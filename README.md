# cRIP and eCLIP workflow

## General description

This pipeline was used for analysing eCLIP and cRIP data to identify RNA binding regions in the Sars-Cov-2 genome. 

By adjusting the parameters of this workflow it can be used to identify RBP regions in other genomes of interest. Beside basic analyses (clipping, trimming, and mapping of fastq reads), this workflow calculates:

- the genomic location of RBP regions (peaks)
- corresponding covalent linkage sites (cl-sites)
- determine differential binding affinities of peaks and cl-sites between different treatments
- normalization of read coverages to vizualize the enrichment of IP over SMI


## Current publication using this workflow

### SND1 binds SARS-CoV-2 negative-sense RNA and promotes viral RNA synthesis through NSP9

Nora Schmidt, Sabina Ganskih, Yuanjie Wei, Alexander Gabel, Sebastian Zielinski, Hasmik Keshishian, Caleb A. Lareau, Liv Zimmermann, Jana Makroczyova, Cadence Pearce, Karsten Krey, Thomas Hennig, Sebastian Stegmaier, Lambert Moyon, Marc Horlacher,Simone Werner, Jens Aydin, Marco Olguin-Nava, Ramya Potabattula, Anuja Kibe, Lars Dölken, Redmond P. Smyth, Neva Caliskan, Annalisa Marsico, Christine Krempl, Jochen Bodem, Andreas Pichlmair, Steven A. Carr, Petr Chlanda, Florian Erhard, and Mathias Munschauer

## Short workflow description

- UMI clipping and storing performed by DeBarCode.jar (see java directory)
- adapter clipping and quality trimming with cutadapt
- map reads to host and viral genome with star
- removal of PCR duplicates by picard
- separation of mapped reads from host and virus
- Macs2 peak calling and confirmation of significant peaks
- differential binding of peaks between different samples
- extraction of cl-sites
- calculation of enriched cl-sites
- normalization of mapped reads IP vs. SMI for better visualisation  

## Recommended requirements

It is recommended to install [conda] (https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and using conda it is an easy task to install [snakemake] (https://snakemake.github.io/).

All tools needed for the analysis are installed during the run of Snakemake while all specifications for the different analysis are written in separate yaml environment files.

## Directory structure

Please define the following directory structure so that the snakemake workflow is able to detect all necessary files.

```bash
├── config
││   ├── cluster.json
││   ├── config.yaml
││   ├── peak_comparisons.tsv
││   ├── sample_comparisons.yaml
││   └── units.tsv
├── data
├── crip-eclip-workflow
├── Log_Err
├── Log_Out
├── references
│    ├── fasta
│    │   ├── GRCh38_SarsCov2.dna.toplevel.fa
│    │   ├── Homo_sapiens.GRCh38.dna_sm.toplevel.fa
│    │   ├── Sars_cov_2.ASM985889v3.dna.toplevel.fa
│    ├── gtf
│    │   ├── GRCh38_SarsCov2.gtf
│    │   ├── Homo_sapiens.GRCh38.106.gtf
│    │   └── Sars_cov_2.ASM985889v3.101.gtf
│    └── STAR_INDEX
│        ├── homo_sapiens_sars_cov2
└── Snakefile
```
