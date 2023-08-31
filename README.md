[![DOI](https://zenodo.org/badge/682779021.svg)](https://zenodo.org/badge/latestdoi/682779021)

# Analysis of cRIP and eCLIP data for the publication

Nora Schmidt, Sabina Ganskih, Yuanjie Wei, Alexander Gabel, Sebastian Zielinski, Hasmik Keshishian, Caleb A. Lareau, Liv Zimmermann, Jana Makroczyova, Cadence Pearce, Karsten Krey, Thomas Hennig, Sebastian Stegmaier, Lambert Moyon, Marc Horlacher, Simone Werner, Jens Aydin, Marco Olguin-Nava, Ramya Potabattula, Anuja Kibe, Lars Dölken, Redmond P. Smyth, Neva Caliskan, Annalisa Marsico, Christine Krempl, Jochen Bodem, Andreas Pichlmair, Steven A. Carr, Petr Chlanda, Florian Erhard, and Mathias Munschauer. 

**SND1 binds SARS-CoV-2 negative-sense RNA and promotes viral RNA synthesis through NSP9**.

# cRIP and eCLIP workflow

## General description

This pipeline was used for analysing eCLIP and cRIP data to identify RNA binding regions in the Sars-Cov-2 genome. 
It is based on the read preprocessing and mapping procedures presented in the [eCLIP pipeline](https://www.encodeproject.org/pipelines/ENCPL357ADL/) from [Gene Yeo's lab](https://github.com/YeoLab/eCLIP) while the downstream analysis such as peak calling, calculating differential binding affinities, analysing covalent linkage sites, and visualizing read coverages were specifically developed and implemented for the identification of RBP regions in viral genomes.

By adjusting the parameters of this workflow, it can be used to identify RBP regions in other genomes of interest (host and/or virus). Besides the basic analysis steps (e.g. clipping, trimming, and mapping of fastq reads), this workflow calculates:

- the genomic location of RBP regions (peaks)
- corresponding covalent linkage sites (cl-sites)
- determine differential binding affinities of peaks and cl-sites between different treatments
- normalization of read coverages to vizualize the enrichment of IP over SMI


## Short description

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

It is recommended to install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and using conda it is an easy task to install [snakemake](https://snakemake.github.io/).

All tools needed for the analysis are installed during the run of Snakemake while all specifications for the different analysis are written in separate yaml environment files.

## Directory structure

Please define the following directory structure so that the snakemake workflow is able to detect all necessary files.

```bash
├── config
│   ├── cluster.json
│   ├── config.yaml
│   ├── peak_comparisons.tsv
│   ├── sample_comparisons.yaml
│   └── units.tsv
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

# Citation:

Schmidt, Nora;  Ganskih, Sabina; Wei, Yuanjie; Gabel, Alexander; et al. (2023): SND1 binds SARS-CoV-2 negative-sense RNA and promotes viral RNA synthesis through NSP9
