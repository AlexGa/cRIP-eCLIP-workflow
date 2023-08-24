# cRIP and eCLIP workflow

This pipeline was used for analysing eCLIP and cRIP data to identify RNA binding regions in the Sars-Cov-2 genome. By adjusting the parameters of this workflow it can be used to identify RBP regions in other genomes of interest. Besides the basic analysis steps such as clipping, trimming, and mapping of fastq reads, this workflow predicts the genomic location of RBP regions (peaks), their corresponding covalent linkage sites (cl-sites), and is able to determine differential binding affinities on peak and cl-site level between for samples from different treatments. 
As a result of this workflow different normalization methods for peak coverages are produced to vizualize the enrichment of IP oder size matched input.

The cRIP and eCLIP analyses of the following publication are based on this snakemake workflow.


#### SND1 binds SARS-CoV-2 negative-sense RNA and promotes viral RNA synthesis through NSP9

Nora Schmidt<sup>1\*</sup>, Sabina Ganskih<sup>1\*</sup>, Yuanjie Wei<sup>1\*</sup>, Alexander Gabel<sup>1\*</sup>, Sebastian Zielinski<sup>1</sup>, Hasmik Keshishian<sup>2</sup>, Caleb A. Lareau<sup>3</sup>, Liv Zimmermann<sup>4</sup>, Jana Makroczyova<sup>4</sup>, Cadence Pearce<sup>2</sup>, Karsten Krey<sup>5</sup>, Thomas Hennig<sup>6</sup>, Sebastian Stegmaier<sup>1</sup>, Lambert Moyon<sup>7</sup>, Marc Horlacher<sup>7</sup>,Simone Werner<sup>1</sup>, Jens Aydin<sup>1</sup>, Marco Olguin-Nava<sup>1</sup>, Ramya Potabattula<sup>8</sup>, Anuja Kibe<sup>1</sup>, Lars Dölken<sup>6</sup>, Redmond P. Smyth<sup>1</sup>, Neva Caliskan<sup>1</sup>, Annalisa Marsico<sup>7</sup>, Christine Krempl<sup>6</sup>, Jochen Bodem<sup>6</sup>, Andreas Pichlmair<sup>5,9</sup>, Steven A. Carr<sup>2</sup>, Petr Chlanda<sup>4</sup>, Florian Erhard<sup>6,10</sup>, and Mathias Munschauer<sup>1,11§</sup>

**Affiliations:**

<sup>1</sup> Helmholtz Institute for RNA-based Infection Research (HIRI), Helmholtz-Center for Infection Research (HZI), Würzburg, Germany </ br> 
<sup>2</sup> Broad Institute of MIT and Harvard, Cambridge, MA 02142, USA </ br>
<sup>3</sup> School of Medicine, Stanford University, Palo Alto, CA 94305, USA </ br>
<sup>4</sup> Schaller Research Groups, Department of Infectious Diseases, Virology, Heidelberg University Hospital, Heidelberg, Germany </ br>
<sup>5</sup> School of Medicine, Institute of Virology, Technical University of Munich, Munich, Germany </ br>
<sup>6</sup> Institute for Virology and Immunobiology, Julius-Maximilians-University Würzburg, Würzburg, Germany </ br>
<sup>7</sup> Computational Health Center, Helmholtz Center Munich, Germany </ br>
<sup>8</sup> Institute of Human Genetics, Julius Maximilians University, Würzburg, Germany </ br>
<sup>9</sup> German Center for Infection Research (DZIF), Munich Partner Site, Munich, Germany </ br>
<sup>10</sup> Faculty for Computer and Data Science, University of Regensburg, Bajuwarenstr. 4, 93053 Regensburg, Germany </ br>
<sup>11</sup> Faculty of Medicine, University of Würzburg, Würzburg, Germany</ br>
<sup>\*</sup> These authors contributed equally</ br>
<sup>§</sup> Corresponding author

# Short description

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

# Recommended requirements

It is recommended to install [conda] (https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and using conda it is an easy task to install [snakemake] (https://snakemake.github.io/).

All tools needed for the analysis are installed during the run of Snakemake while all specifications for the different analysis are written in separate yaml environment files.

# Directory structure

Please define the following directory structure so that the snakemake workflow is able to detect all necessary files.

```bash
├── config
│		 ├── cluster.json
│		 ├── config.yaml
│		 ├── peak_comparisons.tsv
│		 ├── sample_comparisons.yaml
│		 └── units.tsv
├── data
├── crip-eclip-workflow
├── Log_Err
├── Log_Out
├── references
│		 ├── fasta
│		 │		 ├── GRCh38_SarsCov2.dna.toplevel.fa
│		 │		 ├── Homo_sapiens.GRCh38.dna_sm.toplevel.fa
│		 │		 ├── Sars_cov_2.ASM985889v3.dna.toplevel.fa
│		 ├── gtf
│		 │		 ├── GRCh38_SarsCov2.gtf
│		 │		 ├── Homo_sapiens.GRCh38.106.gtf
│		 │		 └── Sars_cov_2.ASM985889v3.101.gtf
│		 └── STAR_INDEX
│		     ├── homo_sapiens_sars_cov2
└── Snakefile
```
