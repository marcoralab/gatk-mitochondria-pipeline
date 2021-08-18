# gatk-mitochondria-pipeline
Snakemake workflow for implementing [GATK best practices](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline) SNP and INDEL variant calling on the mitochondrial genome ([Laricchia et al 2021. BioRxiv](https://www.biorxiv.org/content/10.1101/2021.07.23.453510v1)).

## Usage


### Input


### Output
Final output files that are retained by the overall workflow include:
- `{SAMPLE_IDS}_mtdnacn.txt`: autosome and mitochondrial coverage and estimated mtDNAcn
- `{SAMPLE_IDS}_chrM.bam`: Subseted mitochondrial bam file from original input files
- `{SAMPLE_IDS}_chrM_per_base_coverage.tsv`: Proivdes coverage at each base - used for generating joint vcf file.
- `{SAMPLE_IDS}_chrM_contamination.raw.txt`: Output of haplochecker providing estimated contamination levels and mitochondrial haplogroups
- `{SAMPLE_IDS}_chrM_metrics.txt`: Coverage meterics of mitochondrial genome from `CollectWgsMetrics`.
- `{SAMPLE_IDS}_chrM_theoretical_sensitivity.txt`: Theoretical Sensitivity of mitochondrial genome from `CollectWgsMetrics`.
- `{SAMPLE_IDS}_chrM.duplicate_metrics`: Output from MarkDuplicates.

Intermediate files can be retained by using `--notemp`
