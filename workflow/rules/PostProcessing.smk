
# Import Python Modules
import os
import pandas as pd

SAMPLES = pd.read_csv(config["INPUT"], index_col='SampleID')
OUTDIR = config["outdir"]

## Post Processing Options
chunk_size = config['params']['chunk_size']
overwrite = config['params']['overwrite']
max_mito_cn = config['params']['max_mito_cn']
minimum_homref_coverage = config['params']['minimum_homref_coverage']
vaf_filter_threshold = config['params']['vaf_filter_threshold']
min_het_threshold = config['params']['min_het_threshold']
min_hom_threshold = config['params']['min_hom_threshold']
min_mito_cn = config['params']['min_mito_cn']

# Gnomad Resources
ARTIFACT_PRONE_SITES = config["resources"]["ARTIFACT_PRONE_SITES"]
variant_context = config["resources"]["variant_context"]
phylotree = config["resources"]["phylotree"]
pon_mt_trna = config["resources"]["pon_mt_trna"]
mitotip = config["resources"]["mitotip"]
mt_dbsnp154 = config["resources"]["mt_dbsnp154"]

rule all:
    input:
        expand("{outdir}/joint/final/sample_vcf.vcf.bgz", outdir = OUTDIR)

rule hail_inputs:
    input:
        coverage = expand("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        vcf = expand("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        mtdnacn = expand("{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        haplocheck = expand("{outdir}/samples/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.raw.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist())
    output:
        coverage_tsv = "{outdir}/joint/sample_coverage_paths.txt",
        vcf_tsv = "{outdir}/joint/sample_vcf_paths.txt",
        meta = "{outdir}/joint/sample_metadata.txt",
    params:
        samples = expand("{SAMPLE_IDS}", SAMPLE_IDS = SAMPLES.index.tolist())
    conda: "../envs/r.yaml"
    log: "{outdir}/joint/logs/hail_inputs.log"
    script: "../scripts/hail_inputs.R"

rule annotate_coverage:
    input:
        input_tsv = rules.hail_inputs.output.coverage_tsv,
    output:
        ht = directory("{outdir}/joint/annotate_coverage/all_chrM_annotate_coverage.ht"),
        mt = directory("{outdir}/joint/annotate_coverage/all_chrM_annotate_coverage.mt"),
        all = "{outdir}/joint/annotate_coverage/all_chrM_annotate_coverage.tsv",
        sample = "{outdir}/joint/annotate_coverage/all_chrM_annotate_coverage_sample_level.txt"
    params:
        temp_dir = "{outdir}/joint/temp",
        chunk_size = chunk_size,
        overwrite = overwrite
    conda: '../envs/hail.yaml'
    log:
        python_logger = "{outdir}/joint/logs/python_annotate_coverage.log",
        hail_logs = "{outdir}/joint/logs/hail_annotate_coverage.log",
    script: '../scripts/annotate_coverage.py'

rule combine_vcfs:
    input:
        vcf_paths = rules.hail_inputs.output.vcf_tsv,
        coverage = rules.annotate_coverage.output.mt,
        ARTIFACT_PRONE_SITES = ARTIFACT_PRONE_SITES,
    output:
        vcf = "{outdir}/joint/combine_vcfs/all_chrM.vcf.bgz",
        mt = directory("{outdir}/joint/combine_vcfs/all_chrM.mt"),
        temp = temp(directory("{outdir}/joint/combine_vcfs/temp"))
    params:
        VCF_COL_NAME = "VCF",
        temp_dir = "{outdir}/joint/combine_vcfs/temp",
        output_bucket = "{outdir}/joint/combine_vcfs",
        chunk_size = chunk_size,
        minimum_homref_coverage = minimum_homref_coverage,
        overwrite = overwrite,
        participants_to_subset = None
    conda: '../envs/hail.yaml'
    log:
        python_logger = "{outdir}/joint/logs/python_combine_vcfs.log",
        hail_logs = "{outdir}/joint/logs/hail_combine_vcfs.log",
    script: '../scripts/combine_vcfs.py'

rule add_annotations:
    input:
        vcf = rules.combine_vcfs.output.vcf,
        mt = rules.combine_vcfs.output.mt,
        meta = rules.hail_inputs.output.meta,
        variant_context = variant_context,
        phylotree = phylotree,
        pon_mt_trna = pon_mt_trna,
        mitotip = mitotip,
        mt_dbsnp154 = mt_dbsnp154
    output:
        mt = directory("{outdir}/joint/final/annotated_combined.mt"),
        ht = directory("{outdir}/joint/final/combined_sites_only.ht"),
        txt = "{outdir}/joint/final/combined_sites_only.txt",
        sites_vcf = "{outdir}/joint/final/combined_sites_only.vcf.bgz",
        sample = "{outdir}/joint/final/sample_annotations.txt",
        sample_vcf = "{outdir}/joint/final/sample_vcf.vcf.bgz",
    params:
        OUTPUT_DIR = "{outdir}/joint/final",
        overwrite = overwrite,
        min_het_threshold = min_het_threshold,
        min_hom_threshold = min_hom_threshold,
        vaf_filter_threshold = vaf_filter_threshold,
        min_mito_cn = min_mito_cn,
        max_mito_cn = max_mito_cn
    conda: '../envs/hail.yaml'
    log:
        python_logger = "{outdir}/joint/logs/python_add_annotations.log",
        hail_logs = "{outdir}/joint/logs/hail_add_annotations.log",
    script: "../scripts/add_annotations.py"
