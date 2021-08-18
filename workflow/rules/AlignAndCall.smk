# GATK Best Practives Variant Calling pipeline
## Adapated from https://github.com/gatk-workflows/gatk4-mitochondria-pipeline
## Laricchia, K. M. et al. (2021). BioRxiv. doi:10.1101/2021.07.23.453510

# Import Python Modules
import os
from itertools import product
import pandas as pd

# Input Files
SAMPLES = pd.read_csv(config["INPUT"], index_col='SampleID')
# ID = ["SM-CJGNJ", "SM-CJEIA", "71992"]
# ID = ["TWDF-037-0006-0106925-0010007709", "SM-CJGNJ"]

# Output Files
OUTDIR = config["outdir"]

## GATK Options
max_read_length = config['gatk']["max_read_length"]
coverage_cap = config['gatk']["coverage_cap"]
max_low_het_sites = config['gatk']["max_low_het_sites"]
max_reads_per_alignment_start = config['gatk']["max_reads_per_alignment_start"]
m2_extra_args = config['gatk']["m2_extra_args"]
m2_extra_filtering_args = config['gatk']["m2_extra_filtering_args"]
max_alt_allele_count = config['gatk']["max_alt_allele_count"]
vaf_filter_threshold = config['gatk']["vaf_filter_threshold"]
f_score_beta = config['gatk']["f_score_beta"]
vaf_cutoff = config['gatk']["vaf_cutoff"]

## Reference files
### hg38 whole genome reference
hg38_ref_fasta = config["resources"]["hg38"]["ref_fasta"]
hg38_ref_dict = config["resources"]["hg38"]["ref_dict"]
hg38_ref_fasta_index = config["resources"]["hg38"]["ref_fasta_index"]
### hg19 whole genome reference
hg19_ref_fasta = config["resources"]["hg19"]["ref_fasta"]
hg19_ref_dict = config["resources"]["hg19"]["ref_dict"]
hg19_ref_fasta_index = config["resources"]["hg19"]["ref_fasta_index"]
### hg38 mt refence files
mt_dict = config["resources"]["mt"]["mt_dict"]
mt_fasta = config["resources"]["mt"]["mt_fasta"]
mt_fasta_index = config["resources"]["mt"]["mt_fasta_index"]
mt_amb = config["resources"]["mt"]["mt_amb"]
mt_ann = config["resources"]["mt"]["mt_ann"]
mt_bwt = config["resources"]["mt"]["mt_bwt"]
mt_pac = config["resources"]["mt"]["mt_pac"]
mt_sa = config["resources"]["mt"]["mt_sa"]
non_control_region_interval_list = config["resources"]["mt"]["non_control_region_interval_list"]
blacklisted_sites = config["resources"]["mt"]["blacklisted_sites"]
blacklisted_sites_index = config["resources"]["mt"]["blacklisted_sites_index"]
### hg38 mt shifted reference files
mt_shifted_dict = config["resources"]["mt_shifted"]["mt_shifted_dict"]
mt_shifted_fasta = config["resources"]["mt_shifted"]["mt_shifted_fasta"]
mt_shifted_fasta_index = config["resources"]["mt_shifted"]["mt_shifted_fasta_index"]
mt_shifted_amb = config["resources"]["mt_shifted"]["mt_shifted_amb"]
mt_shifted_ann = config["resources"]["mt_shifted"]["mt_shifted_ann"]
mt_shifted_bwt = config["resources"]["mt_shifted"]["mt_shifted_bwt"]
mt_shifted_pac = config["resources"]["mt_shifted"]["mt_shifted_pac"]
mt_shifted_sa = config["resources"]["mt_shifted"]["mt_shifted_sa"]
shift_back_chain = config["resources"]["mt_shifted"]["shift_back_chain"]
control_region_shifted_reference_interval_list = config["resources"]["mt_shifted"]["control_region_shifted_reference_interval_list"]

rule all:
    input:
        expand("{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist()),
        expand("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz", outdir = OUTDIR, SAMPLE_IDS = SAMPLES.index.tolist())

def ref_fasta_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_fasta
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_fasta
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def ref_index_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_fasta_index
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_fasta_index
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def ref_dict_input(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return hg19_ref_dict
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return hg38_ref_dict
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

def mito_contig(wildcards):
    if SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 37:
        return 'MT'
    elif SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] == 38:
        return 'chrM'
    else:
        return print(SAMPLES.loc[wildcards.SAMPLE_IDS]['hg'] + ': Build Not Recognised')

# rule CollectWgsMetrics_wgs:
#     input:
#         input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
#         input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
#         ref_fasta = ref_fasta_input,
#         ref_fasta_index = ref_index_input,
#     output:
#         meterics = '{outdir}/{SAMPLE_IDS}/{SAMPLE_IDS}_WGS_metrics.txt',
#         sens = '{outdir}/{SAMPLE_IDS}/{SAMPLE_IDS}_WGS_theoretical_sensitivity.txt'
#     params:
#         read_length = max_read_length,
#         coverage_cap = 100000,
#     conda: '../envs/gatk.yaml'
#     shell:
#         r"""
#         gatk CollectWgsMetrics \
#           -INPUT {input.input_bam} \
#           -VALIDATION_STRINGENCY SILENT \
#           -REFERENCE_SEQUENCE {input.ref_fasta} \
#           -OUTPUT {output.meterics} \
#           -USE_FAST_ALGORITHM true \
#           -READ_LENGTH {params.read_length} \
#           -COVERAGE_CAP {params.coverage_cap} \
#           -THEORETICAL_SENSITIVITY_OUTPUT WGS_theoretical_sensitivity.txt
#         """

rule mosdepth_median:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
    output: "{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.meadian.mosdepth.summary.txt"
    params:
        wd = '{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.meadian'
    conda: '../envs/mosdepth.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_mosdepth_median.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_mosdepth_median.stdout"
    shell: "(mosdepth -n --use-median -t 4 -f {input.ref_fasta} {params.wd} {input.input_bam})  2> {log.stderr} 1> {log.stdout}"

rule mosdepth_mean:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
    output: "{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.mean.mosdepth.summary.txt"
    params:
        wd = '{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}.mean'
    conda: '../envs/mosdepth.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_mosdepth_mean.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_mosdepth_mean.stdout"
    shell: "(mosdepth -n -t 4 -f {input.ref_fasta} {params.wd} {input.input_bam}) 2> {log.stderr} 1> {log.stdout}"

rule mtdnacn:
    input:
        median = rules.mosdepth_median.output,
        mean = rules.mosdepth_mean.output
    output:
        mtdnacn = "{outdir}/samples/{SAMPLE_IDS}/mosdepth/{SAMPLE_IDS}_mtdnacn.txt"
    params:
        id = "{SAMPLE_IDS}"
    conda: '../envs/r.yaml'
    log: "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_mtdnacn.log"
    script: "../scripts/mtdnacn.R"

# "Subsets a whole genome bam to just Mitochondria reads"
rule SubsetBamToChrM:
    input:
        input_bam = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['wgs'],
        input_bam_index = lambda wildcards: SAMPLES.loc[wildcards.SAMPLE_IDS]['index'],
        ref_fasta = ref_fasta_input,
        ref_fasta_index = ref_index_input,
        ref_dict = ref_dict_input,
    params:
        contig_name = mito_contig,
    output:
        bam = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.bam',
        bai = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.bai',
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SubsetBamToChrM.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SubsetBamToChrM.stdout"
    shell:
        r"""
        gatk PrintReads \
          -R {input.ref_fasta} \
          -L {params.contig_name} \
          --read-filter MateOnSameContigOrNoMappedMateReadFilter \
          --read-filter MateUnmappedAndUnmappedReadFilter \
          -I {input.input_bam} \
          --read-index {input.input_bam_index} \
          -O {output.bam} \
          2> {log.stderr} 1> {log.stdout}
        """

rule RevertSam:
    input:
        input_bam = rules.SubsetBamToChrM.output.bam,
        input_bai = rules.SubsetBamToChrM.output.bai
    output:
        unmapped_bam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_unmapped.bam'),
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_RevertSam.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_RevertSam.stdout"
    shell:
        r"""
        gatk RevertSam \
            -INPUT {input.input_bam} \
            -OUTPUT_BY_READGROUP false \
            -OUTPUT {output.unmapped_bam} \
            -VALIDATION_STRINGENCY LENIENT \
            -ATTRIBUTE_TO_CLEAR FT \
            -ATTRIBUTE_TO_CLEAR CO \
            -SORT_ORDER queryname \
            -RESTORE_ORIGINAL_QUALITIES false \
            2> {log.stderr} 1> {log.stdout}
        """

module AlignmentPipeline:
    snakefile: "AlignmentPipeline.smk"

use rule AlignAndMarkDuplicates from AlignmentPipeline as AlignToMt with:
    input:
        unmapped_bam = rules.RevertSam.output.unmapped_bam,
        mt_ref_fasta = mt_fasta
    output:
        fastq = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_unmapped.fastq'),
        aln_sam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_aln.sam'),
        mba_bam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_mba.bam'),
        md_bam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_md.bam'),
        metrics_filename = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.duplicate_metrics',
        sorted_bam = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_sorted.bam',
        sorted_bai = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_sorted.bai'
    params:
        bwa_version = "0.7.17",
        bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_AlignToMt.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_AlignToMt.stdout"

use rule AlignAndMarkDuplicates from AlignmentPipeline as AlignToShiftedMt with:
    input:
        unmapped_bam = rules.RevertSam.output.unmapped_bam,
        mt_ref_fasta = mt_shifted_fasta
    output:
        fastq = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_unmapped.fastq'),
        aln_sam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_aln.sam'),
        mba_bam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_mba.bam'),
        md_bam = temp('{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_md.bam'),
        metrics_filename = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted.duplicate_metrics',
        sorted_bam = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_sorted.bam',
        sorted_bai = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_shifted_sorted.bai'
    params:
        bwa_version = "0.7.17",
        bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_AlignToShiftedMt.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_AlignToShiftedMt.stdout"

rule CollectWgsMetrics:
    input:
        sorted_bam = rules.AlignToMt.output.sorted_bam,
        sorted_bam_index = rules.AlignToMt.output.sorted_bai,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
    output:
        meterics = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_metrics.txt',
        sens = '{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_theoretical_sensitivity.txt'
    params:
        read_length = max_read_length,
        coverage_cap = coverage_cap,
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CollectWgsMetrics.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CollectWgsMetrics.stdout"
    shell:
        r"""
        gatk CollectWgsMetrics \
          -INPUT {input.sorted_bam} \
          -VALIDATION_STRINGENCY SILENT \
          -REFERENCE_SEQUENCE {input.mt_ref_fasta} \
          -OUTPUT {output.meterics} \
          -USE_FAST_ALGORITHM true \
          -READ_LENGTH {params.read_length} \
          -COVERAGE_CAP {params.coverage_cap} \
          -INCLUDE_BQ_HISTOGRAM true \
          -THEORETICAL_SENSITIVITY_OUTPUT {output.sens} \
          2> {log.stderr} 1> {log.stdout}
        """

rule CallMt:
    input:
        sorted_bam = rules.AlignToMt.output.sorted_bam,
        sorted_bai = rules.AlignToMt.output.sorted_bai,
        wgsmetrics = rules.CollectWgsMetrics.output.meterics,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz"),
        vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz.tbi"),
        stats = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.vcf.gz.stats"),
        bamOut = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.bam")),
        baiOut = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt.bai"))
    params:
        max_reads_per_alignment_start = max_reads_per_alignment_start,
        m2_extra_args = m2_extra_args,
        region = "chrM:576-16024"
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CallMt.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CallMt.stdout"
    shell:
        r"""
        gatk Mutect2 \
            -R {input.mt_ref_fasta} \
            -I {input.sorted_bam} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {output.vcf} \
            --bam-output {output.bamOut} \
            {params.m2_extra_args} -L {params.region} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {params.max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            2> {log.stderr} 1> {log.stdout}
        """

rule CallShiftedMt:
    input:
        sorted_bam = rules.AlignToShiftedMt.output.sorted_bam,
        sorted_bai = rules.AlignToShiftedMt.output.sorted_bai,
        wgsmetrics = rules.CollectWgsMetrics.output.meterics,
        mt_ref_fasta = mt_shifted_fasta,
        mt_ref_fasta_index = mt_shifted_fasta_index,
        mt_dict = mt_shifted_dict,
    output:
        vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz"),
        vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz.tbi"),
        stats = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.vcf.gz.stats"),
        bamOut = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.bam")),
        baiOut = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_callMt_shifted.bai")),
    params:
        max_reads_per_alignment_start = max_reads_per_alignment_start,
        m2_extra_args = m2_extra_args,
        region = "chrM:8025-9144"
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CallShiftedMt.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CallShiftedMt.stdout"
    shell:
        r"""
        gatk Mutect2 \
            -R {input.mt_ref_fasta} \
            -I {input.sorted_bam} \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {output.vcf} \
            --bam-output {output.bamOut} \
            {params.m2_extra_args} -L {params.region} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {params.max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            2> {log.stderr} 1> {log.stdout}
        """


rule LiftoverAndCombineVcfs:
    input:
        shifted_vcf = rules.CallShiftedMt.output.vcf,
        vcf = rules.CallMt.output.vcf,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        shift_back_chain = shift_back_chain
    output:
        shifted_back_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.shifted_back.vcf"),
        shifted_back_vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.shifted_back.vcf.idx"),
        rejected_vcf = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.rejected.vcf",
        merged_vcf = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.merged.vcf",
        merged_vcf_index = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}.merged.vcf.idx"
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_LiftoverAndCombineVcfs.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_LiftoverAndCombineVcfs.stdout"
    shell:
        r"""
        gatk LiftoverVcf \
          -I {input.shifted_vcf} \
          -O {output.shifted_back_vcf} \
          -R {input.mt_ref_fasta} \
          -CHAIN {input.shift_back_chain} \
          -REJECT {output.rejected_vcf}

        gatk MergeVcfs \
          -I {output.shifted_back_vcf} \
          -I {input.vcf} \
          -O {output.merged_vcf} \
          2> {log.stderr} 1> {log.stdout}
        """

rule MergeStats:
    input:
        shifted_stats = rules.CallMt.output.stats,
        non_shifted_stats = rules.CallShiftedMt.output.stats
    output:
        "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_raw_combined.stats"
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_MergeStats.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_MergeStats.stdout"
    shell:
        r"""
        gatk MergeMutectStats --stats {input.shifted_stats} --stats {input.non_shifted_stats} -O {output} \
        2> {log.stderr} 1> {log.stdout}
        """

rule InitialFilter:
    input:
        raw_vcf = rules.LiftoverAndCombineVcfs.output.merged_vcf,
        raw_vcf_index = rules.LiftoverAndCombineVcfs.output.merged_vcf_index,
        raw_vcf_stats = rules.MergeStats.output,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        blacklisted_sites = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index
    output:
        FilterMutectCalls_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterMutectCalls.vcf"),
        FilterMutectCalls_vcf_index = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterMutectCalls.vcf.idx"),
        filtered_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.vcf.gz"),
        filtered_vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.vcf.gz.tbi"),
        bamout = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_InitialFilter.bam")) # We need to create these files regardless, even if they stay empty
    params:
        m2_extra_filtering_args = m2_extra_filtering_args,
        max_alt_allele_count = max_alt_allele_count,
        vaf_filter_threshold = vaf_filter_threshold,
        f_score_beta = f_score_beta,
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_InitialFilter.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_InitialFilter.stdout"
    shell:
        r"""
        gatk  FilterMutectCalls \
            -V {input.raw_vcf} \
            -R {input.mt_ref_fasta} \
            -O {output.FilterMutectCalls_vcf} \
            --stats {input.raw_vcf_stats} \
            {params.m2_extra_filtering_args} \
            --max-alt-allele-count {params.max_alt_allele_count} \
            --mitochondria-mode \
            --min-allele-fraction {params.vaf_filter_threshold} \
            --f-score-beta {params.f_score_beta}

        gatk VariantFiltration \
            -V {output.FilterMutectCalls_vcf} \
            -O {output.filtered_vcf} \
            --apply-allele-specific-filters \
            --mask {input.blacklisted_sites} \
            --mask-name "blacklisted_site" \
            2> {log.stderr} 1> {log.stdout}
        """

rule SplitMultiAllelicsAndRemoveNonPassSites:
    input:
        filtered_vcf = rules.InitialFilter.output.filtered_vcf,
        filtered_vcf_idx = rules.InitialFilter.output.filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        split_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_split.vcf"),
        split_vcf_index = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_split.vcf.idx"),
        vcf_for_haplochecker = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_splitAndPassOnly.vcf"),
        vcf_for_haplochecker_index = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_splitAndPassOnly.vcf.idx"),
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SplitMultiAllelicsAndRemoveNonPassSites.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SplitMultiAllelicsAndRemoveNonPassSites.stdout"
    shell:
        r"""
        gatk LeftAlignAndTrimVariants \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.split_vcf} \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac

          gatk SelectVariants \
            -V {output.split_vcf} \
            -O {output.vcf_for_haplochecker} \
            --exclude-filtered \
            2> {log.stderr} 1> {log.stdout}
        """

rule GetContamination:
    input:
        filtered_vcf = rules.SplitMultiAllelicsAndRemoveNonPassSites.output.vcf_for_haplochecker,
    output:
        contamination_tmp = temp("{outdir}/samples/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_contamination.txt"),
        contamination_raw_tmp = temp("{outdir}/samples/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_contamination.raw.txt"),
        contamination = "{outdir}/samples/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.txt",
        contamination_raw = "{outdir}/samples/{SAMPLE_IDS}/haplocheck/{SAMPLE_IDS}_chrM_contamination.raw.txt",
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_GetContamination.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_GetContamination.stdout"
    shell:
        r"""
        java -jar src/haplocheck/haplocheck.jar --raw --out {output.contamination_tmp} {input}
        sed 's/\"//g' {output.contamination_raw_tmp} > {output.contamination_raw}
        sed 's/\"//g' {output.contamination_tmp} > {output.contamination}
        2> {log.stderr} 1> {log.stdout}
        """

rule FilterContamination:
    input:
        filtered_vcf = rules.InitialFilter.output.filtered_vcf,
        filtered_vcf_idx = rules.InitialFilter.output.filtered_vcf_idx,
        raw_vcf_stats = rules.MergeStats.output,
        contamination_raw = rules.GetContamination.output.contamination_raw,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        blacklisted_sites = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index,
    output:
        FilterContaminationMutectCalls_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContaminationMutectCalls.vcf"),
        FilterContaminationMutectCalls_vcf_index = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContaminationMutectCalls.vcf.idx"),
        filtered_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.vcf.gz"),
        filtered_vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.vcf.gz.tbi"),
        bamout = temp(touch("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterContamination.bam"))
    params:
        f_score_beta = f_score_beta,
        m2_extra_filtering_args = m2_extra_filtering_args,
        max_alt_allele_count = max_alt_allele_count,
        vaf_filter_threshold = vaf_filter_threshold,
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterContamination.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterContamination.stdout"
    shell:
        r"""
        hc_contamination=$(awk -F "\t" 'FNR == 2' {input.contamination_raw} | awk '/YES/{{print $3;next}}{{print "0"}}')

        gatk  FilterMutectCalls \
            -V {input.filtered_vcf} \
            -R {input.mt_ref_fasta} \
            -O {output.FilterContaminationMutectCalls_vcf} \
            --stats {input.raw_vcf_stats} \
            {params.m2_extra_filtering_args} \
            --max-alt-allele-count {params.max_alt_allele_count} \
            --mitochondria-mode \
            --f-score-beta {params.f_score_beta}
            --contamination-estimate $hc_contamination

        gatk VariantFiltration \
            -V {output.FilterContaminationMutectCalls_vcf} \
            -O {output.filtered_vcf} \
            --apply-allele-specific-filters \
            --mask {input.blacklisted_sites} \
            --mask-name "blacklisted_site" \
            2> {log.stderr} 1> {log.stdout}
        """

rule FilterNuMTs:
    input:
        filtered_vcf = rules.FilterContamination.output.filtered_vcf,
        filtered_vcf_idx = rules.FilterContamination.output.filtered_vcf_idx,
        autosomal_coverage = rules.mtdnacn.output,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        numt_filtered_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterNuMTs.vcf.gz"),
        numt_filtered_vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_FilterNuMTs.vcf.gz.tbi"),
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterNuMTs.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterNuMTs.stdout"
    shell:
        r"""
        auto_cov=$(awk '{{ if ($2 == "median") {{ print $4}} }}' {input.autosomal_coverage})
        echo $auto_cov
        gatk NuMTFilterTool \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.numt_filtered_vcf} \
          --autosomal-coverage $auto_cov \
          2> {log.stderr} 1> {log.stdout}
        """

rule FilterLowHetSites:
    input:
        filtered_vcf = rules.FilterNuMTs.output.numt_filtered_vcf,
        filtered_vcf_idx = rules.FilterNuMTs.output.numt_filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        final_filtered_vcf = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.vcf.gz"),
        final_filtered_vcf_idx = temp("{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.vcf.gz.tbi"),
    params:
        max_sites = max_low_het_sites,
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterLowHetSites.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_FilterLowHetSites.stdout"
    shell:
        r"""
        gatk MTLowHeteroplasmyFilterTool \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.final_filtered_vcf} \
          --max-allowed-low-hets {params.max_sites} \
          2> {log.stderr} 1> {log.stdout}
        """

rule SplitMultiAllelicSites:
    input:
        filtered_vcf = rules.FilterLowHetSites.output.final_filtered_vcf,
        filtered_vcf_idx = rules.FilterLowHetSites.output.final_filtered_vcf_idx,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
    output:
        split_vcf = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz",
        split_vcf_idx = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM.final.split.vcf.gz.tbi",
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SplitMultiAllelicSites.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_SplitMultiAllelicSites.stdout"
    shell:
        r"""
        gatk LeftAlignAndTrimVariants \
          -R {input.mt_ref_fasta} \
          -V {input.filtered_vcf} \
          -O {output.split_vcf} \
          --split-multi-allelics \
          --dont-trim-alleles \
          --keep-original-ac \
          2> {log.stderr} 1> {log.stdout}
        """

rule CoverageAtEveryBase:
    input:
        bam_regular_ref = rules.AlignToMt.output.sorted_bam,
        bam_regular_ref_index = rules.AlignToMt.output.sorted_bai,
        bam_shifted_ref = rules.AlignToShiftedMt.output.sorted_bam,
        bam_shifted_ref_index = rules.AlignToShiftedMt.output.sorted_bai,
        shift_back_chain = shift_back_chain,
        control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
        non_control_region_interval_list = non_control_region_interval_list,
        mt_ref_fasta = mt_fasta,
        mt_ref_fasta_index = mt_fasta_index,
        mt_dict = mt_dict,
        shifted_mt_ref_fasta = mt_shifted_fasta,
        shifted_mt_ref_fasta_index = mt_shifted_fasta_index,
        shifted_mt_dict = mt_shifted_dict,
    output:
        table = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_per_base_coverage.tsv",
        pbc_non_control_region_tsv = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_non_control_region.tsv",
        pbc_non_control_region_meterics = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_non_control_region.metrics",
        pbc_control_region_shifted_tsv = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_control_region_shifted.tsv",
        pbc_control_region_shifted_meterics = "{outdir}/samples/{SAMPLE_IDS}/gatk/{SAMPLE_IDS}_chrM_control_region_shifted.metrics"
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CoverageAtEveryBase.stderr",
        stdout = "{outdir}/samples/{SAMPLE_IDS}/logs/{SAMPLE_IDS}_CoverageAtEveryBase.stdout"
    shell:
        r"""
        gatk CollectHsMetrics \
          -I {input.bam_regular_ref} \
          -R {input.mt_ref_fasta} \
          -PER_BASE_COVERAGE {output.pbc_non_control_region_tsv} \
          -O {output.pbc_non_control_region_meterics}\
          -TI {input.non_control_region_interval_list} \
          -BI {input.non_control_region_interval_list} \
          -covMax 20000 \
          -SAMPLE_SIZE 1

        gatk CollectHsMetrics \
          -I {input.bam_shifted_ref} \
          -R {input.shifted_mt_ref_fasta} \
          -PER_BASE_COVERAGE {output.pbc_control_region_shifted_tsv} \
          -O {output.pbc_control_region_shifted_meterics} \
          -TI {input.control_region_shifted_reference_interval_list} \
          -BI {input.control_region_shifted_reference_interval_list} \
          -covMax 20000 \
          -SAMPLE_SIZE 1

        Rscript workflow/scripts/CoverageAtEveryBase.R {output.pbc_non_control_region_tsv} {output.pbc_control_region_shifted_tsv} {output.table}

        2> {log.stderr} 1> {log.stdout}
        """
