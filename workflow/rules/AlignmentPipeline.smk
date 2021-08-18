# snakejob -s MitochondriaPipeline.smk -j 10 --use-conda

import os

ID = ["test"]
mt_fasta = "resources/hg38/Homo_sapiens_assembly38.fasta"

rule all:
    input:
        expand('{IDS}_chrM_sorted.bam', IDS = ID),
        # expand("sandbox/output/{SAMPLE_IDS}_variantfiltered.vcf.gz", SAMPLE_IDS = ID),


rule AlignAndMarkDuplicates:
    input:
        unmapped_bam = "{IDS}_unmapped_bam",
        mt_ref_fasta = mt_fasta
    output:
        fastq = temp('{IDS}_chrM_unmapped.fastq'),
        aln_sam = temp('{IDS}_chrM_aln.sam'),
        mba_bam = temp('{IDS}_chrM_mba.bam'),
        md_bam = temp('{IDS}_chrM_md.bam'),
        duplicate_metrics = '{IDS}_chrM.duplicate_metrics',
        sorted_bam = '{IDS}_chrM_sorted.bam'
    params:
        bwa_version = "0.7.17",
        bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",
    conda: '../envs/gatk.yaml'
    log:
        stderr = "{IDS}_AlignAndMarkDuplicates.stderr",
        stdout = "{IDS}_AlignAndMarkDuplicates.stdout"
    shell:
        r"""
        gatk SamToFastq \
          -INPUT {input.unmapped_bam} \
          -FASTQ {output.fastq} \
          -INTERLEAVE true \
          -NON_PF true

        {params.bwa_commandline} {input.mt_ref_fasta} {output.fastq} > {output.aln_sam}

        gatk MergeBamAlignment \
          -VALIDATION_STRINGENCY SILENT \
          -EXPECTED_ORIENTATIONS FR \
          -ATTRIBUTES_TO_RETAIN X0 \
          -ATTRIBUTES_TO_REMOVE NM \
          -ATTRIBUTES_TO_REMOVE MD \
          -ALIGNED_BAM {output.aln_sam} \
          -UNMAPPED_BAM {input.unmapped_bam} \
          -OUTPUT {output.mba_bam} \
          -REFERENCE_SEQUENCE {input.mt_ref_fasta} \
          -PAIRED_RUN true \
          -SORT_ORDER "unsorted" \
          -IS_BISULFITE_SEQUENCE false \
          -ALIGNED_READS_ONLY false \
          -CLIP_ADAPTERS false \
          -MAX_RECORDS_IN_RAM 2000000 \
          -ADD_MATE_CIGAR true \
          -MAX_INSERTIONS_OR_DELETIONS -1 \
          -PRIMARY_ALIGNMENT_STRATEGY MostDistant \
          -PROGRAM_RECORD_ID "bwamem" \
          -PROGRAM_GROUP_VERSION "{params.bwa_version}" \
          -PROGRAM_GROUP_COMMAND_LINE "{params.bwa_commandline} {input.mt_ref_fasta} {output.fastq} > {output.aln_sam}" \
          -PROGRAM_GROUP_NAME "bwamem" \
          -UNMAPPED_READ_STRATEGY COPY_TO_TAG \
          -ALIGNER_PROPER_PAIR_FLAGS true \
          -UNMAP_CONTAMINANT_READS true \
          -ADD_PG_TAG_TO_READS false

        gatk MarkDuplicates \
          -INPUT {output.mba_bam} \
          -OUTPUT {output.md_bam} \
          -METRICS_FILE {output.metrics_filename} \
          -VALIDATION_STRINGENCY SILENT \
          -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
          -ASSUME_SORT_ORDER "queryname" \
          -CLEAR_DT "false" \
          -ADD_PG_TAG_TO_READS false

        gatk SortSam \
          -INPUT {output.md_bam} \
          -OUTPUT {output.sorted_bam} \
          -SORT_ORDER "coordinate" \
          -CREATE_INDEX true \
          -MAX_RECORDS_IN_RAM 300000

          2> {log.stderr} 1> {log.stdout}
        """

# rule AlignAndMarkDuplicates_SamToFastq:
#     input: rules.RevertSam.output.unmapped_bam,
#     output:
#         'sandbox/output/{SAMPLE_IDS}_chrM_unmapped.fastq'
#     conda: 'envs/gatk.yaml'
#     shell:
#         r"""
#         gatk SamToFastq \
#           -INPUT {input} \
#           -FASTQ {output} \
#           -INTERLEAVE true \
#           -NON_PF true
#         """
#
# rule AlignAndMarkDuplicates_BWA:
#     input:
#         fastq = rules.AlignAndMarkDuplicates_SamToFastq.output,
#         mt_ref_fasta = 'raw/reference/chrM/Homo_sapiens_assembly38.chrM.fasta',
#     output:
#         'sandbox/output/{SAMPLE_IDS}_chrM_aln.sam',
#     conda: 'envs/bwa.yaml'
#     shell:
#         r"""
#         bwa mem -K 100000000 -p -v 3 -t 2 -Y {input.mt_ref_fasta} {input.fastq} > {output}
#         """
#
# rule AlignAndMarkDuplicates_MergeBamAlignment:
#     input:
#         aln_sam = rules.AlignAndMarkDuplicates_BWA_mem.output[0],
#         unmapped_bam = rules.RevertSam.output.unmapped_bam,
#         mt_ref_fasta = 'raw/reference/chrM/Homo_sapiens_assembly38.chrM.fasta'
#     output: 'sandbox/output/{SAMPLE_IDS}_chrM_mba.bam'
#     params:
#         bwa_version = "0.7.17",
#         bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y",
#         mt_ref_fasta = 'raw/reference/chrM/Homo_sapiens_assembly38.chrM.fasta',
#         fastq = rules.AlignAndMarkDuplicates_SamToFastq.output,
#     conda: 'envs/gatk.yaml'
#     shell:
#         r"""
#         gatk MergeBamAlignment \
#           -VALIDATION_STRINGENCY SILENT \
#           -EXPECTED_ORIENTATIONS FR \
#           -ATTRIBUTES_TO_RETAIN X0 \
#           -ATTRIBUTES_TO_REMOVE NM \
#           -ATTRIBUTES_TO_REMOVE MD \
#           -ALIGNED_BAM {input.aln_sam} \
#           -UNMAPPED_BAM {input.unmapped_bam} \
#           -OUTPUT {output} \
#           -REFERENCE_SEQUENCE {input.mt_ref_fasta} \
#           -PAIRED_RUN true \
#           -SORT_ORDER "unsorted" \
#           -IS_BISULFITE_SEQUENCE false \
#           -ALIGNED_READS_ONLY false \
#           -CLIP_ADAPTERS false \
#           -MAX_RECORDS_IN_RAM 2000000 \
#           -ADD_MATE_CIGAR true \
#           -MAX_INSERTIONS_OR_DELETIONS -1 \
#           -PRIMARY_ALIGNMENT_STRATEGY MostDistant \
#           -PROGRAM_RECORD_ID "bwamem" \
#           -PROGRAM_GROUP_VERSION "{params.bwa_version}" \
#           -PROGRAM_GROUP_COMMAND_LINE "{params.bwa_commandline} {params.mt_ref_fasta} {params.fastq} > {input.aln_sam}" \
#           -PROGRAM_GROUP_NAME "bwamem" \
#           -UNMAPPED_READ_STRATEGY COPY_TO_TAG \
#           -ALIGNER_PROPER_PAIR_FLAGS true \
#           -UNMAP_CONTAMINANT_READS true \
#           -ADD_PG_TAG_TO_READS false
#         """
#
# rule AlignAndMarkDuplicates_MarkDuplicates:
#     input:
#         mba_bam = rules.AlignAndMarkDuplicates_MergeBamAlignment.output
#     output:
#         md_bam = 'sandbox/output/{SAMPLE_IDS}_chrM_md.bam',
#         metrics_filename = 'sandbox/output/{SAMPLE_IDS}_chrM.meterics'
#     params:
#         bwa_version = "0.7.17",
#         read_name_regex = "(?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$" # if this is used -READ_NAME_REGEX {params.read_name_regex} \
#     conda: 'envs/gatk.yaml'
#     shell:
#         r"""
#         gatk MarkDuplicates \
#           -INPUT {input.mba_bam} \
#           -OUTPUT {output.md_bam} \
#           -METRICS_FILE {output.metrics_filename} \
#           -VALIDATION_STRINGENCY SILENT \
#           -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
#           -ASSUME_SORT_ORDER "queryname" \
#           -CLEAR_DT "false" \
#           -ADD_PG_TAG_TO_READS false
#         """
#
# rule AlignAndMarkDuplicates_SortSam:
#     input: rules.AlignAndMarkDuplicates_MarkDuplicates.output.md_bam
#     output: 'sandbox/output/{SAMPLE_IDS}_chrM_sorted.bam'
#     conda: 'envs/gatk.yaml'
#     shell:
#         r"""
#         gatk SortSam \
#           -INPUT {input} \
#           -OUTPUT {output} \
#           -SORT_ORDER "coordinate" \
#           -CREATE_INDEX true \
#           -MAX_RECORDS_IN_RAM 300000
#         """

####################################
