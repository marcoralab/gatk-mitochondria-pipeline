# Resources
External resources required to run mitchondrial pipeline.

## GATK Genome Reference files

Files are avaliable from the Broad References Public Data google buckets.


**hg38 Reference**

```
# https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0
gsutil cp \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.pac" \
  "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.sa" \
  .
```

**hg19 reference**

```
# https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0
gsutil cp \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.sa" \
  .
```

**chrM**

```
# https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0
gsutil cp \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.pac" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta.sa" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/README" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/ShiftBack.chain" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed.idx" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.bed.idx.old" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.bed.old" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/blacklist_sites.hg38.chrM.shifted_by_8000_bases.fixed.bed.idx" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/chrMWithFinalNuMTs.hg38.interval_list" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/control_region_shifted.chrM.interval_list" \
  "gs://gcp-public-data--broad-references/hg38/v0/chrM/non_control_region.chrM.interval_list" \
  .
```

## gnomAD resources

Files are avaliable from [broadinstitute/gnomad-mitochondria](https://github.com/broadinstitute/gnomad-mitochondria) and are made avaliable in `resources/gnomad/`

```
artifact_prone_sites.bed
transfer files from previous repo
mitotip_scores_08_27_2020.txt
pon_mt_trna_predictions_08_27_2020.txt
rCRS-centered_phylo_vars_final_update.txt
chrM_pos_ref_alt_context_categories.txt
```

dbSNP chrM sites .vcf file was generated from reference=GRCh38.p12 dbSNP_BUILD_ID=154.

```
wget https://ftp.ncbi.nih.gov/snp/archive/b154/VCF/GCF_000001405.38.gz
bcftools view GCF_000001405.38.vcf.gz --regions NC_012920.1 -Ou -o GCF_000001405.38.chrM.vcf
```

## Haplocheck

Haplocheck is required for estimating contamination estimates.

```
# TODO make docker file and use singularity
mkdir src/haplocheck
cd src/haplocheck
wget https://github.com/genepi/haplocheck/releases/download/v1.3.2/haplocheck.zip
unzip haplocheck.zip
./haplocheck --out <out-file> <input-vcf>
```

## Example data

Try the pipeline on a couple of samples from [1000 thousand genomes](https://www.internationalgenome.org/data-portal/population).


```
# NA12878
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.bai

# NA18916
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18916/alignment/NA18916.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18916/alignment/NA18916.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam.bai

```
