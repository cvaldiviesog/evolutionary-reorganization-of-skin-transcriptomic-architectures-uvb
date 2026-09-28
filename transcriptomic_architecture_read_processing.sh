#!/usr/bin/env bash

# ============================================================================
# Read processing and mapping pipeline
# Paper: Evolutionary reorganization of transcriptomic architectures under UVB
#
# This script consolidates the commands contained in the original shell
# scripts supplied with the project.
#
# Samples:
#   OC1, OC2, OC3   = control
#   OUV1, OUV2, OUV3 = UVB-exposed
#
# NOTE:
# The published Methods additionally mention Cutadapt, SAMtools, Fastq-Pair,
# and featureCounts. Explicit commands for those steps were not present in
# the supplied shell scripts, so they are NOT invented here.
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Decompress raw reads
# ----------------------------------------------------------------------------

echo "Decompressing files"

# Original script contained this command commented out:
# gunzip -d *.fastq.gz

echo "All files are decompressed"


# ----------------------------------------------------------------------------
# 2. Quality filtering with fastp
# ----------------------------------------------------------------------------

echo "Starting read filtering with fastp"

echo "Working on OC1"
fastp -q 20 -w 4 \
    -i OC1_1.fastq \
    -I OC1_2.fastq \
    -o OC1_1.filtered.fastq \
    -O OC1_2.filtered.fastq

echo "Working on OC2"
fastp -q 20 -w 4 \
    -i OC2_1.fastq \
    -I OC2_2.fastq \
    -o OC2_1.filtered.fastq \
    -O OC2_2.filtered.fastq

echo "Working on OC3"
fastp -q 20 -w 4 \
    -i OC3_1.fastq \
    -I OC3_2.fastq \
    -o OC3_1.filtered.fastq \
    -O OC3_2.filtered.fastq

echo "Working on OUV1"
fastp -q 20 -w 4 \
    -i OUV1_1.fastq \
    -I OUV1_2.fastq \
    -o OUV1_1.filtered.fastq \
    -O OUV1_2.filtered.fastq

echo "Working on OUV2"
fastp -q 20 -w 4 \
    -i OUV2_1.fastq \
    -I OUV2_2.fastq \
    -o OUV2_1.filtered.fastq \
    -O OUV2_2.filtered.fastq

echo "Working on OUV3"
fastp -q 20 -w 4 \
    -i OUV3_1.fastq \
    -I OUV3_2.fastq \
    -o OUV3_1.filtered.fastq \
    -O OUV3_2.filtered.fastq

echo "Read filtering done"


# ----------------------------------------------------------------------------
# 3. Adapter removal with Cutadapt
# ----------------------------------------------------------------------------
#
# Described in the published Methods:
#   Cutadapt v1.18, default parameters.
#
# No Cutadapt command was present in the supplied shell scripts.
# Add the original project command here if it is available.
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 4. Mapping reads to the reference genome with STAR
# ----------------------------------------------------------------------------

echo "Starting STAR mapping"

# The original script had the OC1 command commented out:
# echo "Working on OC1"
# STAR --runThreadN 5 \
#      --genomeDir . \
#      --outFileNamePrefix OC1 \
#      --outSAMtype BAM SortedByCoordinate \
#      --readFilesIn OC1_1.filtered.fastq OC1_2.filtered.fastq
# echo "OC1 done"

echo "Working on OC2"
STAR --runThreadN 5 \
     --genomeDir . \
     --outFileNamePrefix OC2 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn OC2_1.filtered.fastq OC2_2.filtered.fastq
echo "OC2 done"

echo "Working on OC3"
STAR --runThreadN 5 \
     --genomeDir . \
     --outFileNamePrefix OC3 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn OC3_1.filtered.fastq OC3_2.filtered.fastq
echo "OC3 done"

echo "Working on OUV1"
STAR --runThreadN 5 \
     --genomeDir . \
     --outFileNamePrefix OUV1 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn OUV1_1.filtered.fastq OUV1_2.filtered.fastq
echo "OUV1 done"

echo "Working on OUV2"
STAR --runThreadN 5 \
     --genomeDir . \
     --outFileNamePrefix OUV2 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn OUV2_1.filtered.fastq OUV2_2.filtered.fastq
echo "OUV2 done"

echo "Working on OUV3"
STAR --runThreadN 5 \
     --genomeDir . \
     --outFileNamePrefix OUV3 \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn OUV3_1.filtered.fastq OUV3_2.filtered.fastq
echo "OUV3 done"


# ----------------------------------------------------------------------------
# 5. SAMtools filtering
# ----------------------------------------------------------------------------
#
# Published Methods:
#   -q 30 : mapping quality >= 30
#   -f 2  : properly paired reads
#
# No SAMtools filtering command was present in the supplied shell scripts.
# Add the original command here if available.
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 6. Synchronize paired-end reads with Fastq-Pair
# ----------------------------------------------------------------------------
#
# Published Methods:
#   Fastq-Pair v0.4, default parameters.
#
# No Fastq-Pair command was present in the supplied shell scripts.
# Add the original command here if available.
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 7. Gene-level read counting with featureCounts
# ----------------------------------------------------------------------------
#
# Published Methods:
#   -p : paired-end counting
#   -s 0 : unstranded libraries
#
# No featureCounts command was present in the supplied shell scripts.
# Add the original command here if available.
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 8. RSEM expression calculation (present in the supplied scripts)
# ----------------------------------------------------------------------------
#
# This section was present in the original "calc_ex.sh" file. It is kept here
# for provenance, but it is NOT described in the Methods paragraph above,
# which instead describes featureCounts for gene-level read counts.
#
# Reference preparation:
# ----------------------------------------------------------------------------

echo "Preparing RSEM reference"

perl RSEM/rsem-prepare-reference \
    -gff3 Orestias_ascotanensis.genome.gff3 \
    Orestias_ascotanensis.genome.fasta \
    Orestias


# ----------------------------------------------------------------------------
# RSEM expression calculation for OUV2
# ----------------------------------------------------------------------------
#
# This was the only sample explicitly present in calc_ex.sh.
# ----------------------------------------------------------------------------

echo "Calculating expression for OUV2"

perl RSEM/rsem-calculate-expression \
    -p 5 \
    --star \
    --paired-end \
    OUV2_1.filtered.fastq \
    OUV2_2.filtered.fastq \
    Orestias \
    OUV2

echo "Pipeline commands completed"
