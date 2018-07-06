#!/bin/bash

# get arguments
bam_in="$1"
out_prefix="$2"

# bedtools intersect to get ch21 only
/nfs/goldstein/software/bedtools-2.25.0/bin/bedtools intersect -abam $bam_in -b /nfs/seqscratch11/cw3026/bedfiles/chr21.bed > ${out_prefix}.ch21.bam 

# calmd
/nfs/goldstein/software/samtools-1.3/samtools calmd -b -e ${out_prefix}.ch21.bam /nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa > ${out_prefix}.ch21.calmd.bam

# index
/nfs/goldstein/software/samtools-1.3/samtools index ${out_prefix}.ch21.calmd.bam
