#!/usr/bin/bash
## Extracting longer than expected Sequel reads

PICARD=/scratch/nolson/sw/picard.jar
java -jar ${PICARD} FilterSamReads --version

## HG005
INBAM=../../HG002-data/HG005_PacBio_GRCh37.bam 
OUTBAM=HG005_PacBio_GRCh37_gt70kb.bam
READIDS=HG005_read_id_gt70k.tsv
java -Xmx32g -jar ${PICARD} FilterSamReads \
       I=${INBAM} \
       O=${OUTBAM} \
       READ_LIST_FILE=${READIDS} \
       FILTER=includeReadList

## HG006
INBAM=../../HG002-data/HG006_PacBio_GRCh37.bam 
OUTBAM=HG006_PacBio_GRCh37_gt70kb.bam
READIDS=HG006_read_id_gt70k.tsv
java -Xmx32g -jar ${PICARD} FilterSamReads \
       I=${INBAM} \
       O=${OUTBAM} \
       READ_LIST_FILE=${READIDS} \
       FILTER=includeReadList

## HG007
INBAM=../../HG002-data/HG007_PacBio_GRCh37.bam 
OUTBAM=HG007_PacBio_GRCh37_gt70kb.bam
READIDS=HG007_read_id_gt70k.tsv
java -Xmx32g -jar ${PICARD} FilterSamReads \
       I=${INBAM} \
       O=${OUTBAM} \
       READ_LIST_FILE=${READIDS} \
       FILTER=includeReadList


## Long read characterization
# Index bam for IGV
samtools index HG005_PacBio_GRCh37_gt70kb.bam
samtools index HG006_PacBio_GRCh37_gt70kb.bam
samtools index HG007_PacBio_GRCh37_gt70kb.bam

# Convert bam to sam - view in text editor
samtools view -h HG005_PacBio_GRCh37_gt70kb.bam > HG005_PacBio_GRCh37_gt70kb.sam
samtools view -h HG006_PacBio_GRCh37_gt70kb.bam > HG006_PacBio_GRCh37_gt70kb.sam
samtools view -h HG007_PacBio_GRCh37_gt70kb.bam > HG007_PacBio_GRCh37_gt70kb.sam

# Convert bam to fasta
samtools fasta HG005_PacBio_GRCh37_gt70kb.bam > HG005_PacBio_GRCh37_gt70kb.fasta
samtools fasta HG006_PacBio_GRCh37_gt70kb.bam > HG006_PacBio_GRCh37_gt70kb.fasta
samtools fasta HG007_PacBio_GRCh37_gt70kb.bam > HG007_PacBio_GRCh37_gt70kb.fasta

# Base composition analysis
fastacomposition -s -f HG005_PacBio_GRCh37_gt70kb.fasta > HG005_PacBio_GRCh37_gt70kb_basecomp.txt
fastacomposition -s -f HG006_PacBio_GRCh37_gt70kb.fasta > HG006_PacBio_GRCh37_gt70kb_basecomp.txt
fastacomposition -s -f HG007_PacBio_GRCh37_gt70kb.fasta > HG007_PacBio_GRCh37_gt70kb_basecomp.txt

# Ran repeatMasker on fasta files online