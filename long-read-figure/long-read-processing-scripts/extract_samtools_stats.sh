#!/usr/bin/bash

BAMSTAT=$1
fbname=$(basename "$1" .txt)
SUMMARYSTATS=${fbname}_summary.txt
READLENGTH=${fbname}_read_lengths.txt
COVSTATS=${fbname}_cov.txt

grep "^SN" ${BAMSTAT} | cut -f 2- > extract_samtools_stats_output/${SUMMARYSTATS}
grep "^RL" ${BAMSTAT} | cut -f 2- > extract_samtools_stats_output/${READLENGTH}
grep "^COV" ${BAMSTAT} | cut -f 2- > extract_samtools_stats_output/${COVSTATS}