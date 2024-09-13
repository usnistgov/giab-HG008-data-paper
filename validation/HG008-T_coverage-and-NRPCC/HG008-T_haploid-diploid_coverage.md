# Haploid and diploid coverage and NRPCC for HG008-T
Document generated J.McDaniel 2024-09-13 to document commands and calcuations used by J.Zook and V.Patel to calculate HG008-T haploid and diploid coverage and NRPCC for short and long reads presented in the HG008 Scientific Data paper.

## Haploid Coverage
Because the tumor cell line is aneuploid, we also estimate the mean coverage of apparently diploid and haploid regions, as well as number of reads per tumor chromosomal copy (NRPCC). Single cell data showed that for 120 cells, chromosome 4 was haploid in almost all cells for HG008-T and was used as the basis for calculating haploid coverage. 

**1) 1 Mbp segments of GRCh38 chromosome 4 that don't contain any 
[segdups](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3[…]SegmentalDuplications/GRCh38_notinsegdups_gt10kb.bed.gz) 
were generated.**

`bedtools2.30.0/bin/intersectBed -u -f 1.0 -a chr4_1Mbp_segments_GRCh38.bed -b GRCh38_notinsegdups_gt10kb.bed.gz  > chr4_1Mbp_segments_GRCh38_notinanysegdups.bed`

**2) run [mosdepth](https://github.com/brentp/mosdepth) to generate regions.bed from input .bam for use in coverage calculation**

`mosdepth --by chr4_1Mbp_segments_GRCh38_notinanysegdups.bed <prefix for outputs> <input.bam> --threads 3`


**3) calculate chr4_haploid_mosdepth_mean_coverage, mean (haploid) mosdepth coverage**

`cat <mosdepth .regions.bed> | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$4 }END{print SUM/NR}`


## Diploid Coverage
Using inputs from haploid coverage calculation, multiply calculated mean by 2. 

**calculate chr4_diploid_mosdepth_mean_coverage, mean (diploid) mosdepth coverage**  
`cat <mosdepth .regions.bed> | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$4 }END{print (SUM/NR)*2}`


## NRPCC (number of reads per chromosome copy)
Here the assumption being made is that whole genome doubling has occurred and tumor purity is 1. The diploid and haploid coverages presented above, are estimated for cells without whole genome doubling, whereas NRPCC is the coverage of each copy of the chromosome for cells with whole genome doubling.

**1) calculate tumor ploidy from previously calculated coverages**  
Calculations use depth of autosomal coverage and haploid coverage. Autosome coverage was calculated using different methods as appropriate for short or long reads. 
* for short reads  
`(Mean coverage from GA4GH mean_autosome_coverage / chr4_haploid_mosdepth_mean_coverage) *2`
* for long reads  
`(Mean coverage from Cramino mean_coverage / chr4_haploid_mosdepth_mean_coverage) *2`

**2) NRPCC**  
The NRPCC calculation takes into account tumor ploidy, purity and depth of coverage.  
`chr4_haploid_mosdepth_mean_coverage / 2`
