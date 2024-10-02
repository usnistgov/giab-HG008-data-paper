"""
File: get_bam_stat.py
Author: Noah Spies
Usage: python get_bam_stat.py input.bam output.tsv.gz
Description of outputs from chatGPT4: 

reads a BAM file  and generates statistics about the reads. Explanations of the outputs summarized by chatGPT:

aln_lengthmax: The maximum alignment length for each read. This represents the longest segment of the read that aligns to the reference genome.
aln_count: The number of alignments for each read. This counts how many times a read is aligned to different regions of the reference genome.
aln_lengthsum: The total sum of alignment lengths for each read. This is the cumulative length of all aligned segments for a read.
ref_lengthsum: The total sum of reference lengths covered by each read. This is the cumulative length of the reference genome covered by the read's alignments.
ref_lengthmax: The maximum reference length covered by a single alignment of the read. This represents the longest segment of the reference genome that a read aligns to in one go.
ref_lengthcount: The number of distinct reference alignments for each read. This is similar to aln_count but specific to the reference segments.
bases: The sum of bases for each read, excluding supplementary alignments. This represents the total length of the read sequences that are considered for primary alignments.
"""

import sys
import pandas
import pysam

def get_bam_stats(bam_path, stats_path):
    
    aln_lengths = []
    ref_lengths = []

    bases = []
    read_ids = []

    try:
        bam = pysam.AlignmentFile(bam_path)
    except ValueError:
        print(f"Error reading from {bam_path}...")
        raise

    for read in bam:
        if read.is_unmapped or read.is_secondary: continue
        
        read_ids.append(read.query_name)        
        aln_lengths.append(read.query_alignment_length)
        ref_lengths.append(read.reference_length)

        if not read.is_supplementary:
            bases.append(read.query_length)
        else:
            bases.append(None)
            
    result = pandas.DataFrame({"read_id":read_ids,
                               "aln_length":aln_lengths, 
                               "ref_length":ref_lengths, 
                               "bases":bases})
    result = result.groupby("read_id").agg({"aln_length":[sum, max, "count"], 
                                            "ref_length":[sum, max, "count"], 
                                            "bases":sum})
            
    result.columns = pandas.Series(result.columns.tolist()).apply(pandas.Series).sum(axis=1)
    result = result.rename(columns={"aln_lengthcount":"aln_count", 
                                    "basessum":"bases"})
    
    result.to_csv(stats_path, sep="\t", compression="gzip")


def main():
    if len(sys.argv) == 3:
        bam_path = sys.argv[1]
        out_path = sys.argv[2]
    else:
        print("[ERROR] bam and output path must be provided")


    print("Getting read length stats ...")
    get_bam_stats(bam_path, out_path)

if __name__ == '__main__':
    main()
    
    