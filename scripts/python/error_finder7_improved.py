#!/usr/bin/env python
"""Count "errors" by quality bin in MD-tagged BAM and count total bases"""
import argparse
import numpy as np
import pysam
import sys

def base_to_num(letter):
    if letter == 'A': number = 0
    elif letter == 'T': number = 1
    elif letter == 'G': number = 2
    elif letter == 'C': number = 3
    elif letter == 'N': number = 4
    return number

def rev_comp(b):
    if b == 0: b = 1
    elif b ==1: b = 0
    elif b ==2: b = 3
    elif b == 3: b = 2
    return b

def find_errors(bam, fasta, error_file, count_file):
    with pysam.AlignmentFile(bam, 'rb') as inbam:
        m = np.zeros(shape = (3,4,5))
        #m2 = np.zeros(shape = (3,4,5))
        #m3 = np.zeros(shape = (3,4,5))
        total = 0
        for col in inbam.pileup():
            l = []
            cov, mismatch = 0, 0
            for read in col.pileups:
                # only want high-quality mapping reads
                if (not read.is_del and not read.is_refskip and 
                        read.alignment.mapping_quality > 49):
                    cov += 1
                    if read.alignment.query_sequence[
                            read.query_position] != '=':
                        mismatch += 1
                        strand = 'fwd'
                        #base = read.alignment.query_sequence[
                            #read.query_position]
                        if read.alignment.is_reverse:
                         #   base = rev_comp(base)
                            strand = 'rev'
                        l.append((read.alignment.query_qualities[
                                read.query_position],
                            base_to_num(read.alignment.query_sequence[
                                read.query_position]),
                                strand))
            if cov == 0: continue
            # calculate error rate
            error = mismatch/float(cov)
            if cov > 19 and error <= 0.05:
                total += cov
            # make sure there was at least one mismatch read
            #if mismatch == 0: continue
            # only sites with high coverage and low error (not het)
            if cov > 19 and 0 < error <= 0.05:
                with pysam.FastaFile(fasta) as infasta:
                    ref = base_to_num(infasta.fetch(col.reference_name, 
                        col.pos, col.pos + 1))
                for r in l:
                    if r[2] == 'rev':
                        base_act = rev_comp(r[1])
                        ref_act = rev_comp(ref)
                    else: 
                        base_act = r[1]
                        ref_act = ref

                    if r[0] < 15:
                        m[0][ref_act][base_act] += 1
                    elif 15 <= r[0] < 29:
                        m[1][ref_act][base_act] += 1
                    elif r[0] >= 29:  
                        m[2][ref_act][base_act] += 1
    with open(error_file, 'w') as out:
        for slice_2d in m:
            np.savetxt(out, slice_2d, fmt='%1d', delimiter='\t')
    with open(count_file, 'w') as out2:
        out2.write('{}'.format(total))
    #total1, total2, total3 = np.sum(m1), np.sum(m2), np.sum(m3)
    #prop1, prop2, prop3 = m1/total1, m2/total2, m3/total3
                        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bam', help='Specify bam file')
    parser.add_argument('--fasta', help='Specify reference fasta')
    parser.add_argument('--error_file', help='Specify error ouput file')
    parser.add_argument('--count_file', help='Specify base count file')
    args = parser.parse_args()
    find_errors(args.bam, args.fasta, args.error_file, args.count_file)
