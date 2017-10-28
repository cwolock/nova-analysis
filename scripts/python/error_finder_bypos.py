#!/usr/bin/env python
"""Counts errors by read position, counts 'N's as errors"""
import numpy as np
import pysam
import sys


def switch_pos(qpos, platform):
    if platform == 'HiSeq':
        middle = np.median(range(0,126))
    elif platform == 'NovaSeq':
        middle = np.median(range(0,151))
    distance = qpos - middle
    new_pos = middle - distance
    return int(new_pos)

class FindErrors():
    def __init__(self, bam, outf, plat):
        self.bam = bam
        self.outf = outf
        self.plat = plat
    def run(self):
        if self.plat == 'HiSeq':
            bins = range(0,126)
        elif self.plat == 'NovaSeq':
            bins = range(0,151)
        errs = [0] * len(bins)
        totals = [0] * len(bins)
        with pysam.AlignmentFile(self.bam, 'rb') as inbam:
            for col in inbam.pileup():
                l = []
                cov, mismatch = 0, 0
                for read in col.pileups:
                    status = 0
                    # only want high-quality mapping reads
                    if (not read.is_del and not read.is_refskip and 
                            read.alignment.mapping_quality > 49):
                        cov += 1
                        if read.alignment.query_sequence[
                                read.query_position] != '=':
                            mismatch += 1
                            status = 1
                        if read.alignment.is_reverse:
                            pos = switch_pos(read.query_position, self.plat)
                            #print('{}\t{}\n'.format(read.query_position, pos))
                        else: pos = read.query_position      
                        l.append((status, pos))
                if cov == 0: continue
                # calculate error rate
                error = mismatch/float(cov)
                if cov > 19 and error <= 0.05:
                    for r in l:
                        totals[r[1]] += 1
                        if r[0] == 1:
                            errs[r[1]] += 1 

        with open(self.outf, 'w') as out:
            out.writelines(['{}\t{}\t{}\n'.format(
                str(bins[i]), str(errs[i]), str(totals[i])) 
                for i in range(len(bins))])
                        
if __name__ == '__main__':
    out = FindErrors(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]))
    out.run()
