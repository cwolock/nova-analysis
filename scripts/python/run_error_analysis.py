#!/usr/bin/env python
"""Run novaseq error analysis on a list of files
Input is the original bam file location
"""
import argparse
import subprocess

def generate_file(file_list):
    with open(file_list, 'r') as infile:
        for line in infile:
            line = line.strip()
            prefix = '/nfs/seqscratch11/cw3026/NOVA_Analysis/6.13/' + line.split(
                                                                '/')[-1].split('.')[0]
            in_bam = prefix + '.ch21.calmd.bam'
            out_file = prefix + '.ch21.error.txt'
            subprocess.call(
                '/nfs/goldstein/software/python2.7.7/bin/python ~/github/nova-analysis/scripts/python/error_finder_simple_6.7.18.py --bam {} --fasta /nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5_BWAmem/hs37d5.fa --error_file {}'.format(
                                                    in_bam, out_file), shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--file_list', help='Specify name of input file list')
    args = parser.parse_args()
    generate_file(args.file_list)
