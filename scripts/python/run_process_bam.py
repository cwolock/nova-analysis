#!/usr/bin/env python
"""Generate list of vcf and bam files for rando"""
import argparse
import subprocess

def generate_file(file_list):
    with open(file_list, 'r') as infile:
        for line in infile:
            line = line.strip()
            in_bam = line
            out_prefix = '/nfs/seqscratch11/cw3026/NOVA_Analysis/6.13/' + line.split(
                                                                '/')[-1].split('.')[0]
            subprocess.call(
                '~/github/nova-analysis/scripts/bash/process_bam.sh {} {}'.format(
                                                    in_bam, out_prefix), shell=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--file_list', help='Specify name of input file list')
    args = parser.parse_args()
    generate_file(args.file_list)
