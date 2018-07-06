#!/usr/bin/env python
"""Generate list of vcf and bam files for rando"""
import argparse
import MySQLdb

def generate_file(sample_list, file_list):
    chgvids = []
    with open(sample_list, 'r') as infile:
        for line in infile:
            chgvids.append(line.strip())
    with open(file_list, 'w') as outfile:
        db = MySQLdb.connect(user='sequence_view',passwd='View123!',
                             db='sequenceDB',host='10.73.50.40')
        try:
            for chgvid in chgvids:
                cur = db.cursor()
                cur.execute("""
                            SELECT QC.AlignSeqFileLoc, QC.pseudo_prepid
                            FROM dragen_qc_metrics AS QC
                            INNER JOIN dragen_sample_metadata AS MT on 
                            MT.pseudo_prepid = QC.pseudo_prepid
                            WHERE MT.sample_name = "{}" 
                            """.format(chgvid))
                rows = cur.fetchall()
                for loc, pseudo in rows:
                    bam = '{loc}/{ch}.{ps}/{ch}.{ps}.realn.recal.bam'.format(
                        loc=loc,ch=chgvid,ps=pseudo)
                    outfile.write(bam+'\n')
        finally:
            if db.open:
                db.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sample_list', help='Specify name of sample list file')
    parser.add_argument('--file_list', help='Specify name of output file list')
    args = parser.parse_args()
    generate_file(args.sample_list, args.file_list)
