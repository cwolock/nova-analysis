###### PolyG Analyzer
#
# input: FASTQ
#
# output: 

import gzip
import numpy as np
import sys

class TailAnalysis():
    def __init__(self, infastq, metrics, dist, outf10, outf20, outf30):
        self.infastq = infastq
        self.metrics = metrics
        self.dist = dist
        self.outf10 = outf10
        self.outf20 = outf20
        self.outf30 = outf30
    def run(self):
        with gzip.open(self.outf10, 'wb') as out10, \
             gzip.open(self.outf20, 'wb') as out20, \
             gzip.open(self.outf30, 'wb') as out30, \
             gzip.open(self.infastq, 'rb') as infile:
            polyG10, polyG20, polyG30, reads = 0, 0, 0, 0
            d = {}
            lines = enumerate(infile)
            for i, line in lines:
                reads += 1
                name = line
                seq = lines.next()[1]
                plus = lines.next()[1]
                quals = lines.next()[1]
                tail = 0
                for j in range(len(seq)-2, -1, -1):
                    if seq[j] != 'G':
                        break
                    else:
                        tail += 1
                if tail not in d:
                    d[tail] = 1
                else:
                    d[tail] += 1    
                if tail > 9: 
                    polyG10 += 1
                    out10.write(''.join([name, seq, plus, quals]))
                if tail > 19: 
                    polyG20 += 1
                    out20.write(''.join([name, seq, plus, quals]))
                if tail > 29: 
                    polyG30 += 1
                    out30.write(''.join([name, seq, plus, quals]))
        with open(self.dist, 'w') as distfile:
            for k, v in d.iteritems():
                distfile.write('{}\t{}\n'.format(k, v))
        with open(self.metrics, 'w') as metfile:
            metfile.write('G10={},G20={},G30={},reads={},%10={},%20={},%30={}'.format(
                polyG10,polyG20,polyG30,reads,100*polyG10/np.float64(reads),
                100*polyG20/np.float64(reads),100*polyG30/np.float64(reads)))   
                    
if __name__ == '__main__':
    out = TailAnalysis(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]),
        str(sys.argv[4]), str(sys.argv[5]), str(sys.argv[6]))
    out.run()
