# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:37:52 2018

@author: DXKE
"""

"""
input file format
Chr5	5	+	0	0	CHH	CCA
"""

"""

Usage: python CombineBisGenome.py [--fn1 <FN1>] [--fn2 <FN2>] [-o <outfile>]
"""

import sys

Dict1={}
Dict2={}

def CombineBisGenome(fn1,fn2,out):
#def CombineBisGenome(fn1,fn2):
    f1=open(fn1,'r')
    for line1 in f1:
        chr1, pos1, strand1, me1, ume1, pattern1 ,trinuc1 = line1.strip().split()
        CP1=str(chr1)+'_'+str(pos1)
        Dict1[CP1]=[me1,ume1,pattern1]

    f2=open(fn2,'r')
    for line2 in f2:
        chr2, pos2, strand2, me2, ume2, pattern2 ,trinuc2 = line2.strip().split()
        CP2=str(chr2)+'_'+str(pos2)
        Dict2[CP2]=[me2,ume2]
    
    of=open(out,'w')
    #for CHR,DPOS in Dict1.items():
       # for POS in DPOS.keys():
            #if CHR in Dict2.keys:
               # if POS in Dict2.keys.keys:
    for CHRPOS in Dict1.keys():
        if CHRPOS in Dict2:
            CHR,POS=CHRPOS.split('_')
            patt=str(Dict1[CHRPOS][2])
            m1=int(Dict1[CHRPOS][0])
            um1=int(Dict1[CHRPOS][1])
            m2=int(Dict2[CHRPOS][0])
            um2=int(Dict2[CHRPOS][1])
            #print("%s\t%d\t%s\t%d\t%d\t%d\t%d" %(CHR,POS,Dict1[CHRPOS][2],Dict1[CHRPOS][0],Dict1[CHRPOS][1],Dict2[CHRPOS][0],Dict2[CHRPOS][1]),file=of)
            print("%s\t%s\t%s\t%d\t%d\t%d\t%d" %(CHR,POS,patt,m1,um1,m2,um2),file=of)

#CombineBisGenome(argv[1],argv[2],argv[3])



from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: python CombineBisGenome.py [--fn1 <FN1>] [--fn2 <FN2>] [-o <outfile>]\n"

    #
    parser = OptionParser(usage)
    parser.add_option("--fn1", dest="FN1", help="File1 name, bismark coverage2cytosine output ")
    parser.add_option("--fn2", dest="FN2", help="File2 name, bismark coverage2cytosine output ")
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified")
    (options, args) = parser.parse_args()
    #
    if (options.FN1 is None) :
        parser.print_help()
        exit(-1)
    #
    if (options.outfile is None) :
        sys.stdout = open(options.outfile, 'w')
        #
    #
    CombineBisGenome(options.FN1, options.FN2, options.outfile)
#

# ===========================================
if __name__ == "__main__":
    main()
