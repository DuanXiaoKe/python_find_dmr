# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:37:52 2018

@author: DXKE
"""

"""
Usage: python ChrBinFindDmr_v.py [--fn1 <FN1>] [--fn2 <FN2>] [--Ctype <CTYPE>] [--depth <DEPTH>] [--mCnum <mCNUM>] [--diff <DIFF>] [--pval <PVAL>] [-o <outfile>]\n"

"""

"""
input file1 format
Chr4	1001	CHH	0	0	0	0
Chr4	1002	CHH	0	0	0	0
Chr4	1004	CHG	1	3	0	0
"""

"""
input2 file2 format
Chr1    0   1000 
Chr1    500 1500
Chr1    1000    2000
Chr1    1500    2500
Chr1    2000    3000
"""
"""
#for keys in CHRPOS.keys():
        #print(keys,CHRPOS[keys],sep='\t')

Chr1_13334187	['CHG', 21, 8, 10, 3]
Chr1_13334261	['CHG', 9, 9, 2, 7]
Chr1_13334263	['CHG', 18, 22, 5, 5]
Chr1_13334368	['CHG', 23, 8, 11, 2]

"""

import sys
#import sys
import re
import scipy.stats as stats


#def StoreCombinefile(comfn,context,depth):
    #CHRPOS={}


def FindDmrs(comfn,chrbin,context,depth,cnum,diff,pval,outfn):
    CHRPOS={}
    IN1=open(comfn,'r')
    for line1 in IN1:
        chr1, pos, pattern, mC1, C1, mC2, C2 = line1.strip().split()
        context=str(context)
        pattern=str(pattern)
        depth=int(depth)
        mC1=int(mC1)
        C1=int(C1)
        mC2=int(mC2)
        C2=int(C2)
        if re.match(context,pattern) and mC1+C1>=depth and mC2+C2>=depth:  
            chrpos=str(chr1)+'_'+str(pos)
            CHRPOS[chrpos]=[pattern,mC1, C1, mC2, C2]
    #for keys in CHRPOS.keys():
        #print(keys,CHRPOS[keys],sep='\t')   #可以输出
    print('combine file read...done')
    
    
    
    cnum=int(cnum)
    diff=float(diff)
    pval=float(pval)
    
    i=1
    of=open(outfn,'w')
    IN2=open(chrbin,'r')
    for line2 in IN2:
         chrname,start,end = line2.strip().split()
         start=int(start)
         end=int(end)
         mc1=0
         c1=0
         mc2=0
         c2=0
         ccnum=0
         ilist=[]
         for inum in range(start,end+1):
             chrpos1=str(chrname)+'_'+str(inum)
             if chrpos1 in CHRPOS:
                 mc1+=int(CHRPOS[chrpos1][1])
                 c1+=int(CHRPOS[chrpos1][2])
                 mc2+=int(CHRPOS[chrpos1][3])
                 c2+=int(CHRPOS[chrpos1][4])
                 ccnum+=1
                 ilist.append(chrpos1)
                 #print(chrpos1)
         if mc1+c1>0 and mc2+c2>0 and ccnum>=cnum:
             meth1=mc1/(mc1+c1)
             meth2=mc2/(mc2+c2)
             oddsratio, pvalue = stats.fisher_exact([[mc1,c1],[mc2,c2]])
             #print(mc1,c1,mc2,c2,sep='\t')                  
             if pvalue <= pval and abs(meth1-meth2)>=diff:
                 sortedilist=sorted(ilist,key=lambda d : int(d.split('_')[-1]))
                 chrstpos=sortedilist[0]
                 chredpos=sortedilist[-1]
                 chr2=chrstpos.split('_')[0]
                 spos=int(chrstpos.split('_')[-1])
                 epos=int(chredpos.split('_')[-1])
                 print("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.3e" %(chr2,spos,epos,meth1,meth2,meth1-meth2,pvalue),file=of)
                 #print("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.3e" %(chr2,spos,epos,meth1,meth2,meth1-meth2,pvalue))
                 print("%d" %(i))
                 i+=1
            #else:  
              #continue
#CHRPOS=StoreCombinefile(argv[1],argv[2],argv[3])
#FindDmrs(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8])

from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: python ChrBinFindDmr_v.py [--fn1 <FN1>] [--fn2 <FN2>] [--Ctype <CTYPE>] [--depth <DEPTH>] [--mCnum <mCNUM>] [--diff <DIFF>] [--pval <PVAL>] [-o <outfile>]\n"

    #
    parser = OptionParser(usage)
    parser.add_option("--fn1", dest="FN1", help="File1 name, CombineBisGenome.py output ")
    parser.add_option("--fn2", dest="FN2", help="File2 name, ChrBin.py output ")
    parser.add_option("--Ctype", dest="CTYPE", help="C,CG,CHG or CHH, mC type ")
    parser.add_option("--depth", dest="DEPTH", help="base sequence depth ", type=int)
    parser.add_option("--mCnum", dest="mCNUM", help="mC number ", type=int)
    parser.add_option("--diff", dest="DIFF", help="difference: the absolute difference of methylation level of each bin", type=float)
    parser.add_option("--pval", dest="PVAL", help="pvalue: the pvalue of fisher_exact of each bin", type=float)
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
    FindDmrs(options.FN1, options.FN2, options.CTYPE, options.DEPTH, options.mCNUM, options.DIFF, options.PVAL, options.outfile)

# ===========================================
if __name__ == "__main__":
    main()
