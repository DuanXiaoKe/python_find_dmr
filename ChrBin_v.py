# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:36:22 2018

@author: DXKE
"""

"""
Usage: python ChrBin_v.py [--chrlenfile <CHRLEN>] [--bin <BIN>] [--step <STEP>] [-o <outfile>]

"""

import sys

def chrbin(chrlenfn,binsize,stepsize,outfn):
    binsize=int(binsize)
    stepsize=int(stepsize)
    
    of=open(outfn,'w')
    
    IN2=open(chrlenfn,'r')
    for line2 in IN2:
         chrname, chrlen = line2.strip().split()
         chrlen=int(chrlen)
         chrlenstepsize=chrlen+stepsize
     
         startpos = 0
         startpos=int(startpos)
         for endpos in range( binsize,chrlenstepsize,stepsize ):
            if min( chrlen, endpos ) > startpos:
                zendpos=min( chrlen, endpos )
                print("%s\t%d\t%d" %(chrname,startpos,zendpos),file=of)
                                
            startpos += stepsize

#chrbin(argv[1],argv[2],argv[3],argv[4])



from optparse import OptionParser
#
# ===========================================
def main():
    usage = "Usage: python ChrBin_v.py [--chrlenfile <CHRLEN>] [--bin <BIN>] [--step <STEP>] [-o <outfile>]\n" 

    #
    parser = OptionParser(usage)
    parser.add_option("--chrlenfile", dest="CHRLEN", help="File name, chr length file ")
    parser.add_option("--bin", dest="BIN", help="binsize: fixed size of each bin", type=int)
    parser.add_option("--step", dest="STEP", help="stepsize: fixed size of each bin overlaped", type=int)
    parser.add_option("-o", dest="outfile", default=None,
                      help="To standard output if not specified")
    (options, args) = parser.parse_args()
    #
    if (options.CHRLEN is None) :
        parser.print_help()
        exit(-1)
    #
    if (options.outfile is None) :
        sys.stdout = open(options.outfile, 'w')
        #
    #
    chrbin(options.CHRLEN, options.BIN, options.STEP,options.outfile)
#

# ===========================================
if __name__ == "__main__":
    main()
