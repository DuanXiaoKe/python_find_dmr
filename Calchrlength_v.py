# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 11:30:28 2018

@author: DXKE

"""
"""
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
"""
"""
usage:python Calchrlength.py Chr1 Chr2 Chr3 Chr4 Chr5 Chr1.fna Chr2.fna Chr3.fna Chr4.fna Chr5.fna outfile 
"""


from sys import argv
from Bio import SeqIO

def Calchrlen(chr1,chr2,chr3,chr4,chr5,chr1len,chr2len,chr3len,chr4len,chr5len,outfn):
    of=open(outfn,'w')
    entries = [(chr1, chr1len),
           (chr2, chr2len),
           (chr3, chr3len),
           (chr4, chr4len),
           (chr5, chr5len)]
    for (name, filename) in entries:
        record = SeqIO.read(filename,"fasta")
        #print name, len(record)
        print("%s\t%d" % (name,len(record)),file=of)


Calchrlen(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],argv[10],argv[11])