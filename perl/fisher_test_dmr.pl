#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX;

#perl fisher_test.hyp.pl --t linc3-7 --c_type CHG

my ($treatment ,$c_type, $t,$ty);
GetOptions( 
		"treatment=s" => \$treatment,
		"c_type=s" => \$c_type, 
);
 
$t=$treatment;   ###
$ty=$c_type;   

open OF , ">$t\_$ty\_fisher_test.R" || die $!;                      ###stringsAsFactors = F
print OF "data <- read.table(\"$t\_$ty\_methyl_level_regions.txt\",sep =\"\\t\", header = T ,stringsAsFactors = F)\np<-c()\nn<-length(data\$Region)\n";
print OF "for (i in 1:n) {\n\tf<-fisher.test(rbind(c(data[i,3],data[i,6]) , c(data[i,4],data[i,7])))\n\tp[i]<-f\$p\.value\n}\n";                                                        
print OF "q <- p.adjust(p , method = \"BH\" , n=length(p))\np_q <- cbind(p , q , data\$Control , data\$Treatment , data\$DMR_type ,data\$S,data\$E)\nwrite.table(p_q , file = \"$t\_$ty\_fisher_test_results.txt\" , sep = \"\\t\" , quote = F , row.names = data\$Region , col.names = F)";
close OF;                                                                              ###data\$DMR_type

system("R < $t\_$ty\_fisher_test.R --vanilla");
