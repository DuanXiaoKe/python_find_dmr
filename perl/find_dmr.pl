#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX;
#perl find_dmr_V.pl  --control col0_genome -t ros1-4_genome --depth 4 --c_type C --chr_size arabid_chr_length.txt -w 200 --step_size 50 -e 5  -p 0.05 -q 0.05 --diff_change 0.2
#$window_size usually is 2kb; depth >= 5; $chr_size is a file contains each chromosome's length;
my ($control , $treatment , $window_size , $c_type ,$step_size,$depth , $chr_size , $p_value , $q_value , $diff_change , $effective_c_num , $t,$ty);
my (%chr , %methyl_level , %methyl_region_avg ,%ha, @region , %methyl_c , %unmethyl_c ,@tmp) = ();
GetOptions( 
		"control=s" => \$control,
		"treatment=s" => \$treatment,
        "window_size=s" => \$window_size, 
		"step_size=s" => \$step_size,
        "diff_change=f" => \$diff_change, 
		"effective_c_num=i" => \$effective_c_num,
	#	"norm_reads=f{,}" => \@norm_reads,
        "depth=s" => \$depth, 
		"chr_size=s" => \$chr_size,
		"c_type=s" => \$c_type, 
		"p_value=f" => \$p_value,
		"q_value=f" => \$q_value,
); 

$t=(split/_/,$treatment)[0];   ###
$ty=$c_type;              ###
open CHR , $chr_size || die $!;	#chr size file;
open OF , ">$t\_$ty\_methyl_level_regions.txt" || die $!;

while (<CHR>) {
	chomp;
	my @l = split /\t/ , $_;
	$chr{$l[0]} = $l[1];
}

my @chr = sort {$a <=> $b} keys %chr;
my @file = ($control , $treatment);
foreach my $f (@file) {
	open IF , "$f" || die $!;	#file was generated from coverage2cytosine;
	#calculate the methylation level of sites which read coverage is no less than $depth;
	%methyl_level = ();     ####读入下一个文件前清零。
	while (<IF>) {
		my @tmp = split /\t/ , $_;
		if (($tmp[3] + $tmp[4]) >= $depth && $tmp[0] =~ /Chr(\d)/){
			my $chrnum=$1;
		if ($tmp[5] =~ /$c_type/){     ####筛选类别
			push @{$methyl_level{$chrnum}{$tmp[1]}} , ($tmp[3]/($tmp[3]+$tmp[4]),$tmp[3] , $tmp[4]);
		}
	}
	}
	close IF;
	#calculate the average methylation level and number of methylated Cs of each bin, and write it to a tmp file;
	foreach (@chr) {
		for (my $i=1; $i<=$chr{$_}; $i+=$step_size) {
			my $methyl_region_c = 0;
			#my $methyl_region_uc = 0;
			my $count = 0;
			my @COUNT =();
			my $a;
			my $bin = ceil($i / $step_size);
			last if $i + $step_size - 1 > $chr{$_};
			my $region = $_ . "_" . $bin;
			for (my $j=0; $j<$window_size; $j++) {
				if (defined $methyl_level{$_}{$j + $i}) {
					$methyl_region_c += $methyl_level{$_}{$j + $i}[0];
					#$methyl_region_uc += $methyl_level{$_}{$j + $i}[1];
					$count++;
					if ($methyl_level{$_}{$j + $i}[0] == 0) {
					$unmethyl_c{$region}{$f} += $methyl_level{$_}{$j + $i}[2];
					}else {
					$methyl_c{$region}{$f} += $methyl_level{$_}{$j + $i}[1];
					}
					my $pos = $j + $i;
					$a = $_ . "_" . $pos;
				    push @COUNT , $a;
				}
			}
			if ($count >= $effective_c_num && $methyl_region_c >=0.8) {
				#$methyl_region_avg{$region}{$f} = $methyl_region_c / ($methyl_region_c + $methyl_region_uc);
				$methyl_region_avg{$region}{$f} = $methyl_region_c;
				push @region , $region;
				my ($s,$e)=@COUNT[0,-1];
				push @{$ha{$region}{$f}},($s,$e);
			}
		}
	}
}
#%methyl_level = ();

my @nr_region;
my %count;
my $S;
my $E;
@nr_region = grep {++$count{$_} == 2} @region;
#print "@nr_region";                                                                                         #
print OF "Region\tControl\tMethyl_reads_C\tUnmethyl_reads_C\tTreatment\tMethyl_reads_C\tUnmethyl_reads_C\tDMR_type\tS\tE\n";
foreach (@nr_region) {
	my ($ucc , $cc , $uct , $ct);
	if (defined $unmethyl_c{$_}{$file[0]} && defined $unmethyl_c{$_}{$file[1]} && defined $methyl_c{$_}{$file[0]} && defined $methyl_c{$_}{$file[1]}) {
		$ucc = int $unmethyl_c{$_}{$file[0]};
		$cc = int $methyl_c{$_}{$file[0]};
		$uct = int $unmethyl_c{$_}{$file[1]};
		$ct = int $methyl_c{$_}{$file[1]};
	my ($m , $n , $o , $p);	
	$m=(split/_/,$ha{$_}{$file[0]}[0])[1];
	$n=(split/_/,$ha{$_}{$file[1]}[0])[1];
	$o=(split/_/,$ha{$_}{$file[0]}[1])[1];
	$p=(split/_/,$ha{$_}{$file[1]}[1])[1];
	if($m <= $n){$S=$m;}else{$S=$n;}
	if($o <= $p){$E=$p;}else{$E=$o;}
	if($methyl_region_avg{$_}{$file[0]}!=0 && $methyl_region_avg{$_}{$file[1]} == 0){
		print OF "$_\t$methyl_region_avg{$_}{$file[0]}\t$cc\t$ucc\t$methyl_region_avg{$_}{$file[1]}\t$ct\t$uct\thypo\t$S\t$E\n";
	}elsif($methyl_region_avg{$_}{$file[0]}==0 && $methyl_region_avg{$_}{$file[1]} != 0){
		print OF "$_\t$methyl_region_avg{$_}{$file[0]}\t$cc\t$ucc\t$methyl_region_avg{$_}{$file[1]}\t$ct\t$uct\thyper\t$S\t$E\n";
	}elsif($methyl_region_avg{$_}{$file[0]}!=0 && $methyl_region_avg{$_}{$file[1]} != 0){
	if($methyl_region_avg{$_}{$file[0]} / $methyl_region_avg{$_}{$file[1]} >= $diff_change){
		print OF "$_\t$methyl_region_avg{$_}{$file[0]}\t$cc\t$ucc\t$methyl_region_avg{$_}{$file[1]}\t$ct\t$uct\thypo\t$S\t$E\n";
	}elsif ($methyl_region_avg{$_}{$file[1]} / $methyl_region_avg{$_}{$file[0]} >= $diff_change){
	    print OF "$_\t$methyl_region_avg{$_}{$file[0]}\t$cc\t$ucc\t$methyl_region_avg{$_}{$file[1]}\t$ct\t$uct\thyper\t$S\t$E\n";
	}
   }	
}
}
#(%methyl_region_avg , %methyl_c , %unmethyl_c) = ();

close OF;

#calculate p value and fdr of left regions;
#system("perl fisher_test.hyp.pl");       ###p_q <- cbind(p , q  ,data\$dmr_type)
 system("perl fisher_test.hyp.pl --t $t --c_type $ty");
#p value and fdr;
my (@queue , %p_q);
open IF , "$t\_$ty\_fisher_test_results.txt" || die $!;
while(<IF>){
chomp;
#chomp (my @line = <IF>);
#my (@queue , %p_q);
#foreach (@line) {
	my @tmp = split /\t/ , $_;                ###添加$tmp[3]
	push @{$p_q{$tmp[0]}} , ($tmp[1] , $tmp[2], $tmp[3] ,$tmp[4], $tmp[5],$tmp[6],$tmp[7]);
	push @queue , $tmp[0];
#}
}
close IF;

#open OF1 , ">DMRs.txt" || die $!;
open OF2 , ">$t\_$ty\_DMRs_BH.txt" || die $!;            ####
#print OF1 "Chr\tStart\tEnd\tP_value\tFDR\tmethyl_level_1\tmethyl_level_2\tDMR_type\n";
print OF2 "Chr\tStart\tEnd\tP_value\tFDR\tmethyl_level_1\tmethyl_level_2\tDMR_type\n";
foreach (@queue) {
	my @tmp = split /_/ , $_;
	my $c = $tmp[0];
	#my $s = ($tmp[1] - 1) * $window_size + 1;
	#my $e = $tmp[1] * $window_size;                            ####
	#print OF1  "$c\t$s\t$e\t$p_q{$_}[0]\t$p_q{$_}[1]\t$p_q{$_}[2]\t$p_q{$_}[3]\t$p_q{$_}[4]\n";
	if ($p_q{$_}[0] <= $p_value && $p_q{$_}[1] <= $q_value) {    ####
	#	print OF2 "$c\t$s\t$e\t$p_q{$_}[0]\t$p_q{$_}[1]\t$p_q{$_}[2]\t$p_q{$_}[3]\t$p_q{$_}[4]\n";
	print OF2 "$c\t$p_q{$_}[5]\t$p_q{$_}[6]\t$p_q{$_}[0]\t$p_q{$_}[1]\t$p_q{$_}[2]\t$p_q{$_}[3]\t$p_q{$_}[4]\n";

	}
}

#close OF1;
close OF2;

#system 'awk '{if(FNR==1){next}{print "chr"\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$6}}' DMRs_BH.txt > DMRs_bh_.txt';
#system "grep "hyper" DMRs_bh_.txt > $treatment_$c_type_hyper_DMRs";

unlink qw($t\_$ty\_fisher_test.R   $t\_$ty\_fisher_test_results.txt);
