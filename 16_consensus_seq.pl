#!/usr/bin/perl
use strict;
use File::Basename qw(basename dirname);

my $ref=shift;
open(REF,$ref)||die"can't open the ref file\n";
my ($indir,$key)=@ARGV;
my @file=glob("$indir/*.hg38_multianno.txt");

my (%ref,%info,%sample)=();
$/="\>";<REF>;$/="\n";
while (<REF>) {
        chomp;
        my $chr=$_;
        $/="\>";
        my $seq=<REF>;
        chomp $seq;
        $seq=~s/\s//g;
        $/="\n";
        $ref{$chr}=$seq; 
}
close REF;

foreach my $file (@file) {
        next if($file!~/$key/ || $file=~/\-L/);
        my $name=$1 if($file=~/\w+\-(\S+)\.hg38_multianno/);
        open(FILE,$file)||die"can't open the $key $name file\n";
        while (<FILE>) {
                chomp;
                my @temp=split /\t/;
                next if($temp[5] ne "exonic" || $temp[3]!~/^[ATCG]$/ || $temp[4]!~/^[ATCG]$/ || $_=~/rs/);
                $info{"$temp[0]\t$temp[1]\t$temp[3]\>$temp[4]"}{$name}=1;
                $info{"$temp[0]\t$temp[1]\t$temp[3]\>$temp[4]"}{"Normal"}=0;
                $sample{$name}=1;
                $sample{"Normal"}=1;
        }
        close FILE;
}

foreach my $sample (sort keys %sample){
        print ">$sample\n";
        foreach my $pos (sort keys %info) {
                my @info=split /\t/,$pos;
                my ($ref,$alt)=($1,$2) if($info[2]=~/(\w)\>(\w)/);
                my $qian=substr($ref{$info[0]},$info[1]-21,20);
                my $hou=substr($ref{$info[0]},$info[1],20);
                my $wt=$qian.$ref.$hou;
                my $mt=$qian.$alt.$hou;
                if ($info{$pos}{$sample}==1) {
                        print "$mt"."NNN";
                }else{
                        print "$wt"."NNN";
                }
        }
        print "\n";
}
