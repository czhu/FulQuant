#!/usr/bin/perl
use strict;
use warnings;
#use File::Basename;

my $input = $ARGV[0];

open(my $FILE,$input) or die "Could not open file '$input' $!";
#my @fname = split("\\.",$input);
my $fname = $input =~ s/.bed//gr;

#open(my $left,">",$fname.".left.bed");
#open(my $right,">",$fname.".right.bed");
open(my $ei,">",$fname.".ei.bed");
open(my $ie,">",$fname.".ie.bed");

#print $fname[0],"\n";

while(<$FILE>){
        my @line = split("\t",$_);
        my $strand = $line[5];
        my $chromStart = $line[1];
        my $rname = $line[3];
        my @blockStarts = split(",",$line[11]);
        my @blockSizes = split(",", $line[10]);
        #print $line[10],"\t",scalar @blockStarts,"\n";
        #print $line[1],$line[3],"\n";
        #print $left $line[0],"\t",$chromStart,"\t",$chromStart+1,"\t",$rname,"\t",".","\t",$line[5],"\n";
        #print $right $line[0],"\t",$line[2],"\t",$line[2]+1,"\t",$rname,"\t",".","\t",$line[5],"\n";
        for(my $i = 0; $i< (scalar @blockStarts - 1);$i++) {
                my $eiJunc = $chromStart+$blockStarts[$i]+$blockSizes[$i];
                my $ieJunc = $chromStart + $blockStarts[$i+1];
                ## exon intro junction
                print $ei $line[0],"\t",$eiJunc-1,"\t",$eiJunc,"\t",$rname,":",sprintf("Junc%04d",$i+1),"\t",".","\t",$line[5],"\n";
                ## intron exon plut 1 for bed
                print $ie $line[0],"\t",$ieJunc,"\t",$ieJunc+1,"\t",$rname,":",sprintf("Junc%04d",$i+1),"\t",".","\t",$line[5],"\n";
        }
}

close($FILE);
#close($left);
#close($right);
close($ei);
close($ie);
