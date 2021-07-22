#!/usr/bin/perl
my $file=shift;

open INFILE,$file or die "fileOpenError: unable to open $file\n";
while(<INFILE>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my($type,$description,$GO)=@line[0..2];
    if($type=~/^Biological/){
        if($description=~/metabolic/ || $description=~/biosynthetic/ || $description=~/catabolic/){
            print "Metabolic\t$line\n";
        }
        elsif($description=~/transport/){
            print "Transport\t$line\n";
        }
        else{print "\t$line\n";}
    }
    else{print "\t$line\n";}
}
close INFILE;

#my @metabolic_num=keys %metabolic_hash;
#my $metabolic_num=@metabolic_num;
#print "metabolic_num: $metabolic_num\n";