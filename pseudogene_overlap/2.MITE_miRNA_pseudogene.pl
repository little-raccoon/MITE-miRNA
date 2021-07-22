#!/usr/bin/perl

my @infiles=glob "/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/psedugene/intersect/*out";
my $MITE_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/MITE_miRNA.lst";

open MITE,$MITE_file or die "fileOpenError: unable to open $MITE_file\n";
while(<MITE>){
    chomp;
    my $line=$_;
    $MITE_miRNA_hash{$line}=1;
    #print "$line\n";
}

foreach my $infile(@infiles){
    #print "$infile\n";
    $infile=~/\/(\w{3}_pseudogene.*)_miRNA.intersect.out/;
    #print "$infile\n";
    my $species=$1;
    #print "$species\n";
    my $MITE_count=0;
    my $noMITE_count=0;
    open INFILE,$infile or die "fileOpenError: unable to open $infile\n";
    while(<INFILE>){
        chomp;
        my $miRNA=$_;
        #print "$miRNA\n";
        if(exists $MITE_miRNA_hash{$miRNA}){
            $MITE_count++;
        }
        else{
            $noMITE_count++;
        }
    }
    print "$species\t$MITE_count\t$noMITE_count\n";
    close INFILE;
}