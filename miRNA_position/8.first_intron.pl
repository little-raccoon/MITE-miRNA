#!/usr/bin/perl

my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/MITE_miRNA.lst";
my $intronMiRNA_intersect_path="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/intron_MITE/intersect_out/*_intron.intersect";

open MITE_MIRNA,$MITE_miRNA_file or die "fileOpenError: unable to open $MITE_miRNA_file\n";
while(<MITE_MIRNA>){
    chomp;
    my $line=$_;
    $MITE_miRNA_hash{$line}=1;
}
close MITE_MIRNA;


my @intronMiRNA_intersect_files=glob $intronMiRNA_intersect_path;
foreach my $infile(@intronMiRNA_intersect_files){
    open INFILE,$infile or die "fileOpenError: unable to open $infile\n";
    while(<INFILE>){
        chomp;
        my $line=$_;
        my @line=split "\t",$line;
        my($strand,$gene,$miRNA,$start,$end)=($line[6],$line[8],$line[17],$line[3],$line[4]);
        if($gene=~/ID=([^;]*);/){
            $gene=$1;
            if(exists $MITE_miRNA_hash{$miRNA}){
                print "$strand\t$gene\t$start\t$end\t$miRNA\n";
            }
        }
        else{
            print "error1\n";
        }
    }
    close INFILE;
}