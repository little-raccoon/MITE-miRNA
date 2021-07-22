#!/usr/bin/perl

my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/MITE_miRNA.lst";
my $miRNA_location_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/miRNA_location_in_genome/*miRNA_location_in_genome";


open MITE,$MITE_miRNA_file or die "fileOpenError: unable to open $MITE_miRNA_file\n";
while(<MITE>){
    chomp;
    $MITE_miRNA_hash{$_}=1;
}
close MITE;


my @files=glob($miRNA_location_path);
foreach my $file(@files){
    $file=~/genome\/(\S+).miRNA_location_in_genome/;
    my $species=$1;
    if($file=~/miRNA_location_in_genome/){
        open FILE,$file or die "fileOpenError: unable to open $file\n";
        my $MITE_miRNA_intron=0;
        my $MITE_miRNA_exon=0;
        my $MITE_miRNA_promoter=0;
        my $miRNA_intron=0;
        my $miRNA_exon=0;
        my $miRNA_promoter=0;
        while(<FILE>){
            chomp;
            my $exon=$_;
            my $intron=<FILE>;
            my $promoter=<FILE>;
            chomp $intron;
            chomp $promoter;
            my @exon=split "\t",$exon;
            shift @exon;
            my @intron=split "\t",$intron;
            shift @intron;
            my @promoter=split "\t",$promoter;
            shift @promoter;
            $miRNA_intron=@intron;
            $miRNA_exon=@exon;
            $miRNA_promoter=@promoter;
            foreach my $i(@exon){
                if(exists $MITE_miRNA_hash{$i}){
                    $MITE_miRNA_exon++;
                }
            }
            foreach my $i(@intron){
                if(exists $MITE_miRNA_hash{$i}){
                    $MITE_miRNA_intron++;
                }
            }
            foreach my $i(@promoter){
                if(exists $MITE_miRNA_hash{$i}){
                    $MITE_miRNA_promoter++;
                }
            }
        }
        close FILE;
        print "$species\t$miRNA_exon\t$miRNA_intron\t$miRNA_promoter\t$MITE_miRNA_exon\t$MITE_miRNA_intron\t$MITE_miRNA_promoter\n";
    }
}