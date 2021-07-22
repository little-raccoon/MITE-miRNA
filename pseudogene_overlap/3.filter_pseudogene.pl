#!/usr/bin/perl

my $intersect_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/psedugene/intersect/";
my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/MITE_miRNA.lst";


my @intersect_files=glob($intersect_path."*upstream2k_miRNA.intersect");

open MITE_MIRNA,$MITE_miRNA_file or die "fileOpenError: unable to open $MITE_miRNA_file\n";
while(<MITE_MIRNA>){
    chomp;
    my $line=$_;
    $MITE_miRNA_hash{$line}=1;
}


foreach my $file(@intersect_files){
    open INFILE,$file or die "fileOpenError: unable to open $file\n";
    $file=~/(...)_pseudogene/;
    #print "$file\n";
    my $species=$1;
    #print "$species\n";
    my $count;
    my %uniq_miRNA_hash=();
    my %uniq_MITE_miRNA_hash=();
    #print "$file\n";
    while(<INFILE>){
        chomp;
        my $line=$_;
        my @line=split "\t",$line;
        my($num1,$num2,$num3,$num4,$miRNA)=($line[1],$line[2],$line[7],$line[8],$line[5]);
        my $miRNA_len=$num2-$num1+1;
        my @sort_position=sort ($num1,$num2,$num3,$num4);
        my $match_len=$sort_position[2]-$sort_position[1]+1;
        #print "@sort_position\n";
        if($match_len/$miRNA_len>=0.5){
            $uniq_miRNA_hash{$miRNA}=1;
            if(exists $MITE_miRNA_hash{$miRNA}){
                $uniq_MITE_miRNA_hash{$miRNA}=1;
            }
            #$count++;
        }
    }
    my $pse_count=keys %uniq_miRNA_hash;
    my $MITE_pse_count=keys %uniq_MITE_miRNA_hash;
    print "$species\t$pse_count\t$MITE_pse_count\n";
}


