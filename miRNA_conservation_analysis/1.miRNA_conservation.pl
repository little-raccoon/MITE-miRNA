#!/usr/bin/perl

my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/out/All_MITE-miRNA_AllMITE_filter0.5.tab";
my $miRNA_conservation_file="/home/leili_pkuhpc/lustre1/guozhl/PmiREN/conservation/out";

open MIRNA_CONSERVATION,$miRNA_conservation_file;
my @miRNA;
while(<MIRNA_CONSERVATION>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    if($.==1){
        @miRNA=@line;
    }
    else{
        my $species=$line[0];
        my $count=0;
        foreach my $i(@line){
            if($i==1){
                my $miRNA=$miRNA[$count];
                push @$miRNA,$species;
                #print "$i\n";
            }
            $count++;
        }
    }
}
# foreach my $i(@miRNA){
    # print "$i\t@$i\n";
# }

open MITE_MIRNA,$MITE_miRNA_file;
while(<MITE_MIRNA>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $miRNA=$line[1];
    my $species=$line[0];
    if($miRNA=~/MIR(\d+)/){
        my $miRNA_fam=$1;
        if($#$miRNA_fam!=0){
            print "$species\t$miRNA\t@$miRNA_fam\n";
        }
    }
}