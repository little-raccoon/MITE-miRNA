#!/usr/bin/perl

my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/MITE_miRNA.lst";
my $miRNA_conservation_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/miRNA_conservation/Step3_miRNA_conservation_filter0.3.out";


open MIRNA_CONSERVATION,$miRNA_conservation_file;
while(<MIRNA_CONSERVATION>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my($type,$count,$family)=@line;
    my @family=split ",",$family;
    foreach my $i(@family){
        $hash{$i}=$type;
    }
}

open MITE_MIRNA,$MITE_miRNA_file;
while(<MITE_MIRNA>){
    chomp;
    my $line=$_;
    $line=~/(MIRN?\d+)/;
    my $family=$1;
    if(exists $hash{$family}){
        my $type=$hash{$family};
        print "$family\t$line\t$type\n";
    }
    else{
        #print "$line\t$family\tterror1\n";
    }
}


