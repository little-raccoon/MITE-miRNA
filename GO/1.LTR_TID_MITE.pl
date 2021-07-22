#!/usr/bin/perl


my $MITE_GO_file="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/target_gene/analysis/MITE-miRNAtarget_allgene.out";
my $LTR_GO_file="LTR_GO_Enrichment.out";
my $TID_GO_file="TID_GO_Enrichment.out";

open MITE,$MITE_GO_file or die "fileOpenError: unable to open $MITE_GO_file\n";
open LTR,$LTR_GO_file or die "fileOpenError: unable to open $LTR_GO_file\n";
open TID,$TID_GO_file or die "fileOpenError: unable to open $TID_GO_file\n";

while(<MITE>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $GO=$line[2];
    $MITE_hash{$GO}=$line;
    $all_hash{$GO}=1;
}

while(<LTR>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $GO=$line[2];
    $LTR_hash{$GO}=$line;
    $all_hash{$GO}=1;
}

while(<TID>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $GO=$line[2];
    $TID_hash{$GO}=$line;
    $all_hash{$GO}=1;
}
 
close MITE;
close LTR;
close TID;

$out1="./MITE_only.out";
$out2="./LTR_only.out";
$out3="./TID_only.out";
$out4="./MITE_LTR.out";
$out5="./MITE_TID.out";
$out6="./LTR_TID.out";
$out7="./MITE_LTR_TID.out";

open OUT1,">".$out1;
open OUT2,">".$out2;
open OUT3,">".$out3;
open OUT4,">".$out4;
open OUT5,">".$out5;
open OUT6,">".$out6;
open OUT7,">".$out7;

foreach my $GO(sort keys %all_hash){
    my $num=0;  #MITE, 1; LTR, 2; TID, 4
    if(exists $MITE_hash{$GO}){
        $num+=1;
    }
    if(exists $LTR_hash{$GO}){
        $num+=2;
    }
    if(exists $TID_hash{$GO}){
        $num+=4;
    }
    if($num == 1){
        print OUT1 "$MITE_hash{$GO}\n";
    }
    elsif($num == 2){
        print OUT2 "$LTR_hash{$GO}\n";
    }
    elsif($num == 4){
        print OUT3 "$TID_hash{$GO}\n";
    }
    elsif($num == 3){
        print OUT4 "$MITE_hash{$GO}\n$LTR_hash{$GO}\n";
    }
    elsif($num == 5){
        print OUT5 "$MITE_hash{$GO}\n$TID_hash{$GO}\n";
    }
    elsif($num == 6){
        print OUT6 "LTR_hash{$GO}\n$TID_hash{$GO}\n";
    }
    elsif($num == 7){
        print OUT7 "$MITE_hash{$GO}\n$LTR_hash{$GO}\n$TID_hash{$GO}\n";
    }
    else{
        print "error: $GO not in\n";
    }
}

close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
close OUT7;