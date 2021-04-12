#!/usr/bin/perl

my $step2_out_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/miRNA_conservation/Step2_miRNA_conservation.out";

open STEP2_OUT,$step2_out_file or die "fileOpenError: unable to open $step2_out_file\n";
while(<STEP2_OUT>){
    chomp;
    if($.>=2){
        my $line=$_;
        my @line=split "\t",$line;
        my($miRNA,$count1,$count2,$count3,$count4,$count5,$count6,$count7)=@line;
        my @out;
        #$out[0]=$miRNA;
        foreach my $i(@line[1..5]){
            if($i!=0){
                push @out,"T";
            }
            else{
                push @out,"F";
            }
        }
        if($count6>=5.4){
            push @out,"T";
        }
        else{
            push @out,"F";
        }
        if($count7>=18.6){
            push @out,"T";
        }
        else{
            push @out,"F";
        }
    my $out=join "",@out;
    push @{$hash{$out}},$miRNA;
    #print "$miRNA\t$out\n";
    }
}

print "type\tcount\tmiRNA\n";
foreach my $i(sort keys %hash){
    my $out=join ",",@{$hash{$i}};
    my $count=@{$hash{$i}};
    print "$i\t$count\t$out\n";
}
