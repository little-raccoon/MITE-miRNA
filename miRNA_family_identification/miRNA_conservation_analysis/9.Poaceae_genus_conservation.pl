#!/usr/bin/perl
my $step6_input_file="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/miRNA_conservation/Step6_Poaceae_species.out";

open STEP6,$step6_input_file or die "fileOpenError: unable to open $step6_input_file\n";
while(<STEP6>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    if($. == 1){
        print "miRNA\tBrachypodium\tOryza\tPanicum\tSaccharum\tSetaria\tSorghum\tTriticum\tZea\ttotal\n";
        # my @genus;
        # foreach my $i(@line){
            # my @i=split "_",$i;
            # push @genus,$i[0];
        # }
        # my $line1=join "\t",@genus;
        # print "$line1\n";
    }
    else{
        my $oryza=$line[2]+$line[3]+$line[4]+$line[5];
        my $panicum=$line[6]+$line[7];
        my @line2=($line[1],$oryza,$panicum,$line[8],$line[9],$line[10],$line[11],$line[12]);
        my $count=0;
        for(my $j=0;$j<8;$j++){
            if($line2[$j]!=0){
                $line2[$j]=1;
            }
            $count+=$line2[$j];
        }
        my $line2=join "\t",@line2;
        print "$line[0]\t$line2\t$count\n";
    }
}