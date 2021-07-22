#!/usr/bin/perl

my $species_file="/home/leili_pkuhpc/lustre1/guozhl/PMITE_dir.link/Oryza_sRNA/species.lst";
my @statis_files=glob "/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/Oryza_sRNA/miRNA/statis/*statis";

open SPECIES,$species_file or die "fileOpenError: unable to open $species_file\n";
while(<SPECIES>){
    chomp;
    my @line=split "\t",$_;
    $species_hash{$line[0]}=$line[1];
    $reads_hash{$line[0]}=$line[2];
}
close SPECIES;


foreach my $file(@statis_files){
    my $lib;
    my $species;
    
    if($file=~/(SRR\d*)/){
        $lib=$1;
    }
    else{die "cant find lib id in file $file\n";}
    if(exists $species_hash{$lib}){
        $species=$species_hash{$lib};
    }
    else{die "cant find species name of lib: $lib\n";}
    
    open OUT,">>".$species.".rpm";
    
    
    open FILE,$file or die "fileOpenError: unable to open $file\n";
    while(<FILE>){
        chomp;
        my @line=split "\t",$_;
        my($miRNA,$reads)=@line;
        my $rpm;
        if(exists $reads_hash{$lib}){
            my $total_reads=$reads_hash{$lib};
            $rpm=$reads/$total_reads*1000000;
        }
        else{die "cant find total reads of lib: $lib\n";}
        print OUT "$miRNA\t$rpm\n";
    }
    close FILE;
    close OUT;
}