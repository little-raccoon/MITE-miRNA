#!/usr/bin/perl

my $species_taxon_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/miRNA_conservation/species.txt";
my $miRNA_file_path="/home/leili_pkuhpc/lustre1/guozhl/PmiREN/blastDB/pre/";
my @miRNA_files=glob $miRNA_file_path."*fa";

open SPECIES_TAXON,$species_taxon_file or die "fileOpenError: unable to open $species_taxon_file\n";
while(<SPECIES_TAXON>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my($short_name,$species,$family,$taxon_num)=($line[0],$line[1],$line[-2],$line[-1]);
    $taxon_hash{$species}=$taxon_num;
    #print "$species\t$taxon_num\n";
    
}
close SPECIES_TAXON;

foreach my $miRNA_file(@miRNA_files){
    $miRNA_file=~/\/(\w+_\w+)\.hairpin\.fa/;
    my $species=$1;
    #print "$species\n";
    open MIRNA,$miRNA_file or die "fileOpenError: unable to open $miRNA_file\n";
    while(<MIRNA>){
        chomp;
        if($_=~/^>/){
            my $miRNA=$_;
            $miRNA=~/(MIRN?\d+)/;
            my $miRNA_family=$1;
            #print "$miRNA_family\n";
            ${$miRNA_family}{$species}=1;
            $miRNA_family_hash{$miRNA_family}=1;
        }
    }
    close MIRNA;
}

foreach my $miRNA_family(sort keys %miRNA_family_hash){
    #print "$miRNA_family\n";
    my @species_group=keys %$miRNA_family;
    #print "@species_group\n";
    my $species_count=@species_group;
    print "$miRNA_family\t$species_count\n"
}