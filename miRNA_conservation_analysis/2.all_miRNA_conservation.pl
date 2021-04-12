#!/usr/bin/perl

# miRNA conservation，最好计算个 conservation index
# miRNA conservation index. We estimated the gene conservation index for each of the phylogenetic profiles as the percentage of organism clusters in which the query miRNA is present in the reduced dataset.

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
            ${$species}{$miRNA_family}=1;
            $miRNA_family_hash{$miRNA_family}=1;
        }
    }
    close MIRNA;
}

print "miRNA_family\tcount_1:2\tcount_2:1\tcount_3:1\tcount_4:2\tcount_5:2\tcount_6:18\tcount_7:62\n";
foreach my $miRNA_family(sort keys %miRNA_family_hash){
    %hash=("count_1"=>0,"count_2"=>0,"count_3"=>0,"count_4"=>0,"count_5"=>0,"count_6"=>0,"count_7"=>0);
    foreach my $species(sort keys %taxon_hash){
        my $taxon=$taxon_hash{$species};
        my $count_name="count_".$taxon;
        if(exists ${$species}{$miRNA_family}){
            $hash{$count_name}++;
            ${$count_name}++;
            #print "$miRNA_family\t$count_name\n";
            #print "${$count_name}\n";
        }
        #print "$miRNA_family\t${$count_name}\n";
    }
    print "$miRNA_family";
    foreach my $i(sort keys %hash){
        #print "\t$i";
        print "\t$hash{$i}";
    }
    print "\n";
    #print "$miRNA_family\t$count_1\t$count_2\t$count_3\t$count_4\t$count_5\t$count_6\t$count_7\n";
}

