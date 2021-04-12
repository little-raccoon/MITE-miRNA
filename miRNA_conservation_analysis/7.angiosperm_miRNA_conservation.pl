#!/usr/bin/perl

# 科内miRNA保守性，输入科的名字

my $query_family="Poaceae";

my $species_taxon_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/miRNA_conservation/species.txt";
my $miRNA_file_path="/home/leili_pkuhpc/lustre1/guozhl/PmiREN/blastDB/pre/";
my @miRNA_files=glob $miRNA_file_path."*fa";

open SPECIES_TAXON,$species_taxon_file or die "fileOpenError: unable to open $species_taxon_file\n";
while(<SPECIES_TAXON>){
    chomp;
    if($.>6){
        my $line=$_;
        my @line=split "\t",$line;
        my($short_name,$species,$family,$taxon_num)=($line[0],$line[1],$line[-2],$line[-1]);
        $taxon_hash{$species}="Angiosperm";
    }
}
close SPECIES_TAXON;

foreach my $miRNA_file(@miRNA_files){
    $miRNA_file=~/\/(\w+_\w+)\.hairpin\.fa/;
    my $species=$1;
    #print "$species\n";
    if(exists $taxon_hash{$species}){
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
}

my $species_num=0;
foreach my $species(sort keys %taxon_hash){
    $species_num++;
}

#1,2,3,4分别代表不同保守性等级，1表示只有一个物种中存在，2表示在<50%的物种中存在，3表示在>=50%并且<80%的物种中存在，4表示在>=80%的物种中存在。
foreach my $miRNA_family(sort keys %miRNA_family_hash){
    print "$miRNA_family";
    my $count=0;
    foreach my $species(sort keys %taxon_hash){
        if(exists ${$species}{$miRNA_family}){
            $count++;
        }
    }
    if($count=="1"){
        print "\t1";
    }
    elsif($count/$species_num>=0.8){
        print "\t4";
    }
    elsif($count/$species_num>=0.5 && $count/$species_num<0.8){
        print "\t3";
    }
    elsif($count/$species_num<0.5 && $count!=1){
        print "\t2";
    }
    print "\n";
}

