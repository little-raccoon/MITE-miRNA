#!/usr/bin/perl

# 生成科之间的新miRNA数量的二维矩阵

# 排除的科 Arecaceae  Orchidaceae

# 按照演化关系排列
my $query_family="Malvaceae";

my $species_taxon_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/miRNA_conservation/species.txt";
my $miRNA_file_path="/home/leili_pkuhpc/lustre1/guozhl/PmiREN/blastDB/pre/";
my @miRNA_files=glob $miRNA_file_path."*fa";

open SPECIES_TAXON,$species_taxon_file or die "fileOpenError: unable to open $species_taxon_file\n";
while(<SPECIES_TAXON>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my($short_name,$species,$family,$taxon_num)=($line[0],$line[1],$line[-2],$line[-1]);
    $taxon_hash{$species}=$family;
}
close SPECIES_TAXON;

foreach my $miRNA_file(@miRNA_files){
    $miRNA_file=~/\/(\w+_\w+)\.hairpin\.fa/;
    my $species=$1;
    #print "$species\n";
    if(exists $taxon_hash{$species}){
        my $family=$taxon_hash{$species};
        if($family eq "Arecaceae" || $family eq "Orchidaceae"){
        }
        else{
            open MIRNA,$miRNA_file or die "fileOpenError: unable to open $miRNA_file\n";
            while(<MIRNA>){
                chomp;
                if($_=~/^>/){
                    my $miRNA=$_;
                    $miRNA=~/(MIRN?\d+)/;
                    my $miRNA_family=$1;
                    #print "$miRNA_family\n";
                    ${$family}{$miRNA_family}=1;
                    $total_family{$family}=1;
                    $miRNA_family_hash{$miRNA_family}=1;
                }
            }
            close MIRNA;
        }
    }
    else{
        print STDERR "$species not in hash taxon_hash\n";
    }
}

print "family";

my $lst_file="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/miRNA_conservation/lst";

my @species_lst;
open LST,$lst_file or die "fileOpenError: $lst_file\n";
while(<LST>){
    chomp;
    push @species_lst,$_;
}
close LST;

foreach my $k(@species_lst){
    print "\t$k";
}
print "$\n";

# foreach my $k(sort keys %total_family){
    # print "\t$k";
# }
# print "\n";


foreach my $col_family(@species_lst){
    print "$col_family";
    foreach my $row_family(@species_lst){
        my $total_count=0;
        my $novel_count=0;
        my $novel_rate=0;
        foreach my $i(sort keys %$col_family){
            $total_count++;
            if(!exists $$row_family{$i}){
                $novel_count++;
            }
        }
        my $novel_rate;
        if($total_count!=0){
            $novel_rate=$novel_count/$total_count;
            print "\t$novel_rate";
        }
        else{
            print "error: total_count 0";
        }
    }
   print "\n";
}