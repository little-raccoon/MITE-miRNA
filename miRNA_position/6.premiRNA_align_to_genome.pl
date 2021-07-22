#!/usr/bin/perl

my $path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/miRNA_new_gff/pre_miRNA/";
my $bowtie_index_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/miRNA_new_gff/bowtie_index_path.lst";

##### MAIN #####

parse_bowtie_index_path($bowtie_index_path,\%bowtie_index_hash);
foreach my $species(sort keys %bowtie_index_hash){
    #print "$species\n";
    my $pre_miRNA_file=$path.$species.".hairpin.fa";
    my $index=$bowtie_index_hash{$species};
    my $sam_out=$path.$species.".pre.sam";
    #print "$pre_miRNA\n";
    system("bowtie -p 10 $index -f $pre_miRNA_file -S $sam_out");
}

##### FUNCTION #####

sub parse_bowtie_index_path{
    my($infile,$hash)=@_;
    open INFILE,$infile or die "unable to open $infile\n";
    while(<INFILE>){
        chomp;
        my $line=$_;
        my @line=split "\t",$line;
        my ($species,$index)=@line;
        $$hash{$species}=$index;
        #print "$species\t$index\n";
    }
    close INFILE;
}

