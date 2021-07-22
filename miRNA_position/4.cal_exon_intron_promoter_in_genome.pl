#!/usr/bin/perl
# caculate length of bp in exon, intron, promoter and intergenic.
# Guo Zhonglong
# v1.0 20200416

##### GLOBAL VARIABLE #####

my $all_genome_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_genome_path.txt";
my $data_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/";
my $out="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_exon_intron_promoter_intergenic_len.txt";

##### MAIN #####

open OUT,">".$out;
print OUT "species\texon_len\tintron_len\tpromoter_len\tintergenic_len\n";
parse_genome_path($all_genome_path,\%genome_path_hash);
foreach my $species(keys %genome_path_hash){
    #print "$species\n";
    my $genome_path=$genome_path_hash{$species};
    my $exon_bed=$data_path."exon_gff/sort/".$species."_exon.bed";
	my $intron_bed=$data_path."intron_gff/sort/".$species."_intron.bed";
	my $promoter_bed=$data_path."promoter_gff/sort/".$species."_promoter.bed";
	my $exon_len;
	my $intron_len;
	my $promoter_len;
	cal_len($exon_bed,\$exon_len);
	cal_len($intron_bed,\$intron_len);
	cal_len($promoter_bed,\$promoter_len);
	cal_genome_len($genome_path,\$genome_len);
	my $intergenic_len=$genome_len-$exon_len-$intron_len-$promoter_len;
	#print "$exon_len\n";
	print OUT "$species\t$exon_len\t$intron_len\t$promoter_len\t$intergenic_len\n";
}
close OUT;

#### FUNCTION #####

sub cal_genome_len{
    my($file,$total_len)=@_;
	open INFILE,$file or die "unable to open $file\n";
	while(<INFILE>){
	    chomp;
		my $line=$_;
		if($line!~/^>/){
		    my $len=length $line;
			$$total_len+=$len;
		}
	}
}

sub cal_len{
    my($file,$total_len)=@_;
	open INFILE,$file or die "unable to open $file\n";
	while(<INFILE>){
	    chomp;
		my @line=split "\t",$_;
		my $len=$line[2]-$line[1];
		$$total_len+=$len;
	}
	#print "$total_len\n";
}

sub parse_genome_path{
    my($infile,$hash)=@_;
	open INFILE,$infile or die "errorFileOpen: unable to open file $$infile\n";
	while(<INFILE>){
	    chomp;
		my @line=split "\t",$_;
		my($species,$genome_path)=@line;
		$$hash{$species}=$genome_path;
	}
	close INFILE;
}