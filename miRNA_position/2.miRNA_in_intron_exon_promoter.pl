#!/usr/bin/perl
# find miRNA in intron, exon, promoter
# Guo Zhonglong
# v1.0 2020-04-15
#

##### GLOBAL VARIABLE #####

my $data_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/";
my $all_species_gff3_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_gff3_path";
my $miRNA_gff_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/miRNA_new_gff/pre_miRNA/sorted_gff/";
my $intersect_out="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/intersect_out/";

##### MAIN #####

parse_all_gff($all_species_gff3_path,\%basic_info_hash);
foreach(keys %basic_info_hash){
    chomp;
	my $species=$_;
	my $basic_info_file=$basic_info_hash{$species};
	my $miRNA_gff_file=$miRNA_gff_path.$species."_miRNA.gff";
	#basic2gff($basic_info_file,$miRNA_gff_file); 
	# change miRNA-basic-info file to gff. Beacuse the genome version of Lotus_japonicus has been changed, Lja_pre-miRNAs were mapped to version3.0 genome using bowtie; then I converted sam file to gff manually.
	my $exon_gff=$data_path."exon_gff/".$species."_exon.gff";
	my $intron_gff=$data_path."intron_gff/".$species."_intron.gff";
	my $promoter_gff=$data_path."promoter_gff/".$species."_promoter.gff";
	my $miRNA_gff=$miRNA_gff_path.$species."_miRNA.gff";
	#print "$miRNA_gff\n";
	my $exon_out=$intersect_out.$species."_exon.intersect";
	my $intron_out=$intersect_out.$species."_intron.intersect";
	my $promoter_out=$intersect_out.$species."_promoter.intersect";
	system("~/lustre1/software/bedtools_v2.29.2/bedtools2/bin/bedtools intersect -a $exon_gff -b $miRNA_gff -wo > $exon_out");
	system("~/lustre1/software/bedtools_v2.29.2/bedtools2/bin/bedtools intersect -a $intron_gff -b $miRNA_gff -wo > $intron_out");
	system("~/lustre1/software/bedtools_v2.29.2/bedtools2/bin/bedtools intersect -a $promoter_gff -b $miRNA_gff -wo > $promoter_out");
	
}

##### FUNCTION #####

# sub basic2gff{
    # my($basic_file,$gff_file)=@_;
	# open INFILE,$basic_file or die "unable to open file $basic_file\n";
	# open OUTFILE,">".$gff_file;
	# while(<INFILE>){
	    # chomp;
			# if($.!=1){
			# my $line=$_;
			# my @line=split "\t",$line;
			# my($chr,$start,$end,$strand,$name)=($line[4],$line[5],$line[6],$line[7],$line[0]);
			# print OUTFILE "$chr\tPmiREN\tPre-miRNA\t$start\t$end\t.\t$strand\t.\t$name\n";
		# }
	# }
	# close INFILE;
# }

sub parse_all_gff{
    my($infile,$hash)=@_;
	open INFILE,$infile or die "errorFileOpen: unable to open file $$infile\n";
	while(<INFILE>){
	    chomp;
		my @line=split "\t",$_;
		my($species,$gff_path,$basic_info)=@line;
		$$hash{$species}=$basic_info;
	}
	close INFILE;
}
