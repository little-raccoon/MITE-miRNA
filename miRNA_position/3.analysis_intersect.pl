#!/usr/bin/perl
# find miRNA in intron, exon, promoter
# Guo Zhonglong
# v1.0 2020-04-15
# usage: perl $0

##### GLOBAL VARIABLE #####

my $all_species_gff3_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_gff3_path";
my $intersect_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/intersect_out/";
my $out_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/miRNA_location_in_genome/";

##### MAIN #####

parse_all_gff($all_species_gff3_path,\%basic_info_hash);
foreach(keys %basic_info_hash){
    my $species=$_;
	open OUT,">".$out_path.$species.".miRNA_location_in_genome";
	my $exon_intersect=$intersect_path.$species."_exon.intersect";
	my $intron_intersect=$intersect_path.$species."_intron.intersect";
	my $promoter_intersect=$intersect_path.$species."_promoter.intersect";
	
	undef %exon_miRNA_hash;
	undef %intron_miRNA_hash;
	undef %promoter_miRNA_hash;
	
	parse_intersect($exon_intersect,\%exon_miRNA_hash);
	parse_intersect($intron_intersect,\%intron_miRNA_hash);
	parse_intersect($promoter_intersect,\%promoter_miRNA_hash);
	
	remove_redundant_miRNA(\%intron_miRNA_hash,\%exon_miRNA_hash);
	remove_redundant_miRNA(\%intron_miRNA_hash,\%promoter_miRNA_hash);
	remove_redundant_miRNA(\%promoter_miRNA_hash,\%exon_miRNA_hash);
    
	print OUT "exon";
	foreach(keys %exon_miRNA_hash){
	    print OUT "\t$_";
	}
	print OUT "\n";
	
	print OUT "intron";
	foreach(keys %intron_miRNA_hash){
	    print OUT "\t$_";
	}
	print OUT "\n";

	print OUT "promoter";
	foreach(keys %promoter_miRNA_hash){
	    print OUT "\t$_";
	}
	print OUT "\n";
	
	# print OUT print_miRNA("exon",\%exon_miRNA_hash);
	# print OUT print_miRNA("intron",\%intron_miRNA_hash);
	# print OUT print_miRNA("promoter",\%promoter_miRNA_hash);
	close OUT;
}

##### FUNCTION #####

sub print_hash{
    my($hash)=@_;
	foreach(keys %$hash){
	    print "$_\t$$hash{$_}\n";
	}
}

sub print_miRNA{
    my($feature,$hash)=@_;
	foreach(keys %$hash){
	    $out.="\t$_";
	}
	print "$feature$out\n";
}

sub remove_redundant_miRNA{
    my($hash1,$hash2)=@_;
	foreach(keys %$hash1){
	    my $key1=$_;
		if(exists $$hash2{$key1}){
			#print "$key1\n";
		    if($$hash1{$key1}>=$$hash2{$key1}){
			    delete $$hash2{$key1};
			}
			else{
			    delete $$hash1{$key1};
			}
		}
	}
}

sub parse_intersect{
    my($infile,$hash)=@_;
	open INFILE,$infile or die "errorFileOpen: unable to open file: $infile\n";
	while(<INFILE>){
	    chomp;
		my $line=$_;
		my @line=split "\t",$line;
		my($miRNA,$len)=($line[17],$line[-1]);
		#print "$miRNA\t$len\n";
		if(!exists $$hash{$miRNA}){
		    $$hash{$miRNA}=$len;
			
		}
		else{
		    if($len>$$hash{$miRNA}){
			    $$hash{$miRNA}=$len;
			}
			#print "$miRNA\t$len\n";
		}
	}
	# foreach(keys %$hash){
	    # print "$_\t$$hash{$_}\n";
	# }
}

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