#!/usr/bin/perl
# normalized miRNA numbers in exon, intron, promoter and intergenic.
# Guo Zhonglong
# v1.0 20200416

##### INPUT #####

my $option=shift;  #opetion is -miRNA or -MITE_miRNA

##### GLOBAL VARIABLE #####

my $miRNA_num_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/miRNA_location_in_genome/";
my $miRNA_gff_path="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/miRNA_gff/";
my $MITE_miRNA_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/out/All_MITE-miRNA_AllMITE_filter0.5.tab";
my $length_file="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_exon_intron_promoter_intergenic_len.txt";
my $out="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/normlized_twice_number.".$option;

##### MAIN #####

parse_feature_size($length_file,\%feature_size_hash);
parse_MITE_miRNA($MITE_miRNA_file,\%MITE_miRNA_hash);
# foreach(keys %MITE_miRNA_hash){
    # print "$_\t$MITE_miRNA_hash{$_}\n";
# }
open OUT1,">".$out.".norm1";
open OUT2,">".$out.".norm2";
print OUT1 "species\texon\tintron\tpromoter\tintergenic\n";
print OUT2 "species\texon\tintron\tpromoter\tintergenic\n";
foreach my $species(keys %feature_size_hash){
    my @ratio=split "\t",$feature_size_hash{$species};
	my($exon_ratio,$intron_ratio,$promoter_ratio,$intergenic_ratio)=@ratio;
	
    my $miRNA_file=$miRNA_num_path.$species.".miRNA_location_in_genome";
	if($option eq "-miRNA"){
	    miRNA_num($miRNA_file,\%miRNA_num_hash);
	}
	elsif($option eq "-MITE_miRNA"){
	    MITE_miRNA_num($miRNA_file,\%miRNA_num_hash);
	}
	else{
	    print "error in option: -miRNA or -MITE_miRNA";
	}
	my $exon_miRNA_num=$miRNA_num_hash{"exon"};
	my $intron_miRNA_num=$miRNA_num_hash{"intron"};
	my $promoter_miRNA_num=$miRNA_num_hash{"promoter"};
	
	my $miRNA_gff=$miRNA_gff_path.$species."_miRNA.gff";
	my $total_miRNA_num=all_miRNA_num($miRNA_gff);
	my $intergenic_miRNA_num=$total_miRNA_num-$exon_miRNA_num-$intron_miRNA_num-$promoter_miRNA_num;
	if($exon_ratio>0 && $intron_ratio>0 && $promoter_ratio>0 && $intergenic_ratio>0){
		my $norm1_exon=$exon_miRNA_num/$exon_ratio;
		my $norm1_intron=$intron_miRNA_num/$intron_ratio;
		my $norm1_promoter=$promoter_miRNA_num/$promoter_ratio;
		my $norm1_intergenic=$intergenic_miRNA_num/$intergenic_ratio;
		my $norm1_total=$norm1_exon+$norm1_intron+$norm1_promoter+$norm1_intergenic;
                print OUT1 "$species\t$norm1_exon\t$norm1_inron\t$norm1_promoter\t$norm1_intergenic\t$norm1_total\n";
		if($norm1_total>0){
			my $norm2_exon=$norm1_exon/$norm1_total;
			my $norm2_intron=$norm1_intron/$norm1_total;
			my $norm2_promoter=$norm1_promoter/$norm1_total;
			my $norm2_intergenic=$norm1_intergenic/$norm1_total;
			print OUT2 "$species\t$norm2_exon\t$norm2_intron\t$norm2_promoter\t$norm2_intergenic\n";
		}
		else{
		    print "error2: $species\n";
		}
	}
	else{
	    print "error1: $species\n";
	}
}
close OUT1;
close OUT2;

##### FUNCTION #####
sub parse_MITE_miRNA{
    my($file,$hash)=@_;
	open INFILE,$file or die "unable to open file $file\n";
	while(<INFILE>){
	   chomp;
	   my $line=$_;
	   my @line=split "\t",$line;
	   my($species,$miRNA)=@line;
	   $$hash{$miRNA}="T";
	}
}

sub all_miRNA_num{
    my($file)=@_;
	open INFILE,$file or die "unable to open file $file\n";
	my $miRNA_num;
	while(<INFILE>){
	    chomp;
		$miRNA_num++;
	}
	return $miRNA_num;
}

sub MITE_miRNA_num{
    my($file,$hash)=@_;
	open INFILE,$file or die "unable to open $file\n";
	while(<INFILE>){
	    my $line=$_;
		my @line=split "\t",$line;
		my $count=0;
		my $feature=shift @line;
		foreach my $miRNA(@line){
		    if(exists $MITE_miRNA_hash{$miRNA}){
			    $count++;
			}
		}
		$$hash{$feature}=$count;
	}
	close INFILE;
}

sub miRNA_num{
    my($file,$hash)=@_;
	open INFILE,$file or die "unable to open $file\n";
	while(<INFILE>){
	    my $line=$_;
		my @line=split "\t",$line;
		my $feature=shift @line;
		$$hash{$feature}=@line;
	}
	close INFILE;
}

sub parse_feature_size{
    my($file,$hash)=@_;
	open INFILE,$file or die "unable to open $file\n";
	while(<INFILE>){
	    chomp;
		my $line=$_;
		my @line=split "\t",$line;
		if($.!=1){
		    my($species,$exon,$intron,$promoter,$intergenic)=@line;
			my $total=$exon+$intron+$promoter+$intergenic;
			my $exon_ratio=$exon/$total;
			my $intron_ratio=$intron/$total;
			my $promoter_ratio=$promoter/$total;
			my $intergenic_ratio=$intergenic/$total;
            $$hash{$species}=$exon_ratio."\t".$intron_ratio."\t".$promoter_ratio."\t".$intergenic_ratio;
			#print "$species\t$$hash{$species}\n";
		}
	}
	close INFILE;
}
