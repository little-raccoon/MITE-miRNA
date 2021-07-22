#!/usr/bin/perl -w
# split gene_gff, exon_gff, intron_gff, promoter_gff from GFF3 filefield
# Guo Zhonglong
# v1.0 2020-04-14
# successful running --2020-4-15

##### GLOBAL VARIABLE #####

my $all_species_gff3_path ="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/22_species_gff3_path";
my $out_dir="/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/intron_MITE/data/";

##### MAIN #####

parse_all_gff($all_species_gff3_path,\%gff_path_hash);
foreach(keys %gff_path_hash){
    my $species=$_;
	my $gff_path=$gff_path_hash{$species};
    split_gff("mRNA",$species,$gff_path,$out_dir);
	split_gff("exon",$species,$gff_path,$out_dir);
	promoter_gff(2000,$species,$gff_path,$out_dir);
	my $mRNA_gff=$out_dir."mRNA_gff/".$species."_mRNA.gff";
	my $exon_gff=$out_dir."exon_gff/".$species."_exon.gff";
	my $intron_gff=$out_dir."intron_gff/".$species."_intron.gff";
	system("bedtools subtract -a $mRNA_gff -b $exon_gff > $intron_gff");
}

##### FUNCTION #####

sub promoter_gff{
    my($promoter_len,$species,$gff_path,$out)=@_;
	open INFILE,$gff_path or die "errorFileOpen: unable to open file $$infile\n";
	open OUT,">".$out."promoter_gff/".$species."_promoter.gff";
	while(<INFILE>){
	    chomp;
		my $line=$_;
		if($line!~/^#/){
		    my @line=split "\t",$line;
			my $gff_feature=$line[2];
			my $strand=$line[6];
			my($start,$end)=($line[3],$line[4]);
			if($gff_feature eq "mRNA"){
			    if($strand eq "+"){
				    $line[3]=$start-$promoter_len;
					$line[4]=$start-1;
				}
				elsif($strand eq "-"){
				    $line[3]=$end+1;
					$line[4]=$end+$promoter_len;
				}
				else{
				    print "error: strand error ".$species."\t".$gff_path."\n";
					}
				$line[2]="promoter";
				my $output=join "\t",@line;
				print OUT "$output\n";
			}
		}
	}
}

sub split_gff{
    my($feature,$species,$gff_path,$out)=@_;
	open INFILE,$gff_path or die "errorFileOpen: unable to open file $$infile\n";
	open OUT,">".$out.$feature."_gff/".$species."_".$feature.".gff";
	while(<INFILE>){
	    chomp;
		my $line=$_;
		if($line!~/^#/){
		    my @line=split "\t",$_;
			my $gff_feature=$line[2];
			if($feature eq $gff_feature){
			    print OUT "$line\n";
			}
		}
	}
}

sub parse_all_gff{
    my($infile,$hash)=@_;
	open INFILE,$infile or die "errorFileOpen: unable to open file $$infile\n";
	while(<INFILE>){
	    chomp;
		my @line=split "\t",$_;
		my($species,$gff_path)=@line;
		$$hash{$species}=$gff_path;
	}
}