#/usr/bin/perl
# extract 5' upstream 2k bp region (bed format)
# v1.0
# 2020-07-04

my $usage="
    perl $0\n
";

@pseudoFile=glob "/home/leili_pkuhpc/lustre1/guozhl/MITEs/PMITEdb/psedugene/sort_bed/*pseudogene.sort.bed";

foreach my $infile(@pseudoFile){
    my $outfile=$infile;
    $outfile=~s/pseudogene.sort.bed/pseudogene_upstream2k.sort.bed/g;
    extract_upstream_2k($infile,$outfile);
}


#####################
sub extract_upstream_2k{
    my($infile,$outfile)=@_;
    open INFILE,$infile or die "fileOpenError: unable to open $infile\n";
    open OUT,">".$outfile;
    while(<INFILE>){
        chomp;
        my $line=$_;
        my @line=split "\t",$line;
        my ($chr,$start,$end,$strand)=@line;
        my $up_start;
        my $up_end;
        if($strand eq "+"){
            $up_start=$start-2000;
            if($up_start < 0){
                $up_start=0;
            }
            $up_end=$start-1;
        }
        elsif($strand eq "-"){
            $up_start=$end+1;
            $up_end=$end+2000;
        }
        else{
            print "error: row $. colume 4 (strand) not + or -\n";
        }
        print OUT "$chr\t$up_start\t$up_end\t$strand\n";
    }
    close INFILE;
    close OUT;
}