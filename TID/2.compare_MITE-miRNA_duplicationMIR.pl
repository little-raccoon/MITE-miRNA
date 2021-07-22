#!/usr/bin/perl

my @duplication_files=glob "/home/leili_pkuhpc/lustre1/guozhl/PMITE_dir.link/target_duplication/workdir/blastout/*filterdByOverlap";
my $MITE_miRNA_file="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/out/All_MITE-miRNA_AllMITE_filter0.5.tab";


open MITEMIRNA,$MITE_miRNA_file or die "fileOpenError: unable to open $MITE_miRNA_file\n";
while(<MITEMIRNA>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    $MITE_hash{$line[1]}=$line[4];
}
close MITEMIRNA;

foreach my $file(@duplication_files){
    open INFILE,$file or die "fileOpenError: unable to open $file\n";
    while(<INFILE>){
        chomp;
        my $line=$_;
        # print "$line\n";
        my @line=split "\t",$line;
        my $miRNA=$line[0];
        if(!exists $MITE_hash{$miRNA}){
            print "$line\n";  
        }
    }
    close INFILE;
}


# foreach my $miRNA(sort keys %duplication_hash){
    # if(exists $MITE_hash{$miRNA}){
        # if($duplication_hash{$miRNA}>$MITE_hash{$miRNA}){
            print $duplication_line_hash{$miRNA};
            print "\n";
        # }
    # }
    # else{
        # print $duplication_line_hash{$miRNA};
        # print "\n";
    # }
# }
