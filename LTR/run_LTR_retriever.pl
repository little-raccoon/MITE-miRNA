#!/usr/bin/perl

my $input_files="/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/LTR/input_files.lst";
open INPUT,$input_files or die "fileOpenError: unable to open $input_files\n";
while(<INPUT>){
    chomp;
    my @line=split "\t",$_;
    my($scn_file,$genome_file)=@line;
    system("pkurun-cnlong 1 1 LTR_retriever -genome $genome_file -infinder $scn_file -threads 1");
}
