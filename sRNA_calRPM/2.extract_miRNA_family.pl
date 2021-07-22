#!/usr/bin/perl

my @files=glob "/lustre1/leili_pkuhpc/guozhl/MITEs/PMITEdb/Oryza_sRNA/rpm/*rpm";
my $except_oryza_file="./Except_oryza_miRNA_family.lst";
open ORYZA,$except_oryza_file or die "fileOpenError: unable to open $except_oryza_file\n";
while(<ORYZA>){
    chomp;
    $oryza_hash{$_}=1;
}
close ORYZA;

foreach my $file(@files){
    my $outfile=$file;
    my $species;
    $outfile=~s/.rpm$/\.familylst/g;
    if($file=~/([^\/]*)\.rpm$/){
        $species=$1;
    }
    $species_hash{$species}=1;
    
    open FILE,$file or die "fileOpenError: unable to open $file\n";
    while(<FILE>){
        chomp;
        my $line=$_;
        my @line=split "\t",$line;
        my($miRNA,$rpm)=@line;
        my $family;
        if($rpm){
            if($miRNA=~/(miRN?\d+)/){
                $family=$1;
                $family_hash{$family}=1;
                $$family{$species}=1;
            }
            else{die "regular expression error for miRNA: $miRNA\n";}
            #print "$miRNA\t$family\n";
        }
    }
    close FILE;
}

### print header line ###
print "miRNA_family";
foreach my $i(sort keys %species_hash){
    print "\t$i";
}
print "\n";

###  ###
foreach my $miRNA_family(sort keys %family_hash){
    if(exists $oryza_hash{$miRNA_family}){
        print "$miRNA_family";
        foreach my $species(sort keys %species_hash){
            my $out=0;
            if(exists $$miRNA_family{$species}){
                $out=1;
            }
            print "\t$out"
        }
        print "\n";
    }
}