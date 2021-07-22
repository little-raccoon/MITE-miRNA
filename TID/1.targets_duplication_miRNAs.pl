#!/usr/bin/perl

my $work_dir="/home/leili_pkuhpc/lustre1/guozhl/PMITE_dir.link/target_duplication/workdir/";
my $transcript_file=shift;
my $transcript_gff_file=shift;
my $premiRNA_file=shift;
my $premiRNA_gff_file=shift;

##### mapping pre-miRNA using bowite #####

my $tmp1=$transcript_file;
$tmp1=~s/^.*\///g;
my $bowtie2_index=$work_dir."bowtie2_index/".$tmp1;
my $bed_file=$work_dir."bed_files/".$tmp1.".bed";
#print "$bowtie2_index\n";
system("bowtie2-build $transcript_file $bowtie2_index");
system("bowtie2 -f $premiRNA_file -x $bowtie2_index | samtools view -b | bamToBed > $bed_file");


##### mask transcript files #####

my $masked_transcript=$work_dir."masked_transcript/".$tmp1.".masked.fa";
system("perl /home/leili_pkuhpc/lustre1/guozhl/script/maskBedSeq.pl $transcript_file $bed_file $masked_transcript");


##### blast premiRNA and masked_transcript #####

my $blastdb=$masked_transcript.".blastdb";
my $blastout=$work_dir."blastout/".$tmp1.".blastout";

system("makeblastdb -in $masked_transcript -out $blastdb -dbtype nucl");
system("blastn -query $premiRNA_file -db $blastdb -out $blastout -evalue 1e-10 -outfmt 6 -num_alignments 1");

##### filtered blast output by >50% of miRNA length
my $blastout_filter1_file=$blastout.".filterdByLength";

open MIRNA,$premiRNA_file or die "fileOpenError: unable to open $premiRNA_file\n";
while(<MIRNA>){
    chomp;
    my $line=$_;
    if($.==1){
        $gene=$line;
        $gene=~s/^>//g;
    }
    else{
        if($line=~/^>/){
            my $premiRNA_len=length($seq);
            $premiRNA_len_hash{$gene}=$premiRNA_len;
            $gene=$line;
            $gene=~s/^>//g;
            $seq="";
        }
        else{
            $line=~s/\s*//g;
            $seq.=$line;
        }
    }
}
close MIRNA;

# foreach my $i(sort keys %premiRNA_len_hash){
    # my $k=$premiRNA_len_hash{$i};
    # print "$i\t$k\n";
# }

open BLASTOUT,$blastout or die "fileOpenError: unable to open $blastout\n";
open OUT1,">".$blastout_filter1_file;
while(<BLASTOUT>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $miRNA=$line[0];
    my($start,$end)=@line[6..7];
    my $matched_len=$end-$start+1;
    if(exists $premiRNA_len_hash{$miRNA}){
        my $half_len=$premiRNA_len_hash{$miRNA}*0.5;
        if($matched_len>=$half_len){
            print OUT1 "$line\n";
        }
    }
    else{
        print STDERR "error1: $miRNA not in premiRNA_len_hash\n";
    }
}
close BLASTOUT;
close OUT1;


##### filtered blast output by intersecting transcript GFF3
my $intersect_file=$work_dir."intersect/".$tmp1.".intersect";
my $blastout_filter2_file=$blastout.".filterdByOverlap";

system("bedtools intersect -a $premiRNA_gff_file -b $transcript_gff_file -wa -wb > $intersect_file");
open INTERSECT,$intersect_file or die "fileOpenError: unable to open $intersect_file\n";
while(<INTERSECT>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my $miRNA=$line[8];
    my $mark=$line[11];
    if($mark eq "mRNA"){
        if($line[17]=~/Name=([^;]*);/){
            my $gene=$1;
            #print "$gene\n";
            my $tmp2=$miRNA."_".$gene;
            $overlaped_hash{$tmp2}=1;
        }
        else{
            print STDERR "error2: $intersect_file not contain Name=\n";
        }
    }
}
close INTERSECT;

open OUT2,">".$blastout_filter2_file;
open BLASTOUT_FILTED1,$blastout_filter1_file or die "fileOpenError: unable to open $blastout_filter1_file\n";
while(<BLASTOUT_FILTED1>){
    chomp;
    my $line=$_;
    my @line=split "\t",$line;
    my($miRNA,$gene)=@line;
    my $query=$miRNA."_".$gene;
    if(!exists $overlaped_hash{$query}){
        print OUT2 "$line\n";
    }
}
close BLASTOUT_FILTED1;
close OUT2;
