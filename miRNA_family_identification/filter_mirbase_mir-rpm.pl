#!/bin/perl

# filter mirbase pre
# correct position


$pre_struc = shift;	#pre二级结构
$pre_aln = shift;	#reads在pre上的map结果（额外最后后缀上lib名称）
$known_mir_1line = shift;	# mirbase mature序列文件
$reads_count = shift;	# 各lib的reads计数文件
$pre_bed = shift; #pre的bed文件，由blat得到

%pre_struc_hash=();
%pre_seq_hash=();
%aln_hash=();
%known_mature_hash=();
%reads_count_hash=();

%bed_hash=();

%max_hash=();

open(STRUC,$pre_struc);
while(<STRUC>)
{
	if(/^>/)
	{	chomp;s/^>//;/^(\S+)/; my $id=$1;
		my $seq=<STRUC>;chomp($seq); $seq=~s/U/T/g;
		my $struc=<STRUC>;chomp($struc);$struc=~s/\s\([\-\d\.]+\)$//;
		
		$pre_struc_hash{$id} = $struc; $pre_seq_hash{$id} = $seq;
	}
}
close(STRUC);

open(READS,$reads_count);
while(<READS>)
{
	chomp; my @a=split "\t"; $reads_count_hash{$a[1]}=$a[2];
}
close(READS);

open(BEDFILE,$pre_bed);
while(<BEDFILE>)
{
	chomp;
	my @a=split "\t";
	$bed_hash{$a[3]}{chr}=$a[0];
	$bed_hash{$a[3]}{pre_beg}=$a[1];
	$bed_hash{$a[3]}{pre_end}=$a[2];
	$bed_hash{$a[3]}{str}=$a[5];
}
close(BEDFILE);

open(ALN,$pre_aln);
while(<ALN>)
{
	chomp; my @a=split "\t"; my $lib=$a[8];
		if($a[1] eq "-"){next;}
		if(length($a[4])<19 || length($a[4])>25){next;}
		$a[0]=~/_x(\d+)/; my $reads_num=$1; my $rpm=$reads_num/$reads_count_hash{$lib}*1000000;
		$a[0]=~s/_x$reads_num/_x$rpm/;
		
		$aln_hash{$a[2]}{$lib}{$a[0]}{posi}=$a[3];
		$aln_hash{$a[2]}{$lib}{$a[0]}{len}=length($a[4]);
		$aln_hash{$a[2]}{$lib}{$a[0]}{rpm}=$rpm;
		
		if($aln_hash{$a[2]}{$lib}{$a[0]}{rpm} > $max_hash{$a[2]}{rpm} && (length($a[4])>=20 && length($a[4]) <= 24))
		{	$max_hash{$a[2]}{rpm}=$aln_hash{$a[2]}{$lib}{$a[0]}{rpm}; 
			$max_hash{$a[2]}{lib}=$lib; 
			$max_hash{$a[2]}{reads}=$a[0];
		}
}
close(ALN);
	#miR1030-isoform2	+	scaffold2896=GSM1770483-scaffold2896_0=1..22	1	TCTGCATCTGCACCTGCACCA	IIIIIIIIIIIIIIIIIIIII	2	1:T>C
	
open(MIR,$known_mir_1line);
while(<MIR>)
{
	chomp; my @a=split "\t";
	$a[0]=~s/^>//; $a[0]=~s/miR/MIR/;
	$a[0]=~/(\w+-MIR\d+\w*)/; my $id=$1;
	$a[0]=~/^(\S+)/; my $mir = $1;
	my $index = index($pre_seq_hash{$id}, $a[1]);
	$known_mature_hash{$id}{$mir}{posi}=$index;
	$known_mature_hash{$id}{$mir}{len}=length($a[1]);
}
close(MIR);
	
	
	
foreach my $pre(sort keys %pre_seq_hash)
{
	my $good_mir="T";
	
	unless($max_hash{$pre}{reads}){print STDERR "$pre\tNo reads in existing libs!\n"; next;}
	#unless($max_hash{$pre}{reads}){print STDERR "$pre\n"; next}
	
	my $mature_beg=$aln_hash{$pre}{$max_hash{$pre}{lib}}{$max_hash{$pre}{reads}}{posi};
	my $mature_len=$aln_hash{$pre}{$max_hash{$pre}{lib}}{$max_hash{$pre}{reads}}{len};
	my $mature_end=$aln_hash{$pre}{$max_hash{$pre}{lib}}{$max_hash{$pre}{reads}}{posi} + $aln_hash{$pre}{$max_hash{$pre}{lib}}{$max_hash{$pre}{reads}}{len} -1;
	my $mature_lib=$max_hash{$pre}{lib};
	
	my $mature_seq = substr($pre_seq_hash{$pre},$mature_beg,$mature_len);
	my $mature_struc = substr($pre_struc_hash{$pre},$mature_beg,$mature_len);
	
	
	%hash_bp=();
	fill_structure($pre);
	my($star_beg, $star_end)=find_star($mature_beg,$mature_end-2,$mature_len);
	my $star_seq=substr($pre_seq_hash{$pre},$star_beg,$star_end-$star_beg+1);
    my $star_struc=substr($pre_struc_hash{$pre},$star_beg,$star_end-$star_beg+1);
	
	# 筛选duplex结构
		unless(($mature_struc =~ /[\(\.]+/ && $star_struc =~ /[\)\.]+/) || ($mature_struc =~ /[\)\.]+/ && $star_struc =~ /[\(\.]+/)){print STDERR "$pre\t$mature_lib! Bifurcation or additional stemloop in duplex! mature_struc=$mature_struc, star_struc=$star_struc, mat_seq=$mature_seq, star_seq=$star_seq, mat_beg=$mature_beg,mat_end=$mature_end, star_beg=$star_beg,star_end=$star_end\n"; next;}
	
	# 筛选duplex错配数量
	my @mature_mismatch=$mature_struc =~/\./g;
	my @star_mismatch=$star_struc =~/\./g;
		if((@mature_mismatch >5) or (@star_mismatch >5) or (abs(@mature_mismatch - @star_mismatch)>3)){print STDERR "$pre\t$mature_lib! Too many mismatches in duplex! mature_len=$mature_len, mature_beg=$mature_beg, mature_end=$mature_end, $max_hash{$pre}{reads}, mature_struc=$mature_struc, star_struc=$star_struc, mat_seq=$mature_seq, star_seq=$star_seq\n"; next;}
		
	# 筛选pre长度上限
		if(abs($mature_beg-$star_end) > 300 || abs($star_beg-$mature_end) > 300){print STDERR "$pre\t$mature_lib! Too long precursor!\n"; next;}
		
	# 筛选reads比例
	my $total_rpm = 0; my $star_rpm = 0; $mature_rpm = 0;
		# 以mature+-1和star+-1作为stemloop的边界，只计算内部reads数量
		my $f_boundary=0; my $t_boundary=0;
		if($mature_beg > $star_beg){$f_boundary=$star_beg-1;$t_boundary=$mature_end+1}
		else{$f_boundary=$mature_beg-1;$t_boundary=$star_end+1}
	for my $reads(sort keys %{$aln_hash{$pre}{$mature_lib}})
	{
		# 只计算stemloop上的reads（mature/star之间的reads）
		unless($aln_hash{$pre}{$mature_lib}{$reads}{posi} >= $f_boundary && $aln_hash{$pre}{$mature_lib}{$reads}{posi} + $aln_hash{$pre}{$mature_lib}{$reads}{len}-1 <= $t_boundary){next;}
		
			#计算两端都在+-1nt以内的reads
		if( abs($aln_hash{$pre}{$mature_lib}{$reads}{posi} - $mature_beg) <=1 && abs($aln_hash{$pre}{$mature_lib}{$reads}{posi}+$aln_hash{$pre}{$mature_lib}{$reads}{len}-1 - $mature_end) <=1 ){$mature_rpm+=$aln_hash{$pre}{$mature_lib}{$reads}{rpm};}
		if( abs($aln_hash{$pre}{$mature_lib}{$reads}{posi} - $star_beg) <=1 && abs($aln_hash{$pre}{$mature_lib}{$reads}{posi}+$aln_hash{$pre}{$mature_lib}{$reads}{len}-1 - $star_end) <=1 ){$star_rpm+=$aln_hash{$pre}{$mature_lib}{$reads}{rpm};}
		$total_rpm+=$aln_hash{$pre}{$mature_lib}{$reads}{rpm};
	}
	
	# 筛选mature+star超过75%的reads
		if($mature_rpm + $star_rpm <= 0.75 * $total_rpm){my $m_reads=$mature_rpm/1000000*$reads_count_hash{$mature_lib}; my $s_reads=$star_rpm/1000000*$reads_count_hash{$mature_lib};my $t_reads=$total_rpm/1000000*$reads_count_hash{$mature_lib}; print STDERR "$pre\t$mature_lib! Too few reads in duplex! ratio=".(($mature_rpm + $star_rpm)/$total_rpm).", mat_reads=$m_reads, star_reads=$s_reads, total_reads=$t_reads\n"; next;}
	
	# 筛选23/24nt的mir的mature-rpm阈值
		if($mature_len == 23 || $mature_len == 24){if($mature_rpm < 10){print STDERR "$pre\t$mature_lib! Too few reads in mature for 23/24nt miRNA! mature_reads=".$mature_rpm."\n"; next;}}
	
	# 判断当前的mature是否与某个已有的注释一致
	my $tag = "";
	for my $mir(sort keys %{$known_mature_hash{$pre}})
	{
		if($known_mature_hash{$pre}{$mir}{posi} == $mature_beg && $known_mature_hash{$id}{$a[0]}{len} == $mature_len)
		{$tag = $mir;}
	}
	
	
	# 整理输出结果需要的mature/pre的位置信息
	my $out_mature_beg=0;
	my $out_mature_end=0;
	my $out_pre_beg=0;
	my $out_pre_end=0;
	my $out_pre_seq="";
	
	
	# 修正bug，考虑正负链时的不同情况
	if((($mature_beg < $star_beg) && ($bed_hash{$pre}{str} eq "+"))){$out_pre_beg = $bed_hash{$pre}{pre_beg} + $mature_beg; $out_pre_end = $bed_hash{$pre}{pre_beg} + $star_end; $out_pre_seq=substr($pre_seq_hash{$pre},$mature_beg,$star_end-$mature_beg+1);$out_mature_beg = $bed_hash{$pre}{pre_beg} + $mature_beg;$out_mature_end = $bed_hash{$pre}{pre_beg} + $mature_end;}
	
	elsif((($mature_beg > $star_beg) && ($bed_hash{$pre}{str} eq "-"))){$out_pre_beg = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_end); $out_pre_end = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg}  - $star_beg); $out_pre_seq=substr($pre_seq_hash{$pre},$star_beg,$mature_end-$star_beg+1);$out_mature_beg = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_end);$out_mature_end = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_beg);}
	
	elsif((($mature_beg < $star_beg) && ($bed_hash{$pre}{str} eq "-"))){$out_pre_beg = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $star_end); $out_pre_end = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_beg); $out_pre_seq=substr($pre_seq_hash{$pre},$mature_beg,$star_end-$mature_beg+1);$out_mature_beg = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_end);$out_mature_end = $bed_hash{$pre}{pre_beg} + ($bed_hash{$pre}{pre_end}-$bed_hash{$pre}{pre_beg} - $mature_beg);}
	
	elsif((($mature_beg > $star_beg) && ($bed_hash{$pre}{str} eq "+"))){$out_pre_beg = $bed_hash{$pre}{pre_beg} + $star_beg; $out_pre_end = $bed_hash{$pre}{pre_beg} + $mature_end; $out_pre_seq=substr($pre_seq_hash{$pre},$star_beg,$mature_end-$star_beg+1);$out_mature_beg = $bed_hash{$pre}{pre_beg} + $mature_beg;$out_mature_end = $bed_hash{$pre}{pre_beg} + $mature_end;}
	
	
	$pre =~ /MIR(\d+)/; my $conserve_tag = $1;
	
	
	#10rpm以上的pre认为confidence level = 3
	my $level=2; if($total_rpm>=10){$level=3;}
	
	if($good_mir eq "T"){print "$bed_hash{$pre}{chr}\t$bed_hash{$pre}{str}\t$max_hash{$pre}{reads}\t${mature_lib}-${pre}-$max_hash{$pre}{reads}\t$out_mature_beg\.\.$out_mature_end\t$out_pre_beg\.\.$out_pre_end\t$mature_seq\t$out_pre_seq\tconserved=$conserve_tag\t$level\tmirbase_good=${mature_lib}-${pre}-$max_hash{$pre}{reads}=$pre\n";}
}
	
	
	
sub fill_structure{

    #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired
	
	my $id=shift;
	
    my $struct=$pre_struc_hash{$id};
    my $lng=length $struct;

    #local stack for keeping track of basepairings
    my @bps;

    for(my $pos=1;$pos<=$lng;$pos++){
	my $struct_pos=excise_struct($struct,$pos,$pos,"+");

	if($struct_pos eq "("){
	    push(@bps,$pos);
	}

	if($struct_pos eq ")"){
	    my $pos_prev=pop(@bps);
	    $hash_bp{$pos_prev}=$pos;
	    $hash_bp{$pos}=$pos_prev;
	}
    }
    return;
}

sub find_star{

    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions

    #the -2 is for the overhang
    #my $mature_beg=$mature_beg;
    my($mature_beg,$mature_end1,$mature_lng)=@_;
	
	#my $mature_end1=$mature_beg + $mature_len -1 -2;
    #my $mature_lng=$mature_end1-$mature_beg+1;

    my $offset_star_beg=0;
    my $offset_beg=0;
		
    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(!$offset_star_beg and $offset_beg<$mature_lng){
	if($hash_bp{$mature_end1-$offset_beg}){
	    $offset_star_beg=$hash_bp{$mature_end1-$offset_beg};
	}else{
	    $offset_beg++;
	}
    }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg=$offset_star_beg-$offset_beg;

    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while(!$offset_star_end and $offset_end<$mature_lng){
	if($hash_bp{$mature_beg+$offset_end}){
	    $offset_star_end=$hash_bp{$mature_beg+$offset_end};
	}else{
	    $offset_end++;
	}
    }
    #the +2 is for the overhang
    my $star_end=$offset_star_end+$offset_end+2;

    return($star_beg,$star_end);
}

sub excise_struct{

    #excise sub structure

    my($struct,$beg,$end,$strand)=@_;
    my $lng=length $struct;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $subject_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($struct)){return 0;}

    #if excising relative to minus strand, positions are reversed
    if($strand eq "-"){($beg,$end)=rev_pos($beg,$end,$lng);}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_struct=substr($struct,$beg,$end-$beg+1);
 
    return $sub_struct;
}

sub rev_pos{

#   The blast_parsed format always uses positions that are relative to the 5' of the given strand
#   This means that for a sequence of length n, the first nucleotide on the minus strand base pairs with
#   the n't nucleotide on the plus strand

#   This subroutine reverses the begin and end positions of positions of the minus strand so that they
#   are relative to the 5' end of the plus strand	
   
    my($beg,$end,$lng)=@_;
    
    my $new_end=$lng-$beg+1;
    my $new_beg=$lng-$end+1;
    
    return($new_beg,$new_end);
}

	