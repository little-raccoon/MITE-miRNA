#!/bin/bash

# Merge all miRDP2 output files for one species and output basic-info, cluster, expression tables for PmiREN.
# version: v3t4


prediction_folder=$1	# collect all *_filter_P_prediction files of miRDP2 results into this folder
reads_count_file=$2		# reads number count of each sRNA file (e.g. Solanum_lycopersicum\tGSM100001-WT_1\t101165059)
chr_length_file=$3		# length of each chromosome/scaffold/contig of the genome file (e.g. Solanum_lycopersicum\tSL3.0ch00\t20852292)
species=$4				# species name (no space! same with reads_count_file and chr_length_file, e.g. Solanum_lycopersicum)
genome_file=$5			# genome file used by miRDP2 prediction
uniq_folder=$6			# collect all input reads files (e.g. >read1_x100\nATCG\n) into this folder
miRDP2_out=$7			# collect miRDP2 output folders into this folder
output_folder=$8		# output folder path/name


# other files

src=${0%%/merge_and_rank.bash}

blat_file=$src/blat
plant_isoform=$src/plant_isoform.fa
plant_hairpin=$src/plant-hairpin.1line
plant_mature=$src/plant-mature.1line
PNRD_stemloop=$src/PNRD_stemloop.1line
filter_mirbase_perl=$src/filter_mirbase_mir-rpm.pl


mkdir $output_folder

#预处理
	#按照GSM-id筛选
	for i in `ls $prediction_folder/*_P_prediction`; do perl -e '$file=shift;$id=shift;open(FILE,$file);open(ID,$id);$species=shift; $file=~/((?:GSM|SRR)\d+\S+)_filter_P_prediction/;$name=$1; while(<ID>){chomp;@a=split "\t";if($a[0] eq $species){$hash{$a[1]}=$a[2];}} if($hash{$name}=~/\d+/){ $series=1; while(<FILE>){@a=split "\t"; s/$a[3]/${name}-$a[3]-$a[2]-$series/; $a[2]=~/_x(\d+)/;$reads=$1; my $rpm=$reads/$hash{$name}*1000000; s/_x$reads/_x$rpm/g; print; $series++} }' $i $reads_count_file $species >> $output_folder/all-mod_fp_prediction; done
		
#合并并去重所有lib的输出结果
	perl -e 'while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_beg+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$pre_beg\t$pre_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $output_folder/all-mod_fp_prediction > $output_folder/all-mod_fp_prediction-tmp.bed
	
	perl -e 'while(<>){chomp;@a=split; $hash{$a[3]}{chr}=$a[0]; $a[0]=~/SL3\.0ch(\d+)/;$hash{$a[3]}{chrid}=$1; $hash{$a[3]}{beg}=$a[1];$hash{$a[3]}{end}=$a[2];$hash{$a[3]}{str}=$a[5];  }  foreach $k(sort {$hash{$a}{chr} cmp $hash{$b}{chr} || $hash{$a}{beg}<=>$hash{$b}{beg} || $hash{$a}{end}<=>$hash{$b}{end} } keys %hash){print "$hash{$k}{chr}\t$hash{$k}{beg}\t$hash{$k}{end}\t$k\t\.\t$hash{$k}{str}\n"} ' $output_folder/all-mod_fp_prediction-tmp.bed > $output_folder/all-mod_fp_prediction-tmp.sort.bed
	
	# 存在正负链都预测为miRNA的情况
	bedtools merge -d 0 -c 4,5,6 -o collapse,distinct,distinct -i $output_folder/all-mod_fp_prediction-tmp.sort.bed > $output_folder/merge_output.bed
	
	#这里筛选掉了只在一个lib中出现，没有跟其他mir给merge到一起的单独的mir
	perl -e '$merge=shift;$all=shift; open(MERGE,$merge);open(ALL,$all); while(<MERGE>){chomp;my @a=split "\t"; my @b=split ",",$a[3];if($#b==0){next;} my $highest=0;my $high_id=""; foreach my $id(@b){$id=~/_x(\d+\.?\d*)-/;if($1>$highest){$high_id=$id;$highest=$1;}} $high_hash{$high_id}="T";} while(<ALL>){$seq=$_;my @a=split "\t";if($high_hash{$a[3]} eq "T"){print $seq;}}' $output_folder/merge_output.bed $output_folder/all-mod_fp_prediction > $output_folder/all-mod_fp_prediction-nr
	
#取出预测到的mature+-1的序列
	perl -e 'while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_beg+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$mature_beg\t$mature_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $output_folder/all-mod_fp_prediction-nr > $output_folder/all-mod_fp_prediction-nr-mature.bed
	perl -e '$scaffold=shift;open(SCAFFOLD,$scaffold); %hash=();  $name=<SCAFFOLD>; chomp($name);$name=~s/^>(\S+).*$/$1/; $seq=();  while(<SCAFFOLD>)  {chomp;  s/\s+$//;  if(m/^[A-Za-z\*]+$/){$seq.=$_;}   elsif(m/^>/){$hash{$name}=$seq; $name=$_;$name=~s/^>(\S+).*$/$1/; $seq="";} elsif(m/^$/){next;}  else{print STDERR "err! $_\n";} } $hash{$name}=$seq; $bed=shift;open(BED,$bed); $upstream=shift|0;$downstream=shift|0; while(<BED>){chomp;my @t=split "\t"; if(/^#/){next;} my $ori_seq=substr($hash{$t[0]},$t[1],$t[2]-$t[1]);if($t[5] eq "-"){$ori_seq=~tr [ATCGatcg] [TAGCtagc];$ori_seq=reverse($ori_seq);} if($t[5] eq "+"){$beg=$t[1]-$upstream>=0?($t[1]-$upstream):0;if($beg==0){$offset=$upstream-$t[1];}else{$offset=0;} } else{$beg=$t[1]-$downstream>=0?($t[1]-$downstream):0;if($beg==0){$offset=$downstream-$t[1];}else{$offset=0;}} $length=$upstream+$downstream+$t[2]-$t[1]-$offset; $out=substr($hash{$t[0]},$beg,$length);if($t[5] eq "-"){$out=~tr [ATCGatcg] [TAGCtagc];$out=reverse($out);} $region_beg=index($out,$ori_seq);$region_end=$region_beg+length($ori_seq); print ">$t[0]=$t[3]=$region_beg..$region_end\t$ori_seq\n$out\n"; } ' $genome_file $output_folder/all-mod_fp_prediction-nr-mature.bed 1 1 > $output_folder/all-mod_fp_prediction-nr.mature_1nt.fa
	

#用plant_isoform文件注释预测结果
	mkdir $output_folder/index
	bowtie-build -f $output_folder/all-mod_fp_prediction-nr.mature_1nt.fa $output_folder/index/mature-1nt
	bowtie -a -v 2 $output_folder/index/mature-1nt -f $plant_isoform > $output_folder/plant_isofomr-predicted_mature-v2.aln
	
	#如果对应多个已知mir，取序号较小的家族
	perl -e '$file=shift;$id=shift;open(FILE,$file);open(ID,$id); while(<ID>){chomp;@a=split "\t"; @b=split "=",$a[2];$a[0]=~/miR(\d+)-/;$hash{$b[1]}{$1} ="T";} while(<FILE>){chomp;$seq=$_; @a=split "\t";if(keys %{$hash{$a[3]}} > 1){my $a = join ",", sort {$a<=>$b} keys %{$hash{$a[3]}};$a.=","; if($a=~/156,/ && ($a=~/159,/ || $a=~/319,/)){$a=~s/156,//;} my @tmp = split ",",$a;  print STDERR "$a[3]\t$a\n"; print "$seq\tconserved=$tmp[0]\n"} elsif(keys %{$hash{$a[3]}} == 1){my @tmp=keys %{$hash{$a[3]}};print "$seq\tconserved=$tmp[0]\n"} else{print "$seq\tnon_conserved\n";} }' $output_folder/all-mod_fp_prediction-nr $output_folder/plant_isofomr-predicted_mature-v2.aln > $output_folder/all-mod_fp_prediction-nr-annotated 2> $output_folder/multimir_pre
	
	
#筛选非保守的mir
	perl -e 'while(<>){chomp;$seq=$_;@a=split "\t";  if($a[9]=~/non_conserved/){my @tmp=split ",",$a[8];if(length($a[6])>=21 && length($a[6])<=22){if(length($a[7]) >= 56){print $seq}} else{$a[2]=~/_x(\d+\.?\d*)/; if($#tmp>=1 && length($a[7]) >= 56 && $1>=10){print "$_\t3\tmiRDP2\n"}} }  else{print "$seq\t3\tmiRDP2=$a[3]\n"}}' $output_folder/all-mod_fp_prediction-nr-annotated > $output_folder/ORI-filtered-all-nr-anno
	
#和mirbase结果合并
	perl -e 'while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_beg+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$pre_beg\t$pre_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $output_folder/ORI-filtered-all-nr-anno > $output_folder/TMP-filtered-all-nr-anno.bed

	
#筛选mirbase-mir，并加入filtered-all-nr-anno之中
	#和mirbase数据合并 - part1
	#mirbase pre seq BLAT to genome
		perl -e '$file=shift;$species=shift;open(FILE,$file); $species=~/^(\w)\w*_(\w\w)/;$abbv="$1$2";$abbv=lc($abbv); while(<FILE>){if(/^>$abbv/){s/\t/\n/;print} } ' $plant_hairpin $species > $output_folder/mirbase_pre.fa
		$blat_file -noHead $genome_file $output_folder/mirbase_pre.fa $output_folder/BLAT_out-mirbase.psl
	
	#筛选BLAT结果
		#psl文件先转变成bed文件；之后找overlap；之后对于每个query输出一个最好且不与同家族已输出的subject重叠的subject
		#最好的标准是query里不能有gap，精确匹配和匹配长度都尽量等于query长度，subject中尽量没有gap/短gap
		perl -e 'while(<>){chomp;@a=split "\t";print "$a[13]\t$a[15]\t$a[16]\t$a[9]\tiden=$a[0],alnlen=".($a[12]-$a[11]).",qlen=$a[10],subgap=$a[6],sgaplen=$a[7]\t$a[8]\n";}' $output_folder/BLAT_out-mirbase.psl > $output_folder/BLAT_out-mirbase.all_line.bed
		bedtools intersect -a $output_folder/BLAT_out-mirbase.all_line.bed -b $output_folder/BLAT_out-mirbase.all_line.bed -wa -wb -f 0.8 > $output_folder/BLAT_out-mirbase.all_line.f0.8.inter
		
		#筛选psl的align时有预筛选，如果2*qlen-iden-alnlen < 15或0.1*qlen才会考虑，否则认为没有较好的匹配
		perl -e '$bed=shift;$inter=shift;open(BED,$bed);open(INTER,$inter); while(<INTER>){chomp;my @a=split "\t"; my $l1="$a[3]=$a[0]=$a[1]=$a[2]=$a[5]";my $l2="$a[9]=$a[6]=$a[7]=$a[8]=$a[11]"; $inter_hash{$l1}{$l2}="T";$inter_hash{$l2}{$l1}="T";} while(<BED>){chomp;my @a=split "\t";my $id="$a[3]=$a[0]=$a[1]=$a[2]=$a[5]";my $fam=$1 if($a[3]=~/(MIR\d+)/);my ($iden,$alnlen,$qlen,$subgap,$sgaplen)=($1,$2,$3,$4,$5) if($a[4]=~/iden=(\d+),alnlen=(\d+),qlen=(\d+),subgap=(\d+),sgaplen=(\d+)/); $hash{$id}{iden}=$iden;$hash{$id}{alnlen}=$alnlen;$hash{$id}{qlen}=$qlen;$hash{$id}{subgap}=$subgap;$hash{$id}{sgaplen}=$sgaplen;$hash{$id}{line}=$_; $name{$fam}{$a[3]}{$id}=$id;} my %out=(); foreach my $fam(sort keys %name){foreach my $query(sort keys %{$name{$fam}}){my @tmp=sort {($hash{$name{$fam}{$query}{$a}}{qlen}*2-$hash{$name{$fam}{$query}{$a}}{alnlen}-$hash{$name{$fam}{$query}{$a}}{iden}) <=> ($hash{$name{$fam}{$query}{$b}}{qlen}*2-$hash{$name{$fam}{$query}{$b}}{alnlen}-$hash{$name{$fam}{$query}{$b}}{iden}) || $hash{$name{$fam}{$query}{$b}}{subgap}<=>$hash{$name{$fam}{$query}{$b}}{subgap} || $hash{$name{$fam}{$query}{$b}}{sgaplen}<=>$hash{$name{$fam}{$query}{$b}}{sgaplen}} keys %{$name{$fam}{$query}}; my $count=0; my $tag="F"; my $score0=$hash{$name{$fam}{$query}{$tmp[0]}}{qlen}*2-$hash{$name{$fam}{$query}{$tmp[0]}}{alnlen}-$hash{$name{$fam}{$query}{$tmp[0]}}{iden};if($score0<=15 || $score0<=0.1*$hash{$name{$fam}{$query}{$tmp[0]}}{qlen}){while($count<=$#tmp){ foreach my $o(sort keys %{$out{$fam}}){if($inter_hash{$name{$fam}{$query}{$tmp[$count]}}{$o} eq "T"){my $now=$hash{$name{$fam}{$query}{$tmp[$count]}}{qlen}*2-$hash{$name{$fam}{$query}{$tmp[$count]}}{alnlen}-$hash{$name{$fam}{$query}{$tmp[$count]}}{iden};my $old=$hash{$o}{qlen}*2-$hash{$o}{alnlen}-$hash{$o}{iden};if($now>$old || ($now==$old && ($hash{$name{$fam}{$query}{$tmp[$count]}}{subgap}>=$hash{$o}{subgap} || $hash{$name{$fam}{$query}{$tmp[$count]}}{sgaplen}>=$hash{$o}{sgaplen})) ){$tag="T";last;} else{delete($out{$fam}{$o});} }} if($tag eq "T"){$count++;next;} else{last;}} } else{print STDERR "no good align!\t$hash{$tmp[0]}{line}\n";} $out{$fam}{$name{$fam}{$query}{$tmp[$count]}}="T" if($tag eq "F"); } }  foreach my $f(sort keys %out){foreach my $id(sort keys %{$out{$f}}){my @tmp=split "=",$id;print "$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[0]\tiden=$hash{$id}{iden},alnlen=$hash{$id}{alnlen},qlen=$hash{$id}{qlen},subgap=$hash{$id}{subgap},sgaplen=$hash{$id}{sgaplen}\t$tmp[4]\n";}}' $output_folder/BLAT_out-mirbase.all_line.bed $output_folder/BLAT_out-mirbase.all_line.f0.8.inter > $output_folder/mirbase-blat.bed 2> $output_folder/not_aligned_mirbase_pre
		#perl -e 'while(<>){unless(/^\d+/){next;} @a=split "\t"; if($a[17] != 1){next;} $str=$a[8];$id=$a[9];$chr=$a[13];$beg=$a[15];$end=$a[16];if($a[0]>$hash{$id}{len}){$hash{$id}{len}=$a[0];$hash{$id}{line}="$chr\t$beg\t$end\t$id\t\.\t$str\n"} } foreach my $k(sort keys %hash){print $hash{$k}{line};}' $output_folder/BLAT_out-mirbase.psl > $output_folder/mirbase-blat.bed
	
	# 筛选mirbase-pre，只保留与miRDP2结果不重叠的部分，用于筛选和合并
		bedtools intersect -a $output_folder/mirbase-blat.bed -b $output_folder/TMP-filtered-all-nr-anno.bed -wa -wb > $output_folder/TMP-mirbase-TMP_all_nr_anno.inter
		perl -e '$fasta=shift;$inter=shift;open(FASTA,$fasta);open(INTER,$inter); while(<INTER>){chomp;@a=split "\t"; $hash{$a[3]}="T";} while(<FASTA>){if(/^>(\S+)/){if($hash{$1} eq "T"){next} else{print;$seq=<FASTA>;print "$seq";}}}' $output_folder/mirbase_pre.fa $output_folder/TMP-mirbase-TMP_all_nr_anno.inter > $output_folder/non_overlap-mirbase_pre.fa
	
	mkdir $output_folder/filter_mirbase_tmp
	# 整理二级结构
		cat $output_folder/non_overlap-mirbase_pre.fa | RNAfold --noPS > $output_folder/non_overlap-mirbase_pre.struc
	# 整理aln
		bowtie-build -f $output_folder/non_overlap-mirbase_pre.fa $output_folder/filter_mirbase_tmp/non_overlap-mirbase_pre
	
	# 整理mapping结果
	for gsm in `ls $uniq_folder | grep -E "\.fa$" `
	do
		bowtie -a -v 0 $output_folder/filter_mirbase_tmp/non_overlap-mirbase_pre -f $uniq_folder/$gsm > $output_folder/filter_mirbase_tmp/$gsm.v0.aln
		perl -e '$file=shift;$id=shift;$id=~s/\.fa$//; open(FILE,$file);while(<FILE>){chomp;print "$_\t$id\n";}' $output_folder/filter_mirbase_tmp/$gsm.v0.aln $gsm >> $output_folder/filter_mirbase_tmp/all.aln
	done
	
	#筛选mirbase-mir并加入filtered-all-nr-anno
		perl $filter_mirbase_perl $output_folder/non_overlap-mirbase_pre.struc $output_folder/filter_mirbase_tmp/all.aln $plant_mature $reads_count_file $output_folder/mirbase-blat.bed > $output_folder/filtered_mirbase_mir 2> $output_folder/filter_mirbase_tmp/bad_mirbase_pre
	
	cat $output_folder/ORI-filtered-all-nr-anno > $output_folder/filtered-all-nr-anno
	cat $output_folder/filtered_mirbase_mir | perl -e 'while(<>){@a=split "\t";;if($a[0] ne "" && $a[1] ne ""){print}}' >> $output_folder/filtered-all-nr-anno
	

	#转bed用于和mirbase结果合并less
	perl -e 'while(<>){chomp;@tmp=split "\t",$_; $tmp[4]=~m/^(\d+)\.\.(\d+)$/;$mature_beg=$1; $tmp[5]=~m/^(\d+)\.\.(\d+)$/;$pre_beg=$1; if($tmp[1] eq "+"){$pre_beg-=2;$mature_beg-=2;} elsif($tmp[1] eq "-"){$pre_beg+=1;$mature_beg+=1;} $mature_len=length($tmp[6]);$pre_len=length($tmp[7]); $mature_end=$mature_beg+$mature_len;$pre_end=$pre_beg+$pre_len; $sign=""; if(($mature_beg==$pre_beg && $tmp[1] eq "+")||($mature_end==$pre_end && $tmp[1] eq "-") ){$sign="5";}elsif(($mature_beg==$pre_beg && $tmp[1] eq "-")||($mature_end==$pre_end && $tmp[1] eq "+")){$sign="3";}else{$sign="A";}  print "$tmp[0]\t$pre_beg\t$pre_end\t$tmp[3]\t$sign\t$tmp[1]\n";} ' $output_folder/filtered-all-nr-anno > $output_folder/filtered-all-nr-anno.bed

#与mirbase数据合并 - part2
	#merge
	cat $output_folder/filtered-all-nr-anno.bed >> $output_folder/tmp.bed
	perl -e 'while(<>){chomp;@a=split; $hash{$a[3]}{chr}=$a[0]; $a[0]=~/SL3\.0ch(\d+)/;$hash{$a[3]}{chrid}=$1; $hash{$a[3]}{beg}=$a[1];$hash{$a[3]}{end}=$a[2];$hash{$a[3]}{str}=$a[5];  }  foreach $k(sort {$hash{$a}{chr} cmp $hash{$b}{chr} || $hash{$a}{beg}<=>$hash{$b}{beg} || $hash{$a}{end}<=>$hash{$b}{end} } keys %hash){print "$hash{$k}{chr}\t$hash{$k}{beg}\t$hash{$k}{end}\t$k\t\.\t$hash{$k}{str}\n"} ' $output_folder/tmp.bed > $output_folder/sort-tmp.bed
	bedtools merge -d 0 -c 4,5,6 -o collapse,distinct,distinct -i $output_folder/sort-tmp.bed -s > $output_folder/merge_output-mirbase.bed
	

#与PNRD数据合并
	#PNRD pre seq BLAT to genome
	perl -e '$file=shift;$species=shift;open(FILE,$file); $species=~/^(\w)\w*_(\w\w)/;$abbv="$1$2";$abbv=lc($abbv); while(<FILE>){if(/^>$abbv/){s/\t/\n/;print} } ' $PNRD_stemloop $species > $output_folder/PNRD_pre.fa
	$blat_file -noHead $genome_file $output_folder/PNRD_pre.fa $output_folder/BLAT_out-PNRD.psl

	
	mkdir $output_folder/basic_info_cluster
#两个序列预测二级结构
	for i in `ls $miRDP2_out/`; do
	perl -e 'my $file=$ARGV[0];$file=~/miRDP2out\/([^\s\/]+)\//;$prefix=$1; while(<>){if(/^>/){s/^>/>$prefix /} print;}' $miRDP2_out/$i/${i}_structures >> $output_folder/basic_info_cluster/ori_strucs
	done
	
	#从miRDP2中间文件的结构数据中提取各mir的结构信息
	#新取struc程序，尽量保留了所有同lib同id的原始pri序列，用index判断哪个该归属于一个pre序列
	perl -e '$pre=shift;$struc=shift;open(PRE,$pre);open(STRUC,$struc); $out_struc_20nt=shift;$out_struc_pre=shift;open(OUT20,">",$out_struc_20nt);open(OUTPRE,">",$out_struc_pre);  while(<STRUC>){if(/^>(\S+)-15-0-10\s(\S+_\d+)\s/){my $l=$1;my $i=$2; my $seq=<STRUC>;chomp($seq);my $struc=<STRUC>;chomp($struc); $seq=~tr [U] [T];$struc=~s/^([\(\)\.]+)\s\(-*\d+\.\d+\)$]/$1/; $hash{$l}{$i}.="$seq=$struc,";} } while(<PRE>){chomp;@a=();@b=();$lib="";$id=""; @a=split "\t";my @b=split "=",$a[10]; if($b[1]=~/^(\S+)-($a[0]_\d+)-\S+_x[\d\.]+/){$lib=$1;$id=$2;} if($lib ne "" && $id ne ""){my ($seq,$struc)=("","");$ori_pre=$a[7]; foreach my $pri(split ",",$hash{$lib}{$id}){if(index($pri,$ori_seq)>=0){$pri=~/([^=\s]+)=([^=\s]+)/;($seq,$struc)=($1,$2);}}unless($ori_pre=~/^[ATCGUatcguNn]+$/){next;} my $index=index($seq,$ori_pre);my $length=length($ori_pre); if($index-20<0 || $index+$length+20>length($seq)){print STDERR "pri too long!$index,$length\t$seq\t$struc\t$ori_pre\t$_\n";} my $sub20_seq=substr($seq,($index-20>=0?$index-20:0),$length+40); my $sub20_struc=substr($struc,($index-20>=0?$index-20:0),$length+40); my $sub_struc=substr($struc,$index,$length);     ($sub20_struc,$sub_struc)=&clean($sub20_struc,$sub_struc);      print OUT20 ">$a[0]=$a[3]=$index..".($index+$length)."\n$sub20_seq\n$sub20_struc\n"; print OUTPRE ">$a[0]=$a[3]=$index..".($index+$length)."\n$ori_pre\n$sub_struc\n"; } else{print STDERR "bad_pattern! $id\t$b[1],$a[0]\n";} } sub clean{my @out=();foreach my $s(@_){chomp($s);my @letter=split "",$s; my $c=0;my @queue=();while($c<=$#letter){if($#queue<0 && $letter[$c] eq "\)"){$letter[$c]=".";} elsif($letter[$c] eq "\(" ){push @queue,$c;} elsif($letter[$c] eq "\)"){pop @queue;} $c++;} foreach my $n(@queue){$letter[$n]=".";} push @out,join("",@letter);} return @out}' $output_folder/filtered-all-nr-anno $output_folder/basic_info_cluster/ori_strucs $output_folder/basic_info_cluster/stemloop_20nt.struc $output_folder/basic_info_cluster/stemloop.struc 2> $output_folder/basic_info_cluster/pri_too_long.err

	#单独整理mirbase-good的结构信息
	perl -e '$filter=shift;$struc=shift;open(FILTER,$filter);open(STRUC,$struc); $error=shift;open(ERRORFILE,">",$error); while(<STRUC>){if(/^>(\S+)/){my $id=$1;$seq=<STRUC>;chomp($seq);$seq=~s/U/T/g; my $s=<STRUC>;chomp($s);$s=~s/\s\([\d\-\.]+\)$//; $hash{$id}{seq}=$seq;$hash{$id}{struc}=$s; }} while(<FILTER>){chomp;@a=split "\t";@b=split "=",$a[-1]; my $index=index($hash{$b[2]}{seq},$a[7]);my $len=length($a[7]); print STDERR ">a[0]=$a[3]=mirbase_good\n".substr($hash{$b[2]}{seq},$index,$len)."\n".substr($hash{$b[2]}{struc},$index,$len)."\n";  if($index<20 || $index+$len+20>length($hash{$b[2]}{seq})){print ERRORFILE "too long pre! $_; pri=$hash{$b[2]}{seq}\n";} {print ">$a[0]=$a[3]=mirbase_good\n".substr($hash{$b[2]}{seq},$index-20>0?$index-20:0,$len+40)."\n".substr($hash{$b[2]}{struc},$index-20>0?$index-20:0,$len+40)."\n";}}' $output_folder/filtered_mirbase_mir $output_folder/non_overlap-mirbase_pre.struc $output_folder/bad-mirbase_good_struc >> $output_folder/basic_info_cluster/stemloop_20nt.struc 2>> $output_folder/basic_info_cluster/stemloop.struc

#整理bastc-info表格	
	perl -e '$file=shift;$structwenty=shift;$strucstem=shift;$bed=shift;$chr_length=shift;$species=shift;  open(BED,$bed);while(<BED>){chomp;my @a=split "\t";$bed_hash{$a[3]}{beg}=$a[1];$bed_hash{$a[3]}{end}=$a[2]; } open(STRUCSTEM,$strucstem);while(<STRUCSTEM>){if(/^>/){my $id=$_;my $seq=<STRUCSTEM>;my $struc=<STRUCSTEM>;chomp($struc); my @b=split "=",$id;$seq=uc($seq);$seq=~s/U/T/g;$struc=~s/\s\([\d\.\-]+\)$//;$hash_stem{$b[1]}{struc}=$struc; } }open(STRUCTWENTY,$structwenty);while(<STRUCTWENTY>){if(/^>/){my $id20=$_;my $seq20=<STRUCTWENTY>;my $struc20=<STRUCTWENTY>;chomp($struc20);chomp($seq20); my @c=split "=",$id20;$seq20=uc($seq20);$seq20=~s/U/T/g;$struc20=~s/\s\([\d\.\-]+\)$//;$hash_stem_20{$c[1]}{struc20}=$struc20;$hash_stem_20{$c[1]}{seq}=$seq20;}} open(CHRLENGTH,$chr_length);while(<CHRLENGTH>){chomp;my @abc=split "\t";if($abc[0] eq $species){$chr_len{$abc[1]}=$abc[2];}}close(CHRLENGTH); open(FILE,$file);$novel_fam=1;$accession=1;%conserve=();%non_conserve=(); my @nihil=split "_",$species;$nihil[0]=~/^(\w)/;$abbv=uc($1);$nihil[-1]=~/^(\w\w)/;$abbv.=$1; while(<FILE>){chomp;my @d=split "\t";if($d[8] eq "non_conserved"){$non_conserve{$d[3]}=$_;} elsif($d[8]=~/conserved=(\d+)/){$conserve{$1}{$d[3]}=$_;}else{print STDERR "bad conservation! $d[3]\n"}} foreach my $fcon(sort keys %conserve){my $member=1;foreach my $k(sort keys %{$conserve{$fcon}}){my $a=int($member/26)>0?int($member/26)+96:0;my $b=$member%26+96;my $tag="";if($a==0){$tag=chr($b)} else{$tag=chr($a).chr($b)} my $miR_id="${abbv}-MIR${fcon}$tag";$miR_id=~s/\s//; my @tmp=split "\t", $conserve{$fcon}{$k};$species=~s/_/ /g; my $star_seq="";my $star_beg=0;my $star_end=0; my $ind=index($tmp[7],$tmp[6]);my $confidence=$tmp[9];my $source=$tmp[10]; if($ind>=0 && abs($ind)<=2){$star_seq=substr($tmp[7],-21,21)} elsif(abs($ind-(length($tmp[7])-length($tmp[6])-1))<=2){$star_seq=substr($tmp[7],0,21)} else{print STDERR "$tmp[3],$ind,".(length($tmp[7])-length($tmp[6])-1).",$tmp[1],$tmp[6],$tmp[7]\n";last} "$tmp[4]..$tmp[5]"=~/(\d+)\.\.(\d+)\.\.(\d+)\.\.(\d+)/;if($1 == $3){$star_beg=$bed_hash{$tmp[3]}{end}-20;$star_end=$bed_hash{$tmp[3]}{end};}elsif($2 == $4){$star_beg=$bed_hash{$tmp[3]}{beg}+1;$star_end=$bed_hash{$tmp[3]}{beg}+21;}$star_seq=uc($star_seq);  my $out=sprintf "$miR_id\tPmiREN%09d\tMIR${fcon}\t$species\t${abbv}-$tmp[0]\t".($bed_hash{$tmp[3]}{beg}+1-20>=1?$bed_hash{$tmp[3]}{beg}+1-20:1)."\t".($bed_hash{$tmp[3]}{end}+20<=$chr_len{$tmp[0]}?$bed_hash{$tmp[3]}{end}+20:$chr_len{$tmp[0]})."\t$tmp[1]\t$hash_stem_20{$tmp[3]}{seq}\t$hash_stem_20{$tmp[3]}{struc20}\t".uc($tmp[7])."\t$hash_stem{$tmp[3]}{struc}\t$tmp[5]\t${abbv}-miR${fcon}$tag\t$tmp[4]\t".uc($tmp[6])."\t${abbv}-miR${fcon}${tag}*\t$star_seq\t$star_beg\t$star_end\t$confidence\t$source\n", $accession; $out=~s/(\d+)\.\.(\d+)/$1\t$2/g;print $out; $member++;$accession++;} } my $member=1;foreach my $fnon(sort keys %non_conserve){ my $miR_id="${abbv}-MIRN$member";$miR_id=~s/\s//; my @tmp=split "\t", $non_conserve{$fnon};$species=~s/_/ /g; my $star_seq="";my $star_beg=0;my $star_end=0; my $ind=index($tmp[7],$tmp[6]);my $confidence=$tmp[9];my $source=$tmp[10]; if(abs($ind)<=2){$star_seq=substr($tmp[7],-21,21)} elsif(abs($ind-(length($tmp[7])-length($tmp[6])-1))<=2){$star_seq=substr($tmp[7],0,21)} "$tmp[4]..$tmp[5]"=~/(\d+)\.\.(\d+)\.\.(\d+)\.\.(\d+)/;if($1 == $3){$star_beg=$bed_hash{$tmp[3]}{end}-20;$star_end=$bed_hash{$tmp[3]}{end};}elsif($2 == $4){$star_beg=$bed_hash{$tmp[3]}{beg}+1;$star_end=$bed_hash{$tmp[3]}{beg}+21;}$star_seq=uc($star_seq);  my $out=sprintf "$miR_id\tPmiREN%09d\tMIRN$member\t$species\t${abbv}-$tmp[0]\t".($bed_hash{$tmp[3]}{beg}+1-20>=1?$bed_hash{$tmp[3]}{beg}+1-20:1)."\t".($bed_hash{$tmp[3]}{end}+20<=$chr_len{$tmp[0]}?$bed_hash{$tmp[3]}{end}+20:$chr_len{$tmp[0]})."\t$tmp[1]\t$hash_stem_20{$tmp[3]}{seq}\t$hash_stem_20{$tmp[3]}{struc20}\t".uc($tmp[7])."\t$hash_stem{$tmp[3]}{struc}\t$tmp[5]\t${abbv}-miRN$member\t$tmp[4]\t".uc($tmp[6])."\t${abbv}-miRN${member}*\t$star_seq\t$star_beg\t$star_end\t$confidence\t$source\n", $accession; $out=~s/(\d+)\.\.(\d+)/$1\t$2/g;print $out; $member++;$accession++;}' $output_folder/filtered-all-nr-anno $output_folder/basic_info_cluster/stemloop_20nt.struc $output_folder/basic_info_cluster/stemloop.struc $output_folder/filtered-all-nr-anno.bed $chr_length_file $species > $output_folder/basic_info_cluster/primary-${species}-basic_info


	perl -e '$info=shift;$inter=shift;open(INFO,$info);open(INTER,$inter); while(<INTER>){chomp;my @a=split "\t"; $inter_hash{$a[9]}.="$a[3],"; $occupied_hash{uc($a[3])}="T";} while(<INFO>){chomp;$line=$_;my @a=split "\t";my @tmp=split "=",$a[23]; if($tmp[0] eq "miRDP2"){$info_hash{$a[2]}{$tmp[1]}=$line;} elsif($tmp[0] eq "mirbase_good"){$info_hash{$a[2]}{$tmp[1]}=$line;$occupied_hash{uc($tmp[2])}="T";} else{print STDERR "bad line: $id\n";} } foreach my $fam(sort keys %info_hash){my $member=1; foreach my $id(sort keys %{$info_hash{$fam}}){if(exists($inter_hash{$id})){my @known_id=split ",",$inter_hash{$id}; if($info_hash{$fam}{$id}=~/MIRN/){$info_hash{$fam}{$id}=~/^(\S+)\t/;print "$1\t$info_hash{$fam}{$id}=$known_id[0]\n";my $o=join ",",@known_id;print STDERR "unassigned_known\t$id\t$o\n";next;} my @nouse=split "",$known_id[0];$nouse[0]=uc($nouse[0]);$known_id[0]=join "",@nouse; print "$known_id[0]\t$info_hash{$fam}{$id}=$known_id[0]\n"; if($known_id[1]=~/MIR\d+/){my $o=join ",",@known_id;print STDERR "multi_mir\t$id\t$o\n";} next;} elsif($info_hash{$fam}{$id}=~/MIRN\d+/){$info_hash{$fam}{$id}=~/^(\S+)\t/;print "$1\t$info_hash{$fam}{$id}\n";next;} my @nihil=split "\t",$info_hash{$fam}{$id};$nihil[0]=~/^(\w\w\w)-/;$s=$1; my $a=int($member/26)>0?int($member/26)+96:0;my $b=$member%26+96;my $tag="";if($a==0){$tag=$s."-".$nihil[2].chr($b)} else{$tag=$s."-".$nihil[2].chr($a).chr($b);} while($occupied_hash{uc($tag)} eq "T"){$member++;$a=int($member/26)>0?int($member/26)+96:0;$b=$member%26+96;$tag="";if($a==0){$tag=$s."-".$nihil[2].chr($b)} else{$tag=$s."-".$nihil[2].chr($a).chr($b);} } print "$tag\t$info_hash{$fam}{$id}\n";$member++; }}' $output_folder/basic_info_cluster/primary-${species}-basic_info $output_folder/TMP-mirbase-TMP_all_nr_anno.inter > $output_folder/basic_info_cluster/tmp_info 2>$output_folder/basic_info_cluster/basic_info_id.err
	cat $output_folder/basic_info_cluster/tmp_info | cut -f 1,3- | perl -e 'while(<>){chomp;@a=split "\t"; if($a[0]=~/(MIR\d+)/){$new_fam=$1; if($a[2] ne $new_fam){$a[2]=$new_fam} my $name=$a[0];$name=~s/MIR/miR/;$a[14]=$name;$a[18]=$name."*";} my $out=join "\t",@a; print "$out\n";}' | sort > $output_folder/basic_info_cluster/secondary-${species}-basic_info
	
	# 整理inter中多个pre都overlap到同一个mirbase-mir的情况
	perl -e 'while(<>){chomp;@a=split "\t";$mir_hash{$a[3]}.="$a[9]";} foreach my $mir(sort keys %mir_hash){my @tmp=split ",",$mir_hash{$mir}; if($tmp[1]=~/_x\d+/){my $o=join ",",@tmp;print "multi_prediction\t$o\t$mir\n";}}' $output_folder/TMP-mirbase-TMP_all_nr_anno.inter >> $output_folder/basic_info_cluster/basic_info_id.err

	#由于mirbase-good没有跟pnrd比对过，所以这里先从basic-info里提取包含了mirbase-good的所有pre的bed信息，去和pnrd做overlap
	perl -e 'while(<>){chomp;@a=split "\t"; @tmp=split "=",$a[23];$a[4]=~/\w\w\w-(\S+)/;$chr=$1; print "$chr\t".($a[5]-1)."\t$a[6]\t$chr=$a[5]=$a[6]=$tmp[1]\t\.\t$a[7]\n";}' $output_folder/basic_info_cluster/secondary-${species}-basic_info > $output_folder/mirdp_mirbasegood.pre.bed
	#整理PNRD-pre位置
	perl -e 'while(<>){chomp;@a=split "\t";print "$a[13]\t$a[15]\t$a[16]\t$a[9]\tiden=$a[0],alnlen=".($a[12]-$a[11]).",qlen=$a[10],subgap=$a[6],sgaplen=$a[7]\t$a[8]\n";}' $output_folder/BLAT_out-PNRD.psl > $output_folder/BLAT_out-PNRD.all_line.bed
	bedtools intersect -a $output_folder/BLAT_out-PNRD.all_line.bed -b $output_folder/BLAT_out-PNRD.all_line.bed -wa -wb -f 0.8 > $output_folder/BLAT_out-PNRD.all_line.f0.8.inter
	perl -e '$bed=shift;$inter=shift;open(BED,$bed);open(INTER,$inter); while(<INTER>){chomp;my @a=split "\t"; my $l1="$a[3]=$a[0]=$a[1]=$a[2]=$a[5]";my $l2="$a[9]=$a[6]=$a[7]=$a[8]=$a[11]"; $inter_hash{$l1}{$l2}="T";$inter_hash{$l2}{$l1}="T";} while(<BED>){chomp;my @a=split "\t";my $id="$a[3]=$a[0]=$a[1]=$a[2]=$a[5]";my $fam=$1 if($a[3]=~/(MIR\d+)/);my ($iden,$alnlen,$qlen,$subgap,$sgaplen)=($1,$2,$3,$4,$5) if($a[4]=~/iden=(\d+),alnlen=(\d+),qlen=(\d+),subgap=(\d+),sgaplen=(\d+)/); $hash{$id}{iden}=$iden;$hash{$id}{alnlen}=$alnlen;$hash{$id}{qlen}=$qlen;$hash{$id}{subgap}=$subgap;$hash{$id}{sgaplen}=$sgaplen;$hash{$id}{line}=$_; $name{$fam}{$a[3]}{$id}=$id;} my %out=(); foreach my $fam(sort keys %name){foreach my $query(sort keys %{$name{$fam}}){my @tmp=sort {($hash{$name{$fam}{$query}{$a}}{qlen}*2-$hash{$name{$fam}{$query}{$a}}{alnlen}-$hash{$name{$fam}{$query}{$a}}{iden}) <=> ($hash{$name{$fam}{$query}{$b}}{qlen}*2-$hash{$name{$fam}{$query}{$b}}{alnlen}-$hash{$name{$fam}{$query}{$b}}{iden}) || $hash{$name{$fam}{$query}{$b}}{subgap}<=>$hash{$name{$fam}{$query}{$b}}{subgap} || $hash{$name{$fam}{$query}{$b}}{sgaplen}<=>$hash{$name{$fam}{$query}{$b}}{sgaplen}} keys %{$name{$fam}{$query}}; my $count=0; my $tag="F"; my $score0=$hash{$name{$fam}{$query}{$tmp[0]}}{qlen}*2-$hash{$name{$fam}{$query}{$tmp[0]}}{alnlen}-$hash{$name{$fam}{$query}{$tmp[0]}}{iden};if($score0<=15 || $score0<=0.1*$hash{$name{$fam}{$query}{$tmp[0]}}{qlen}){while($count<=$#tmp){ foreach my $o(sort keys %{$out{$fam}}){if($inter_hash{$name{$fam}{$query}{$tmp[$count]}}{$o} eq "T"){my $now=$hash{$name{$fam}{$query}{$tmp[$count]}}{qlen}*2-$hash{$name{$fam}{$query}{$tmp[$count]}}{alnlen}-$hash{$name{$fam}{$query}{$tmp[$count]}}{iden};my $old=$hash{$o}{qlen}*2-$hash{$o}{alnlen}-$hash{$o}{iden};if($now>$old || ($now==$old && ($hash{$name{$fam}{$query}{$tmp[$count]}}{subgap}>=$hash{$o}{subgap} || $hash{$name{$fam}{$query}{$tmp[$count]}}{sgaplen}>=$hash{$o}{sgaplen})) ){$tag="T";last;} else{delete($out{$fam}{$o});} }} if($tag eq "T"){$count++;next;} else{last;}} } else{print STDERR "no good align!\t$hash{$tmp[0]}{line}\n";} $out{$fam}{$name{$fam}{$query}{$tmp[$count]}}="T" if($tag eq "F"); } }  foreach my $f(sort keys %out){foreach my $id(sort keys %{$out{$f}}){my @tmp=split "=",$id;print "$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[0]\tiden=$hash{$id}{iden},alnlen=$hash{$id}{alnlen},qlen=$hash{$id}{qlen},subgap=$hash{$id}{subgap},sgaplen=$hash{$id}{sgaplen}\t$tmp[4]\n";}}' $output_folder/BLAT_out-PNRD.all_line.bed $output_folder/BLAT_out-PNRD.all_line.f0.8.inter > $output_folder/PNRD-blat.bed 2> $output_folder/not_aligned_PNRD_pre

	cat $output_folder/PNRD-blat.bed >> $output_folder/mirdp_mirbasegood.pre.bed
	#sort bed文件
	perl -e 'while(<>){chomp;@a=split; $hash{$a[3]}{chr}=$a[0]; $a[0]=~/SL3\.0ch(\d+)/;$hash{$a[3]}{chrid}=$1; $hash{$a[3]}{beg}=$a[1];$hash{$a[3]}{end}=$a[2];$hash{$a[3]}{str}=$a[5];  }  foreach $k(sort {$hash{$a}{chr} cmp $hash{$b}{chr} || $hash{$a}{beg}<=>$hash{$b}{beg} || $hash{$a}{end}<=>$hash{$b}{end} } keys %hash){print "$hash{$k}{chr}\t$hash{$k}{beg}\t$hash{$k}{end}\t$k\t\.\t$hash{$k}{str}\n"} ' $output_folder/mirdp_mirbasegood.pre.bed > $output_folder/sorted-mirdp_mirbasegood.pre.bed
	#merge混合pre文件，寻找mir-mirbase-PNRD相互重合的情况
	bedtools merge -d 0 -c 4,5,6 -o collapse,distinct,distinct -i $output_folder/sorted-mirdp_mirbasegood.pre.bed > $output_folder/merge-sorted-mirdp_mirbasegood.pre.bed

	#继续按照有重叠的标准寻找PNRD注释和miRDP2注释的重叠情况
	#	根据重叠情况给每个mir加上mirbase和PNRD的链接
	perl -e '$pnrd=shift;$mirbase=shift;$basic=shift;open(PNRD,$pnrd);open(MIRBASE,$mirbase);open(BASIC,$basic); while(<PNRD>){chomp;my @a=split "\t";my @id=split ",",$a[3]; my @info=();my $p="";my $tag="F"; foreach my $i(@id){if($i=~/(GSM|SRR)\d+/){push @info,$i;} if($i=~/MIR[f]?\d+/){if($p eq ""){$p=$i;next;} elsif($p=~/MIRf\d+/ && $i=~/MIR\d+/){$p=$i} elsif($p=~/MIR(\d+)/){my $p_num=$1; if($i=~/MIR(\d+)/ && $1<$p_num){$p=$i}} $tag="T";} } print STDERR "multi_info! ".(join ",",@a)."\n" if($#info>=1); print STDERR "multi_PNRD! ".(join ",",@a)."\n" if($tag eq "T"); foreach my $in(@info){$pnrd_hash{$in}=$p;} } while(<MIRBASE>){if(/^>(\S+)\s+(MI\S+)/){$mirbase_hash{$1}=$2; }} while(<BASIC>){chomp;my @a=split "\t";my $link_m="";my $link_p=""; my @tmp=split "=",$a[23];if($tmp[2]=~/MIR\d+/){$tmp[2]=~/^([A-Za-z])/;my $nihil=lc($1);$tmp[2]=~s/^([A-Za-z])/$nihil/; $link_m="http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=$mirbase_hash{$tmp[2]}";} $a[4]=~/\w\w\w-(\S+)/;$chr=$1; my $id="$chr=$a[5]=$a[6]=$tmp[1]";if($pnrd_hash{$id}=~/MIR[f]?\d+/){$link_p="http://structuralbiology.cau.edu.cn/PNRD/result_detail.php?ncRNA_name=$pnrd_hash{$id}";} my $tag_m="×";my $tag_p="×";if($link_m=~/http/){$tag_m="√"} if($link_p=~/http/){$tag_p="√"} print((join "\t",@a[0..21])."\t\tThis is $a[0]\t$a[22]\t$tag_m/$tag_p\t$link_m\t$link_p\n")} ' $output_folder/merge-sorted-mirdp_mirbasegood.pre.bed $output_folder/mirbase_pre.fa $output_folder/basic_info_cluster/secondary-${species}-basic_info > $output_folder/basic_info_cluster/${species}-basic_info 2> $output_folder/basic_info_cluster/final.err
	

	
#basic_info提取mature序列
	cat $output_folder/basic_info_cluster/${species}-basic_info | cut -f 15,18 | sed 's/^/>/' | sed 's/\t/\n/' > $output_folder/basic_info_cluster/${species}.mature.fa
	
#整理cluster序列
	perl -e '$basic_info=shift;$species=shift;open(INFO,$basic_info);while(<INFO>){@a=split "\t";$hash{$a[4]}{"$a[0]=$a[1]"}{beg}=$a[5];$hash{$a[4]}{"$a[0]=$a[1]"}{end}=$a[6];$hash{$a[4]}{"$a[0]=$a[1]"}{str}=$a[7];$hash{$a[4]}{"$a[0]=$a[1]"}{acc}=$a[1];} $cluster_id=1; foreach my $chr(sort keys %hash){my @list=sort {$hash{$chr}{$a}{beg}<=>$hash{$chr}{$b}{beg} } keys %{$hash{$chr}};$num=1; while($num<=$#list){if(abs($hash{$chr}{$list[$num-1]}{end}-$hash{$chr}{$list[$num]}{beg}) <= 10000 && $hash{$chr}{$list[$num-1]}{str} eq $hash{$chr}{$list[$num]}{str}){if($cluster_hash{$cluster_id}{mir}{$list[$num-1]} ne "T"){$cluster_hash{$cluster_id}{beg}=$hash{$chr}{$list[$num-1]}{beg}; } $cluster_hash{$cluster_id}{mir}{$list[$num-1]}="T";$cluster_hash{$cluster_id}{mir}{$list[$num]}="T";$cluster_hash{$cluster_id}{chr}=$chr;$cluster_hash{$cluster_id}{str}=$hash{$chr}{$list[$num]}{str}; } else{if(exists($cluster_hash{$cluster_id}) && $cluster_hash{$cluster_id}{mir}{$list[$num-1]} eq "T"){$cluster_hash{$cluster_id}{end}=$hash{$chr}{$list[$num-1]}{end};$cluster_id++}  } $num++;} if(exists($cluster_hash{$cluster_id}) && $cluster_hash{$cluster_id}{mir}{$list[$num-1]} eq "T"){$cluster_hash{$cluster_id}{end}=$hash{$chr}{$list[$num-1]}{end};$cluster_id++;} } $species=~/^(\w)\w+_(\w\w)\w+$/;$abbv=$1.$2; foreach my $c(sort keys %cluster_hash){my $content=join ",",sort keys %{$cluster_hash{$c}{mir}}; print "${abbv}-CL$c\t$species\t$content\t$cluster_hash{$c}{chr}\t$cluster_hash{$c}{str}\t$cluster_hash{$c}{beg}\t$cluster_hash{$c}{end}\n";}' $output_folder/basic_info_cluster/${species}-basic_info $species > $output_folder/basic_info_cluster/${species}-cluster
	
	
#整理exp表格
	mkdir $output_folder/expression
	mkdir $output_folder/expression/index
	mkdir $output_folder/expression/mature
	mkdir $output_folder/expression/star
	mkdir $output_folder/expression/stemloop
	
	#从basic-info里提取star-1nt
	perl -e 'while(<>){chomp;@a=split "\t";$stem=$a[8];$mature=$a[17];$star=$a[19]; $ind1=index($stem,$mature);$out1=substr($stem,($ind1-1>=0?$ind1-1:0),length($mature)+2);$len=length($out1);print ">$len=$a[0]=$a[1]\n$out1\n"; $ind2=index($stem,$star);$out2=substr($stem,($ind2-1>=0?$ind2-1:0),length($star)+2);$len=length($out2);print STDERR ">$len=$a[0]*=$a[1]\n$out2\n";}' $output_folder/basic_info_cluster/${species}-basic_info >$output_folder/expression/mature-1nt.fa 2> $output_folder/expression/star-1nt.fa

	bowtie-build -f $output_folder/expression/mature-1nt.fa $output_folder/expression/index/mature-1nt
	bowtie-build -f $output_folder/expression/star-1nt.fa $output_folder/expression/index/star-1nt
	#bowtie-build -f $output_folder/basic_info_cluster/stemloop.fa $output_folder/expression/index/stemloop
	
	#mapping
	for gsm in `ls $uniq_folder | grep -E "\.fa$"`
	do
		bowtie -a -v 0 $output_folder/expression/index/mature-1nt -f $uniq_folder/$gsm > $output_folder/expression/mature/$gsm.mature.aln
		perl -e '$file=shift;$id=shift;open(FILE,$file);open(ID,$id);$species=shift; $file=~/((?:GSM|SRR)\d+\S+)\.fa\.mature\.aln/; $name=$1; while(<ID>){chomp;my @a=split "\t";if($a[0] eq $species){$hash{$a[1]}=$a[2];}} if($hash{$name}=~/\d+/){ while(<FILE>){my @a=split "\t"; $a[0]=~/_x(\d+)/;$reads=$1; my $rpm=$reads/$hash{$name}*1000000; $a[2]=~/^(\d+)=/;if($a[3]>2 || $a[3]+length($a[4]) < ${1}-2){next;} $a[2]=~s/^\d+=//;$out_hash{$a[2]}+=$rpm;} my @list=sort keys %out_hash;$header=join "\t",@list; $name=~/^((?:GSM|SRR)\d+)\S*/;print "$1"; foreach my $mir(@list){print "\t$mir,$out_hash{$mir}"} print "\n"; }' $output_folder/expression/mature/$gsm.mature.aln $reads_count_file $species >> $output_folder/expression/TMP-${species}-mature.exp
	
		bowtie -a -v 0 $output_folder/expression/index/star-1nt -f $uniq_folder/$gsm > $output_folder/expression/star/$gsm.star.aln
		perl -e '$file=shift;$id=shift;open(FILE,$file);open(ID,$id);$species=shift; $file=~/((?:GSM|SRR)\d+\S+)\.fa\.star\.aln/;$name=$1; while(<ID>){chomp;my @a=split "\t";if($a[0] eq $species){$hash{$a[1]}=$a[2];}} if($hash{$name}=~/\d+/){ while(<FILE>){my @a=split "\t"; $a[0]=~/_x(\d+)/;$reads=$1; my $rpm=$reads/$hash{$name}*1000000; $a[2]=~/^(\d+)=/;if($a[3]>2 || $a[3]+length($a[4]) < ${1}-2){next;} $a[2]=~s/^\d+=//;$out_hash{$a[2]}+=$rpm;} my @list=sort keys %out_hash;$header=join "\t",@list; $name=~/^((?:GSM|SRR)\d+)\S*/;print "$1"; foreach my $mir(@list){print "\t$mir,$out_hash{$mir}"} print "\n"; }' $output_folder/expression/star/$gsm.star.aln $reads_count_file $species >> $output_folder/expression/TMP-${species}-star.exp
		
	done

	perl -e '$aln=shift;$mir_info=shift;$gsm_info=shift;$species=shift;open(ALN,$aln);open(LIST,$mir_info);open(GSM,$gsm_info);@libs=(); while(<GSM>){chomp;my @a=split "\t"; $a[2]=~s/\s/_/g;$a[0]=~s/\s+$//g;$a[0]=~s/^\s+//g;$a[4]=~s/\s+$//g;$a[4]=~s/^\s+//g;$a[1]=~s/\s+$//g;$a[1]=~s/^\s+//g;$a[2]=~s/\s+$//g;$a[2]=~s/^\s+//g; if($a[2] ne $species){next;} $lib_hash{$a[0]}=$a[4]; if($a[1]=~/^SRR/){$lib_hash{$a[1]}=$a[4];} $tissue_hash{$a[4]}{count}++;} while(<ALN>){chomp;@a=split "\t";$lib=shift @a;push @libs,$lib; foreach my $t(@a){my @nihil=split ",",$t;$hash{$nihil[0]}{$lib}=$nihil[1]; $tissue_hash{$lib_hash{$lib}}{$nihil[0]}+=$nihil[1];} } my $header=join "\t",@libs;print "\t$header\n"; my $tissue=join "\t",sort keys %tissue_hash; while(<LIST>){chomp;@a=split "\t";$id="$a[14]=$a[1]";$id=~s/miR/MIR/;print "$id"; foreach my $l(@libs){my $out=$hash{$id}{$l}>0?$hash{$id}{$l}:0; print "\t$out";} print "\n";} ' $output_folder/expression/TMP-${species}-mature.exp $output_folder/basic_info_cluster/${species}-basic_info $gsm_info $species > $output_folder/expression/${species}-mature.exp
	perl -e '$aln=shift;$mir_info=shift;$gsm_info=shift;$species=shift;open(ALN,$aln);open(LIST,$mir_info);open(GSM,$gsm_info);@libs=(); while(<GSM>){chomp;my @a=split "\t"; $a[2]=~s/\s/_/g;$a[0]=~s/\s+$//g;$a[0]=~s/^\s+//g;$a[4]=~s/\s+$//g;$a[4]=~s/^\s+//g;$a[1]=~s/\s+$//g;$a[1]=~s/^\s+//g;$a[2]=~s/\s+$//g;$a[2]=~s/^\s+//g; if($a[2] ne $species){next;} $lib_hash{$a[0]}=$a[4]; if($a[1]=~/^SRR/){$lib_hash{$a[1]}=$a[4];} $tissue_hash{$a[4]}{count}++;} while(<ALN>){chomp;@a=split "\t";$lib=shift @a;push @libs,$lib; foreach my $t(@a){my @nihil=split ",",$t;$hash{$nihil[0]}{$lib}=$nihil[1]; $tissue_hash{$lib_hash{$lib}}{$nihil[0]}+=$nihil[1];} } my $header=join "\t",@libs;print "\t$header\n"; my $tissue=join "\t",sort keys %tissue_hash; while(<LIST>){chomp;@a=split "\t";$id="$a[18]=$a[1]";$id=~s/miR/MIR/;print "$id"; foreach my $l(@libs){my $out=$hash{$id}{$l}>0?$hash{$id}{$l}:0; print "\t$out";} print "\n";} ' $output_folder/expression/TMP-${species}-star.exp $output_folder/basic_info_cluster/${species}-basic_info $gsm_info $species > $output_folder/expression/${species}-star.exp

