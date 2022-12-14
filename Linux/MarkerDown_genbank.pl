#!usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;

my $genome;
my $vector;
my $overlap=100;
my $window=5000;
my $insert_vector=0;
my $insert_name;
my $conf;
my $out="Genbank_File";

GetOptions(
	'genome:s'	=>\$genome,
	'vector:s'	=>\$vector,
	'overlap:i' =>\$overlap,
	'win:i'	    =>\$window,
	'in_vector:i' =>\$insert_vector,
	'conf:s'	=>\$conf,
	'o:s'	=>\$out,
);

my $system=$^O;
#print "$system\n";

#######################################################################################################
#the help
if(!$conf || !$genome ){
		print "\n\n perl $0   -genome Model.gb  -conf Insert_Delete.config.txt  -win 5000 -overlap 100  -vector LeveL2_vector.gb -in_vector [default:-1] -insert_name INSERT1 [default:INSERT_GENOME]  -o output_dir\n\n";
		print "\t\t-genome\tGenome file in genbank format.\n\t\t-conf\tConfig file for Inserted or Deleted location on Genome.\n\t\t-win\tThe window size of Fragment[defualt:5000bp]\n\t\t-overlap\tThe overlap length of two fragment[default:100bp]\n\t\t -vector\tThe Vector file in genbank format.\n\t\t-in_vector\tthe insert site on vector[default:-1,means add a genome to the end of the vector].\n\t\t-insert_name\tThe name of fragment showed on the Vector after inserted.\n\t\t-o\tThe output dir.\n\n";
		exit 0;
		
}
my $path=getcwd;
$genome=abs_path($genome);
$conf=abs_path($conf);
$vector=abs_path($vector);
my %Splice;


############################################################################################################
# output files
#
if(!-d $out){
	`mkdir -p $out`;
}
$out=abs_path($out);
my $out_gb="$out/MarkerDown_genbank.gb";
my $out_fa="$out/After_Insert_Delete.fa";
my $out_ge="$out/Genome_Final.fa";
my $out_log="$out/Genome.processed.log";
my $out_loca="$out/Genome.Split.xls";
my $out_ve="$out/Genome_And_Vector.gb";

if(-f $out_gb){
	`rm -rf $out_gb`;
}

open LOG,">$out_log";
print LOG "Input Genome File:\t$genome\n";
print LOG "Output Genome File:\t$out_gb\n";
print LOG "Afater Insert or Delete Fasta File:\t$out_fa\n";
print LOG "The Log Filer:\t$out_log\n";


if($system=~/MSWin32/){
	$genome=~s/\//\\/g;
	$out=~s/\//\\/g;
}


#########################################################################################################################
#get the length and seq for the chr of genome.
my %genome;
my $ge = Bio::SeqIO->new(-file=>$genome)    or die "couldn't create Bio::SeqIO";
while(my $geseq = $ge->next_seq){
    my @features_ge = $geseq->all_SeqFeatures;# sort features by their primary tags
	my $seqg=$geseq->seq;
	my $chr_id=$geseq->id;
	my $len_genome=length($seqg);
	   $genome{$chr_id}=$seqg;
	   #print("$chr_id\t$len_genome\n");

   }




my $total_in=0;
my $total_de=0;
my %Insert;
my %Mark;




&Conf($conf);

################################################################################################################################################
### read the config file which containing:Insert or Delete location.

sub Conf(){
	my $file_con=shift;
	my %tem_insert;
	my %tem_delete;
    my $insert_len=0;
	my $delete_len=0;
    my %Delete;
	my $final_genome;

	print LOG "\nInsert or Delete region in genome:\n";
	
	open IN,$file_con;
	while(<IN>){
		chomp;
		next if /^$/;
		my @A=(split /\s+/,$_);
		if(/ID/){
			next;
		}else{
			if($A[4]=~/INSERT/){
				$tem_insert{$A[0]}{$A[2]}="$A[1]\t$A[5]";
				$Mark{$A[0]}{'INSERT'}{$A[2]}="$A[1]\t$A[5]";
				my $len_i=length($A[5]);
				   $insert_len+=$len_i;
				   $total_in+=$len_i;
				   print LOG "$A[0]\t$A[1]\t$A[2]\t$len_i\n";
			}else{
				$tem_delete{$A[0]}{$A[2]}{$A[3]}=$A[1];
				$Mark{$A[0]}{'Delete'}{$A[2]}="$A[3]\t$A[1]";
				#print "$A[2]\t$A[1]\t$A[0]\n";
				my $len_d=abs($A[3] - $A[2] + 1);
				  $delete_len+=$len_d;
				  $total_de+=$len_d;
				  print LOG "$A[0]\t$A[1]\t$A[2]-$A[3]\t$len_d\n";
			}
		}
	}
	close IN;

	my $fn=0;
	my $fn_st;
	my $fn_ed;
	my %DelSite;
	
	print LOG "\n\nAfter inserted and delete the genome:\n\n";
	print LOG "Chromosome\tFragmentID\tType\tLength\tStart\tEnd\n";
	open OF,">$out_fa";
	open FA,">$out_ge";
	foreach my $chr(sort keys %tem_delete){
		foreach my $chr_start(sort {$a<=>$b} keys %{$tem_delete{$chr}}){
				foreach my $chr_end(sort {$a<=>$b} keys %{$tem_delete{$chr}{$chr_start}}){
					$fn++;
					$DelSite{$fn}="$chr_start\t$chr_end";
					
					
					if($fn==1){
						if($chr_start>1){
							$fn_st=1;
							$fn_ed=$chr_start - 1;
							#print "Fragment$fn\t$fn_st-$fn_ed\n";
							$Delete{$fn}="$fn_st\t$fn_ed";
						}else{
							$fn_st=$chr_end + 1;
						}
					}else{
						my $dn=$fn - 1;
						if(!$fn_ed){
							$fn_ed=$chr_start - 1;
							#print "Fragment$dn\t$fn_st-$fn_ed\n";
							$Delete{$dn}="$fn_st\t$fn_ed";
						}else{
							if(exists $DelSite{$dn}){
								my @D=(split /\t/,$DelSite{$dn});
								$fn_st=$D[1] + 1;
								$fn_ed=$chr_start - 1;
								#print "Fragment$fn\t$fn_st-$fn_ed\n";
								$Delete{$fn}="$fn_st\t$fn_ed";
								$fn_st=$chr_end + 1;
							}
						}
					}


				}

		}
		my $chr_len=length($genome{$chr});
		$fn=$fn+1;
		#print "Fragment$fn\t$fn_st-End\n";
		$Delete{$fn}="$fn_st\t$chr_len";
		#my $flag=0;
		my $Fran=0;
		foreach my $ochr(sort {$a<=>$b} keys %Delete){
				my @loc=(split /\t/,$Delete{$ochr});
				my $flag=0;
				   
				if(exists $tem_insert{$chr}){
					foreach my $insite(sort {$a<=>$b} keys %{$tem_insert{$chr}}){
							if(@loc>1){
								if($loc[0] < $insite && $loc[1] > $insite){
									   $ochr=$ochr+$Fran;
									my $tem_n=$ochr+1;
									my $tem_n2=$ochr+2;
									my $tem_n2_start=$insite + 1;
									my $len_t1=$insite - $loc[0] + 1;
									my @in=split(/\t/,$tem_insert{$chr}{$insite});
									my $len_t2=length($in[1]);
									my $len_t3=$loc[1] - $tem_n2_start + 1;
									my $fa1=substr($genome{$chr},$loc[0]-1,$len_t1);
									my $fa2=substr($genome{$chr},$tem_n2_start-1,$len_t3);
									   $final_genome.=$fa1;
									   $final_genome.=$in[1];
									   $final_genome.=$fa2;
									print OF ">$chr.Fragment$ochr\tGenome\t$len_t1\t$loc[0]\t$insite\n$fa1\n";
									print OF ">$chr.Fragment$tem_n\t$in[0]\t$len_t2\t$insite\t$tem_n2_start\n$in[1]\n";
									print OF ">$chr.Fragment$tem_n2\tGenome\t$len_t3\t$tem_n2_start\t$loc[1]\n$fa2\n";
									$Splice{$chr}{$ochr}=$fa1;
									$Splice{$chr}{$tem_n}=$in[1];
									$Splice{$chr}{$tem_n2}=$fa2;

									print LOG "$chr\tFragment$ochr\tGenome\t$len_t1\t$loc[0]\t$insite\n";
									print LOG "$chr\tFragment$tem_n\t$in[0]\t$len_t2\t$insite\t$tem_n2_start\n";
									print LOG "$chr\tFragment$tem_n2\tGenome\t$len_t3\t$tem_n2_start\t$loc[1]\n";
									$flag=1;
									$Fran=$Fran + 2;
								}
							}
					}
				}
			if($flag==0){
				my $ochr_new=$ochr+$Fran;
				my $len_t4=$loc[1] - $loc[0] + 1;
				if($ochr_new==1){
					$loc[0]=1;
				}
				my $fa3=substr($genome{$chr},$loc[0]-1,$len_t4);
				   $final_genome.=$fa3;
				   $Splice{$chr}{$ochr_new}=$fa3;
				print OF ">$chr.Fragment$ochr_new\tGenome\t$len_t4\t$Delete{$ochr}\n$fa3\n";
				print LOG  "$chr\tFragment$ochr_new\tGenome\t$len_t4\t$Delete{$ochr}\n";

			}
		}
		my $fin_len=$chr_len + $insert_len - $delete_len;
		print LOG "\n$chr\tTotal Length:$chr_len\tInsert Length:$insert_len\tDelete Length:$delete_len\tFinal Length:$fin_len\n\n";
		print FA ">$chr.Final_genome\t$fin_len\n$final_genome\n";
	}

	close OF;
	close FA;

}

my %Insert_genome;
my $total_genome=0;

foreach my $chrs(sort keys %Mark){
	my $start_i=1;
	my $nu=0;
	my $chr_seq;
	foreach my $typ(keys %{$Mark{$chrs}}){	
		if($typ eq 'INSERT'){
				$nu++;
			foreach my $i(sort {$a<=>$b} keys %{$Mark{$chrs}{$typ}}){
				my $in_seq=(split /\t/,$Mark{$chrs}{$typ}{$i})[1];
				my $end_i=$i;
				my $sub_len=$end_i - $start_i +1;
				my $s_i=$start_i - 1;
				my $sub_i=substr($genome{$chrs},$s_i,$sub_len);
				$chr_seq.=$sub_i.$in_seq;
				#print ">$chrs.Fragmen$nu\t$s_i-$i\n$sub_i\n";
				$nu++;   
				$start_i=$i+1;

			}
		}
	}
	if($start_i>1){
		#$nu++;
		my $si=$start_i -1;
		my $len_si=length($genome{$chrs});
		   $len_si=$len_si - $start_i + 1;
		my $subi=substr($genome{$chrs},$si,$len_si);
		   $chr_seq.=$subi;
		   #print ">$chrs.Fragmen$nu\t$si-$len_si\n$subi\n";
	}
	$Insert_genome{$chrs}=$chr_seq;
	my $lss=length($chr_seq);
	   $total_genome+=$lss;
	#print ">F $lss\n$chr_seq\n";
}

#####################################################################################################################
# out put the genome file after insert to the vector.
my $vn=0;
open OG,">$out_ve";
open INS,"<$vector";
while(<INS>){
	chomp;
	next if /^$/;
	if(/(\d+) bp/){
		my $new_ve=$1+$total_genome;
		$_=~s/(\d+) bp/$new_ve bp/;
	}
	if(/bases 1 to (\d+)/){
		my $new_ve=$1+$total_genome;
		$_=~s/bases 1 to (\d+)/bases 1 to $new_ve/g
	}
	if($vn==0){
		print OG "$_\n";
	}
	if(/FEATURES/){
		#print OG "$_\n";
		$vn++;
	}

}
close INS;

my $ve=Bio::SeqIO->new(-file=>$vector);

while(my $veseq= $ve->next_seq){
		my $vid=$veseq->id;
		my $vseq=$veseq->seq;
		my $len_v=length($vseq);
		my $in_seq;
		my $start_in=0;

		foreach my $inv(sort keys %Insert_genome){
			$in_seq.=$Insert_genome{$inv};
		}
		my $in_seq_len=length($in_seq);

		if($insert_vector > $len_v){
			print "The insert_vector must less than the length of the vector.\n\n";
			exit 0;
		}else{
			if($insert_vector == $len_v || $insert_vector == -1){
				$vseq=$vseq.$in_seq;
				$start_in=$len_v+1;
			}else{
				if($insert_vector==0){
					$vseq=$in_seq.$vseq;
				}else{
					my $sub_v1=substr($vseq,0,$insert_vector);
					my $sub_v2_len=$len_v - $insert_vector + 1;
					my $sub_v2=substr($vseq,$insert_vector-1);
					$vseq=$sub_v1.$in_seq.$sub_v2;
					$start_in=$insert_vector+1;
				}
			}
		}
	   my %ve_fa=&FA($vseq);
	   my $new_ve_len=length($vseq);
	   
	   my @ve_feature=$veseq->all_SeqFeatures;
	   foreach my $vf(@ve_feature){
				my $vle1=$vf->primary_tag;
				my @vlocation=$vf->location;
				my @vtag=$vf->get_all_tags;
				my @vector_vf;
				
				foreach my $vlo(@vlocation){
						my $vstart=$vlo->start;
						my $vend=$vlo->end;
						my $vstrand=$vlo->strand;
						if($vstart > $insert_vector && $insert_vector >= 0 ){
								$vstart=$vstart + $in_seq_len;
						}
						if($vend>$insert_vector && $insert_vector >= 0){
								$vend=$vend+$in_seq_len;
						}

						my $tem_v="  $vle1            1\.\.$new_ve_len";
						if($vle1!~/source/){
								$tem_v="  $vle1            $vstart\.\.$vend";
								if($vstrand eq "-1"){
									$tem_v="  $vle1            complement($vstart\.\.$vend)";

								}
							
						}
						print OG "$tem_v\n";
				}
				foreach my $vv(@vtag){
					my @vval=$vf->get_tag_values($vv);
					my $vtem;
					foreach my $vva(@vval){
							$vtem="               /$vv=\"$vva\"";
							print OG "$vtem\n";
					}
				}

	   }
	    my $end_in=$in_seq_len + $start_in - 1;
		my $tem_in="  misc_feature            $start_in\.\.$end_in";
		my $tem_in_2="               /note=\"INSERT_GENOME\"";
		print OG "$tem_in\n$tem_in_2\n";
		
		print OG "ORIGIN\n";

	    foreach my $fv(sort {$a<=>$b} keys %ve_fa){
			print OG "\t$fv\t$ve_fa{$fv}\n";
		}
		print OG "\/\/\n";
	    close OG;


       

}



###################################################################################################
# out put the mark genbank file

my $new_len=0;

$ge = Bio::SeqIO->new(-file=>$genome)    or die "couldn't create Bio::SeqIO";
while(my $geseq = $ge->next_seq){
	my @features_ge = $geseq->all_SeqFeatures;# sort features by their primary tags
	my $chr_id=$geseq->id;
	my $seqg=$Insert_genome{$chr_id};
	my $len_g=length($seqg);
	my $new_len_g=$len_g;
	   $new_len=$len_g;
	my @ge;

	foreach my $fg (@features_ge){
		my $glevel1=$fg->primary_tag;
		my @location_g=$fg->location;
		my @tagg=$fg->get_all_tags;
		my $nut=0;
		my $frow;
		my @genome;
	
		foreach my $gg(@tagg){
				$nut++;
			my @value=$fg->get_tag_values($gg);
			my $temp;
			foreach my $va(@value){
				    $temp="               /$gg=\"$va\"";
				if($nut==1){
						$temp="               /$gg=\"$va\"";
				}
				push @genome,$temp;
				#print "$temp\n";
          }
		}
	
		$frow=join("\n",@genome);
		foreach my $tag2(@location_g){
			my $start2=$tag2->start;
			my $end2=$tag2->end;
			my $strand=$tag2->strand;
			if($glevel1=~/source/){
					$len_g=$end2;
					$new_len_g=$len_g;
					my $sorr="  $glevel1            1\.\.$new_len\n";
					push @ge,$sorr;

			}
			my $new_start2=&ChangeLoc($chr_id,$start2) ;
			my $new_end2=&ChangeLoc($chr_id,$end2) ;
			my $tem_loc="$new_start2\.\.$new_end2\n";
			if($strand eq "-1"){
				$tem_loc="complement($new_start2\.\.$new_end2)\n";
			}
			#print "$glevel1\t$start2\t$end2\tNew\t$new_start2\t$new_end2\n";
			my $grr1="  $glevel1            $tem_loc" if $glevel1!~/source/;
			my $grr2="$frow\n";
			push @ge,$grr1;
			push @ge,$grr2;
       }
	 }

	
	 if(exists $Mark{$chr_id}){
		 foreach my $types(sort keys %{$Mark{$chr_id}}){
			foreach my $ts(sort {$a<=>$b} keys %{$Mark{$chr_id}{$types}}){
					my @TS=(split /\t/,$Mark{$chr_id}{$types}{$ts});
					if($types eq 'INSERT'){
						my $starts=$ts;
						   $starts=&ChangeLoc($chr_id,$starts) + 1;
						my $end=length($TS[1]) + $starts - 1;
						my $in_row="  misc_feature            $starts..$end\n               \/note=\"$TS[0]\"\n";
						push @ge,$in_row;
				    }else{
						my $startd=$ts;
						   $startd=&ChangeLoc($chr_id,$startd);
						my $endd=$TS[0];
						   $endd=&ChangeLoc($chr_id,$endd);
						my $de_row="  misc_feature            $startd\.\.$endd\n               \/note=\"$TS[1]\"\n";
						push @ge,$de_row;

					}
			}
		 }

	 }
	
	
   if(!$ge[1]){
	   $ge[1]="    ";
   }
   my $grow=join("",@ge);



	my $n_row=0;
	open IN,$genome;
	open O,">>$out_gb";
	while(<IN>){
		chomp;
		if($n_row==0){
			$_=~s/(\d+) bp/$new_len bp/g;
			$_=~s/bases 1 to (\d+)/bases 1 to $new_len/g;
			$_=~s///g;
			print O "$_\n";
		}
		$n_row++ if /FEATURES/;
	}
	close IN;

	print O "$grow\n";
	print O "ORIGIN\n";

	my %fa=&FA($Insert_genome{$chr_id});
	foreach my $fn(sort {$a<=>$b} keys %fa){
		print O "\t$fn\t$fa{$fn}\n";
	}
	print O "\/\/\n";
	close O;
	
}


#########################################################################################################
#output the final splice fragment by window 5000bp and overlap 100bp.

open O,">$out_loca";
print O "Chr\tFragmentID\tLength\tSeq\n";
foreach my $sch(sort keys %Splice){
		my $sch_len=0;
		foreach my $sn(sort {$a<=>$b} keys %{$Splice{$sch}}){
				my $nu=0;
				my $start_b=1;
				my $end_b=1;
				my $spf=length($Splice{$sch}{$sn});
				my @splice=&Splice($Splice{$sch}{$sn});
				my $ID_start=$sch_len + 1;
				my $ID_end=$spf+$ID_start - 1;
				   $sch_len+=$spf;
				
				foreach my $sp(@splice){
						   $nu++;
						my $f_len=length($sp);
						my $tagid="$sch\tFragment$sn\_$ID_start:$ID_end.$nu\t$f_len\t$sp";
						print O "$tagid\n";
				}
		}
}
close O;



sub FA(){
	my $fa=shift;
	my $row=1;
	my $nrow=0;
	my $ttag=1;
	my $len=length($fa);
	my $fn=($len % 10);
	my $maxn=$len/10;
	if($fn>0){
		$maxn=$maxn+1;
	}
	my %FA;
	my %seq;
	for(my $fi=1;$fi<=$maxn;$fi++){
		   $nrow++;
		my $st=($fi-1)*10;
		my $subf=substr($fa,$st,10);
		if($fi==$maxn){
			$subf=substr($fa,$st);
		}
		push @{$FA{$ttag}},$subf;
		if($nrow==6){
			$row++;
			$ttag=($row-1)*60+1;
			$nrow=0;
		}

	}
	foreach my $n(sort {$a<=>$b} keys %FA){
			my $fas=join(" ",@{$FA{$n}});
			   $seq{$n}=$fas
	}
	return(%seq);

}




sub ChangeLoc(){
	my $chrs=shift;
    my $loc=shift;
	my $strand=shift;
	my $new_loc=0;
	my $in_len=0;
	my $in_nums=0;
	my %tem;
    
		if(exists $Mark{$chrs}){
			foreach my $type(sort keys %{$Mark{$chrs}}){
				foreach my $ins(sort {$a<=>$b} keys %{$Mark{$chrs}{$type}}){
						my @tem_in=(split /\t/,$Mark{$chrs}{$type}{$ins});
						if($type eq 'INSERT'){
							   $in_nums++;
							my $len_ins=length($tem_in[1]);
							   $in_len+=$len_ins;
							   if($loc > $ins){
									$new_loc=$loc + $in_len;
									$tem{$in_nums}=$len_ins;
									#print "$loc\t$new_loc\n";
							   }
						}
				}
			}
		}

	if($new_loc==0){
				$new_loc=$loc;
	}
	
	return($new_loc);
}



sub Splice(){
	my $fas=shift;
	my $fas_len=length($fas);
	my $mid_len=int($fas_len/2);
	my $max_len=2*$window;
	my $n_k=$fas_len/$window;
	my $n_k2=int($n_k);
	my $start_i=0;
	if($n_k > $n_k2){
		$n_k=$n_k2 + 1;
	}
	my @fas;
	if($n_k==1){
			push @fas,$fas;
	}else{
			if($n_k==2){
				my $len_s1=$mid_len+$overlap;
				my $len_s2=int($mid_len/2);
				if($overlap > $mid_len){
					$len_s1=$mid_len+$len_s2;
				}
				my $sub1=substr($fas,0,$len_s1);
				my $st2=$len_s1 - $overlap;
				if($overlap > $mid_len){
				   $st2=$len_s1 - $len_s2;
				}
				my $sub2=substr($fas,$st2);
				push @fas,$sub1;
				push @fas,$sub2;
			}else{
				my $nn=$n_k -2;
				for(my $i=1;$i<=$nn;$i++){
					$start_i=($i - 1)*($window - $overlap);
					my $subi=substr($fas,$start_i,$window);
					push @fas,$subi;
					
				}
				my $strl=($n_k - 2)*($window - $overlap);
				my $last_len=length(substr($fas,$strl));
				my $fd=$last_len/$window;
				if($fd>=1.5){
					my $lsub1=substr($fas,$strl,$window);
					my $lsub2=substr($fas,-$window);
					push @fas,$lsub1;
					push @fas,$lsub2;
				}else{
					if($fd<=1){
						my $subl1=substr($fas,$strl);
						push @fas,$subl1;
					}else{
						my $half=int($last_len/2);
						my $half2=int($half/2);
						my $slen=$half+$overlap;
						if($overlap>$half){
							$slen=$half+$half2;
						}
						my $hsub=substr($fas,$strl,$slen);
						my $hsub2=substr($fas,-$slen);
						push @fas,$hsub;
						push @fas,$hsub2;
					}
					
				}

			}
	}
		
	return(@fas);

}









