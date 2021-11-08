#!/usr/bin/perl
use strict;
use warnings;
use Tk;
use Cwd;
use Cwd 'abs_path';
use Encode;
use Bio::SeqIO;
use Getopt::Long;


my ($input_mw,$vector,$output_mw,$width,$overlap,$arm,$left,$right,$min,$max,$opt,$tm,$varm,$primer3);
my $path=getcwd;
my $main_primer3="$path/primer3_core.exe";
my $system=$^O;
if($system =~/MSWin32/){
  $main_primer3=~s/\//\\/g;
}


my $mw = MainWindow->new;
   $mw->geometry("900x800");
   $mw->title("Genome_Spliced_Primer_Prediction");
my $frame0 = $mw->Frame->pack(-side=>'top');  
   $frame0->Button(-text => 'Run',-command =>\&button_para)->pack(-side=>'left');
   $frame0->Button(-text => 'Exit',-command =>sub{exit})->pack(-side=>'right');

my $frame1 = $mw->Frame->pack(-side=>'top',-pady=>2,-padx=>2);
my $button1 = $frame1->Button(-text => "GenBankFile of Genome", -command => \&button1_input, -bg => "grey")->pack();  
my $entry1 = $frame1->Entry(-width=>50)->pack();
my $button2 = $frame1->Button(-text => "GenBankFile of Vector", -command => \&button2_input, -bg => "grey")->pack();
my $entry2 = $frame1->Entry(-width=>50)->pack();
my $button3 = $frame1->Button(-text => "Outdir", -command => \&button_output, -bg => "grey")->pack();
my $entry3 = $frame1->Entry(-width=>50)->pack();

my @help=("Width: The length of genome spliced"," ","Overlap: The length of overlap between genome fragment"," ","Length of arm: The length of homologous arm"," ","Length of left range genome: The range length at genome for left primer predicting"," ","Length of right range genome: The range length at genome for right primer predicting"," ","min_primer: The min length of primer"," ","max_primer: The max length of primer"," ","opt_primer: The opted length of primer"," ","Diff_Tm: PRIMER_PAIR_MAX_DIFF_TM"," ","arm of Vector: The arm of Vector,default picked up from Vector genome"," ","path of primer3: The path of primer3_core"," ");
## 
my $frame2 = $mw->Frame->pack(-side=>'left',-pady=>3,-padx=>3);
my $width_lab = $frame2->Label(-text=>"Width:",-pady=>4)->pack();   
my $overlap_lab = $frame2->Label(-text=>"Overlap:",-pady=>4)->pack();   
my $arm_lab = $frame2->Label(-text=>"Length of arm:",-pady=>4)->pack();   
my $left_lab = $frame2->Label(-text=>"Length of left range genome:",-pady=>4)->pack();   
my $right_lab = $frame2->Label(-text=>"Length of right range genome:",-pady=>4)->pack();   
my $min_lab = $frame2->Label(-text=>"min_primer:",-pady=>4)->pack();   
my $max_lab = $frame2->Label(-text=>"max_primer:",-pady=>4)->pack();   
my $opt_lab = $frame2->Label(-text=>"opt_primer:",-pady=>4)->pack();   
my $tm_lab = $frame2->Label(-text=>"Diff_Tm:",-pady=>4)->pack();   
my $varm_lab = $frame2->Label(-text=>"arm of Vector:",-pady=>4)->pack();   
my $primer3_lab = $frame2->Label(-text=>"path of primer3:",-pady=>4)->pack();   

my $frame3 = $mw->Frame->pack(-side=>'left',-pady=>3,-padx=>3);
my $width_ent = $frame3->Entry()->pack(-pady=>3);
   $width_ent->insert('0',5000);
   $width_ent->xview('end');
my $overlap_ent = $frame3->Entry()->pack(-pady=>3);
   $overlap_ent->insert('0',100);
   $overlap_ent->xview('end');
my $arm_ent = $frame3->Entry()->pack(-pady=>3);
   $arm_ent->insert('0',20);
   $arm_ent->xview('end');
my $left_ent = $frame3->Entry()->pack(-pady=>3);
   $left_ent->insert('0',40);
   $left_ent->xview('end');
my $right_ent = $frame3->Entry()->pack(-pady=>3);
   $right_ent->insert('0',40);
   $right_ent->xview('end');
my $min_ent = $frame3->Entry()->pack(-pady=>3);
   $min_ent->insert('0',18);
   $min_ent->xview('end');
my $max_ent = $frame3->Entry()->pack(-pady=>3);
   $max_ent->insert('0',30);
   $max_ent->xview('end');
my $opt_ent = $frame3->Entry()->pack(-pady=>3);
   $opt_ent->insert('0',20);
   $opt_ent->xview('end');
my $tm_ent = $frame3->Entry()->pack(-pady=>3);
   $tm_ent->insert('0',5);
my $varm_ent = $frame3->Entry()->pack(-pady=>3);
my $primer3_ent = $frame3->Entry()->pack(-pady=>3);
   $primer3_ent->insert('0',$main_primer3);
   $primer3_ent->xview('end');

my $frame4 = $mw->Frame->pack(-side=>'left',-pady=>2,-padx=>2);
my $width_h = $frame4->Listbox(-bg=>"white",-height=>17,-width=>90)->pack();
   $width_h->insert('end',@help);


sub button_para	 {
	$width = $width_ent->get();
	$overlap = $overlap_ent->get();
	$arm = $arm_ent->get();
	$left = $left_ent->get();
	$right = $right_ent->get();
	$min = $min_ent->get();
	$max = $max_ent->get();
	$opt = $opt_ent->get();
	$tm = $tm_ent->get();
	$varm = $varm_ent->get();
	$primer3 = $primer3_ent->get();
	
	if($min >= $max || $opt >= $max){
		print ("The PRIMER_MAX_SIZE is less than PRIMER_MIN_SIZE or PRIMER_OPT_SIZE\n");
	}else{
		#system("perl $main_bin   -genbank $input_mw -width $width -overlap $overlap -arm $arm -left $left -right $right -min_primer $min -max_primer $max -opt_primer $opt -diff_tm $tm -varm $varm -vector $vector -primer3 $primer3 -outdir $output_mw ");
		&main_bin($input_mw,$width,$output_mw,$overlap,$arm,$left,$right,$vector,$varm,$primer3,$min,$max,$opt,$tm);
	}
	
	exit(0);
}	


sub button1_input {
    my $file = $mw->getOpenFile();  #得到输入文件的路径
	   $file = encode("gb2312",$file);
	   if($system=~/MSWin32/){
			$file=~s/\//\\/g;
	   }
	   $input_mw = $file;  
	   $entry1->delete('0','end');
	   $entry1->insert('0',$file);
	   $entry1->xview('end');
	   #print "Your input genome file's "."$file\n";
	
}
sub button2_input {
    my $file = $mw->getOpenFile();  #得到输入文件的路径
	   $file = encode("gb2312",$file);
	   if($system=~/MSWin32/){
			$file=~s/\//\\/g;
	   }
	   $vector = $file; 
       $entry2->delete('0','end');
	   $entry2->insert('0',$file);
	   $entry2->xview('end');
	   #print "Your input vector file's "."$file\n";
}	

sub button_output {
	#my $file = $mw->getSaveFile();  #得到输出路径
	my $dir = $mw->chooseDirectory;
	   $dir = encode("gb2312",$dir);
	   $output_mw = abs_path($dir);
	   $output_mw = "$output_mw/Genome_Spliced_Primer_result";
	   if($system=~/MSWin32/){
			$dir=~s/\//\\/g;
			$output_mw=~s/\//\\/g;
	   }
	   if(!-d $output_mw){
		`mkdir  $output_mw`;
	   }
	if (defined $dir and $dir ne '') {
		$entry3->delete('0','end');
		$entry3->insert('0',$dir);
		$entry3->xview('end');
	}
	#print "Your output dir's "."$dir\n";
}


MainLoop;


my %Splice;
my %gene;
my %arm;
my %Splice_Genome;
my $left_range;
my $right_range;
my $overlap_len;
my $arm_len;
my $PRIMER_MIN_SIZE;
my $PRIMER_MAX_SIZE;
my $PRIMER_OPT_SIZE;
my $PRIMER_PAIR_MAX_DIFF_TM;
my $out;
my $gb;
my $vector_gb;
my $primer_bin;
my %Kmer_seq;
	
my $start=0;
my $end=0;
my $sequ;
my $len=0;
my @order_Frag;
my $genome_size=0;
my $right_arm;
my $left_arm;
my $num_chr=0;

sub main_bin{
	   $gb=shift;
	   $width=shift;
	   $out=shift;
	   $overlap_len=shift;
	   $arm_len=shift;
	   $left_range=shift;
	   $right_range=shift;
	   $vector_gb=shift;
	   $varm=shift;
	   $primer_bin=shift;
	   $PRIMER_MIN_SIZE=shift;
	   $PRIMER_MAX_SIZE=shift;
	   $PRIMER_OPT_SIZE=shift;
	   $PRIMER_PAIR_MAX_DIFF_TM=shift;
	my $start_time=localtime();
    
	if(!-d $out){
		`mkdir $out`;
	}
	if($system=~/MSWin32/){
		$out=~s/\//\\/g;
	}

	if(!-f $primer_bin){
      $primer_bin=$main_primer3;
	}


	my $gene_info="$out/Gene_info.xls";
	my $genome_new="$out/Genome_Spliced_Linked.fa";
	my $seq_spliced1="$out/Genome_Spliced_Fragment.xls";
	my $seq_spliced2="$out/Genome_Spliced_Fragment.fa";
	my $seq_spliced3="$out/Genome_Spliced_size_segment.fa";
	my $primer_conf="$out/Primer.input";
	my $log="$out/genome.log";

	if($system=~/MSWin32/){
		$gene_info=~s/\//\\/g;
	    $genome_new=~s/\//\\/g;
		$seq_spliced1=~s/\//\\/g;
		$seq_spliced2=~s/\//\\/g;
		$seq_spliced3=~s/\//\\/g;
		$primer_conf=~s/\//\\/g;
		$log=~s/\//\\/g;
	}


	my $misc_varm;
	my $misc_source;
    ## pick 20bp in misc_feature from Vector genome
	if(-f $vector_gb){
	    my $vector = Bio::SeqIO->new(-file=>$vector_gb) or die "couldn't find $vector_gb";
	    my $vseq = $vector->next_seq or die "couldn't find a sequence in $vector_gb";
	    my @features_v = $vseq->all_SeqFeatures;
	    for my $vf(@features_v){
	        my $vlevel1=$vf->primary_tag;
	        my @location_v=$vf->location;
	        foreach my $vl(@location_v){
	                my $start_v=$vl->start;
	                my $end_v=$vl->end;
	                my $strand_v=$vl->strand;
	                my $seq_v=$vseq->subseq($start_v,$end_v);
	                my $varm_tem=substr($seq_v,0,$arm_len);
	                if($vlevel1=~/source/){
                          $misc_source=$varm_tem;
	                }
			}
	  }
	}
    
	if(!$varm){
		if($misc_varm){
			$varm=$misc_varm;
		}else{
			$varm=$misc_source;
		}
	}

    my $io = Bio::SeqIO->new(-file=>$gb)    or die "couldn't create Bio::SeqIO";
	my $seq = $io->next_seq    or die "couldn't find a sequence in the file";
	my @features = $seq->all_SeqFeatures;# sort features by their primary tags
	my %sorted_features;

    
	for my $f (@features) {
		my @tag;
	    my @name;
		my @location;
		my $level1 = $f->primary_tag;
		   @tag = $f->get_all_tags;
		   @location = $f->location;

		if($level1 eq "source"){
			next;
		}else{
			if($level1 eq "misc_feature" || $level1 eq "CDS" || $level1 eq "gene"){
				@name = $f->get_tag_values('note');
				foreach my $l(@location){
					my $start = $l->start;
					my $end   = $l->end;
					my $strand= $l->strand;
					my $exon_seq = $seq->subseq($start,$end);
					my $len = abs($end - $start + 1);
					my $row = "$start\t$end\t$len";
	                if($name[0] eq "DELETE"){
						$Splice{$start}{$end}=$strand;
					}else{
						$gene{$row}{$name[0]}{$level1}=$strand;
                    }
					 #print()"$name[]0]\t$start\t$end\t$len\t$strand\t$level1\n");
				}
			}else{
				next;
			}
		}
	}
	
	my %order;
	foreach my $location (keys %gene){
		my @tem=(split /\t/,$location);
	    foreach my $ge(keys %{$gene{$location}}){
	          my @feature_type=keys %{$gene{$location}{$ge}};
	          my $feature_type=join(";",sort(@feature_type));
	          my $strand = $gene{$location}{$ge}{$feature_type[0]};
	          my $delete = "KEEP";
	          if( exists $Splice{$tem[0]}){
	              $delete = "DELETE";
	          }
	          my @temp=(split /\t/,$location);
	             $order{$temp[0]}{$temp[1]}{$strand}="$ge\t$location\t$strand\t$delete\t$feature_type";
	             #print "$ge\t$location\t$strand\t$delete\t$feature_type\n";
	    }
	}
	open O,">$gene_info";
	print O "Gene\tStart\tEnd\tLength\tStrand\tNote\tFeatureType\n";
	foreach my $k1 (sort {$a<=>$b} keys %order){
		foreach my $k2 (keys %{$order{$k1}}){
		    foreach my $str(keys %{$order{$k1}{$k2}}){
	                print O ("$order{$k1}{$k2}{$str}\n");
	        }
	   }
	}
	close O;


	### Splice Genome by DELETE or SCAR
	
	my $nu=0;
	
	
	
	open O,">$seq_spliced1";
	open O1,">$seq_spliced2";
	open O2,">$genome_new";
	print O "Fragment_ID\tStart\tEnd\tOrign_Start\tOrign_End\tLength\n";

	foreach my $st(sort {$a <=> $b } keys %Splice ){
	       $nu++;
		my $tem_seq;
		foreach my $ed(keys %{$Splice{$st}}){
		     #print "$nu\t$st\t$ed\t$Splice{$st}{$ed}\n";
		      if($nu==1){
				 $start=$st-1;
				 $end=$ed + 1;
				 $tem_seq=$seq->subseq(1,$start);
				 $sequ=$tem_seq;
				 $len=length($tem_seq);
				 $left_arm=substr($tem_seq,0,$arm_len);
				 $right_arm=substr($tem_seq,-$arm_len);
				 $arm{"Fragment$nu"}{'left'}=$left_arm;
			     $arm{"Fragment$nu"}{'right'}=$right_arm;
				 #print "Fragment$nu\_1:$start\t$left_arm\t$right_arm\n";
				 print O "Fragment$nu\_1:$start\t1\t$start\t1\t$start\t$len\n";
				 print O1 ">Fragment$nu\_1:$start\tLength:$len\n$tem_seq\n";
				 $Splice_Genome{"Fragment$nu\_1:$start"}=$tem_seq;
				 push @order_Frag,"Fragment$nu\_1:$start";
			}else{
				 $start=$end;
				 $end=$st - 1;
				 $tem_seq=$seq->subseq($start,$end);
				 $sequ.=$tem_seq;
				 my $tem_len=length($tem_seq);
				 my $new_start=$len + 1;
				 my $new_end=$new_start + $tem_len - 1;
				    $len=$len+$tem_len;
				    $left_arm=substr($tem_seq,0,$arm_len);
				    $right_arm=substr($tem_seq,-$arm_len);
				if(!exists $arm{"Fragment$nu"}{'left'}){
			        $arm{"Fragment$nu"}{'left'}=$left_arm;
				}
				$arm{"Fragment$nu"}{'right'}=$right_arm;
			    $Splice_Genome{"Fragment$nu\_$start:$end"}=$tem_seq;
				#print "Fragment$nu\_$start:$end\t$left_arm\t$right_arm\n";
			   print O "Fragment$nu\_$start:$end\t$new_start\t$new_end\t$start\t$end\t$tem_len\n";
			   print O1 ">Fragment$nu\_$start:$end\tLength:$tem_len\n$tem_seq\n";
			   push @order_Frag,"Fragment$nu\_$start:$end";
			   $end = $ed + 1;
			}
		}
	}

	$nu++;
	$num_chr=$nu;
	$start=$end;
	$end = $seq->length;
	my $tem_seq=$seq->subseq($start,$end);
	   $sequ.=$tem_seq;
	my $tem_len=length($tem_seq);
	my $new_start=$len + 1;
	my $new_end=$new_start + $tem_len -1;
	   $left_arm=substr($tem_seq,0,$arm_len);
	   $right_arm=substr($tem_seq,-$arm_len);
	if(!exists $arm{"Fragment$nu"}{'left'}){
	   $arm{"Fragment$nu"}{'left'}=$left_arm;
	}
	$arm{"Fragment$nu"}{'right'}=$right_arm;

	#print "Fragment$nu\_$start:$end\t$left_arm\t$right_arm\n";
	print O "Fragment$nu\_$start:$end\t$new_start\t$new_end\t$start\t$end\t$tem_len\n";
	print O1 ">Fragment$nu\_$start:$end\tLength:$tem_len\n$tem_seq\n";
	$Splice_Genome{"Fragment$nu\_$start:$end"}=$tem_seq;
	push @order_Frag,"Fragment$nu\_$start:$end";
	$len=length($sequ);
	$genome_size=$len;
	my $name=$seq->keywords;
	   $name=~s/\s+/_/g;
	print O2 ">$name\_Spliced_Genome  Length:$len\n$sequ\n";
	close O;
	close O1;
	close O2;
	my $num_kmer=0;
	### splice genome as 5kb segment.
	open O,">$seq_spliced3";
	open PR,">$primer_conf";

	foreach my $fr(@order_Frag){
        my %Kmer=&Kmer($fr,$Splice_Genome{$fr},$width,$overlap_len);
        foreach my $k(keys %Kmer){
	        foreach my $kn(sort keys %{$Kmer{$k}}){
	               $num_kmer++;
				my $pid="$k\.$kn";
                my $len_k=length($Kmer{$k}{$kn});
                my $prim;
                if($kn==0){
						print O ">$k\t$len_k\n$Kmer{$k}{$kn}\n";
						$prim=&Primer($k,$Kmer{$k}{$kn});
						$pid=$k;
				}else{
					    print O ">$k.$kn\t$len_k\n$Kmer{$k}{$kn}\n";
					    $prim=&Primer($pid,$Kmer{$k}{$kn});
               }
				$Kmer_seq{$pid}=$Kmer{$k}{$kn};
				print PR "$prim\n";
			}
		}
	}
	close O;
	close PR;

	system("$primer_bin  -format_output  -io_version=4  < $primer_conf >$out/primer3.out ");

	if(-f "$out/primer3.out"){
		&Parse_Primer("$out/primer3.out");
	}

	my $end_time=localtime();
	my $para=<<PP;

Begin at:       $start_time

The Genome File:        $gb
The Vector File:        $vector_gb
The Primer3:            $primer_bin
The OutputDir:          $out

Genome_size:                    $genome_size
Kmer_size:                  $width
Kmer From the same Fragment Overlap:        $overlap_len
Homologous_arm:                 $arm_len
Total_kmers:                    $num_kmer

The Left range:     $left_range
The Right range:    $right_range
The min primer size:    $PRIMER_MIN_SIZE
The max primer size:    $PRIMER_MAX_SIZE
The opt primer size:    $PRIMER_OPT_SIZE
The diff primer Tm: $PRIMER_PAIR_MAX_DIFF_TM

Finish at:      $end_time

PP

	print "$para\n";
	open O,">$log";
	print O "$para\n";
	close O;




}


sub Kmer(){
      my($ID,$sequence,$width,$overlap)=@_;
      my $step=$width - $overlap;
	  my $seq_len=length($sequence);
	  my %hash;

	  #print "$ID\t$width\t$overlap\t$seq_len\n";
	  chomp($sequence);
	  if($seq_len <= $width){
		  #return($ID,$sequence);
		 my $tem=0;
		    $hash{$ID}{$tem}=$sequence;
		    
	  }else{
			my $num_kmer=int(($seq_len - $overlap)/$step);
			my $last_len=$step * $num_kmer - $overlap * ($num_kmer -1);
			my $last_star=$last_len + ($num_kmer - 1)*$overlap;
			   #print "$seq_len\t$overlap\t$step\t$num_kmer\t$last_len\t$last_star\n";
			   #my $last_seq=substr($sequence,-$width);
			my $last_seq=substr($sequence,$last_star,);
			my $l_id=$num_kmer + 1;
			my $last_id="$ID\.$l_id";
			   $hash{$ID}{$l_id}=$last_seq;
			for(my $i=1;$i<$l_id;$i++){
				#print "$ID\t$i\t$num_kmer\t$l_id\n";
				 my $pos=($i - 1)*$step;
				 my $subseq=substr($sequence,$pos,$width);
				 my $id_tag="$ID\.$i";
				    $hash{$ID}{$i}=$subseq;
				 my $ktem=$i+1;
				 if($ktem==$l_id){
				    my $temID="$ID\.$ktem";
				    my $klen=length($hash{$ID}{$ktem});
				    if($klen<1000){
						my $laslen=$klen+$width;
						my $half=int($laslen/2);
					       $half=$half+$overlap;
					    my $seq1=substr($sequence,$pos,$half);
					       $hash{$ID}{$i}=$seq1;
						   $pos=$pos+$half-$overlap;
						my $seq2=substr($sequence,$pos,);
					       $hash{$ID}{$l_id}=$seq2;
					}
				}
			}
		}
		return(%hash);
}



sub Primer(){
    my @tem=@_;
	my $seq_len=length($tem[1]);
    my $left_arm=substr($tem[1],0,$left_range);
	my $right_arm=substr($tem[1],-$right_range);
	#print "Left arm:\t$left_arm\n";

    my $left_gc=&GC($left_arm);
    my $right_gc=&GC($right_arm);
    my $avg=int(($left_gc + $right_gc)/2);
	my $min_gc=int($avg/2);
	my $max_gc=$avg + 10;



	#my $start=$overlap_len ;
	my $start=$left_range ;
	if($seq_len < $start){
		  $start=30;
    }
	my $tem=int($start/2);
	#my $end=$seq_len - $start - $tem;
	#my $end=$seq_len - $start - ($right_range * 2);
	my $end=$seq_len - $left_range - $right_range;
	
	my $primer_conf=<<PP;
		
SEQUENCE_INCLUDED_REGION=0,$seq_len
SEQUENCE_TARGET=$start,$end
SEQUENCE_FORCE_LEFT_START=0
PRIMER_NUM_RETURN=5
PRIMER_SALT_MONOVALENT=20
PRIMER_MIN_GC=0
PRIMER_MAX_GC=100
PRIMER_MIN_TM=35
PRIMER_MAX_TM=90
PRIMER_OPT_TM=65
PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE
PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE
PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE
PRIMER_PAIR_MAX_DIFF_TM=$PRIMER_PAIR_MAX_DIFF_TM
PRIMER_PRODUCT_SIZE_RANGE=100-$seq_len
PRIMER_PICK_ANYWAY=1
=
PP
	chomp($primer_conf);
	my $seq=$tem[1];
	chomp($seq);

    my $par="SEQUENCE_ID=$tem[0]\nSEQUENCE_TEMPLATE=$seq$primer_conf";

	    return($par);
	    #print "$tem[]0]\t$seq_len\tSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,$start,$end,$start;0,$start2,$end,$start\n";
		
		
		
}


sub GC(){
    my $seq_tag=shift;
	my $len=length($seq_tag);
	my @tag=(split//,$seq_tag);
	my $nu_gc=0;
	my $gc=0;
	foreach my $w(@tag){
	        if($w=~/C|G/i){
			    $nu_gc++;
			}
	}    
	   $gc=sprintf("%0.4f",$nu_gc/$len) if $len;
	   $gc=$gc*100;
	   $gc=int($gc*100);
	return($gc);
}
											
sub Tm(){
	 my $seq_tem=shift;
	 my $len_tem=length($seq_tem);
	 my $gcp=&GC($seq_tem);
	 my $gcn=$len_tem * $gcp/100;
	 my $atn=$len_tem-$gcn;
	 my $Tm=0;
	 #PCR添加剂、辅助溶剂和修饰核苷酸的使用会使得Tm降低5~6度。
	 if($len_tem < 24){
			$Tm=4*$gcn + 2*$atn - 5;
	 }else{
		    $Tm=85 + 16.6*log(0.3) + 0.4*$gcp - (675/$len_tem);
	 }
	 $Tm=sprintf("%.2f",$Tm);
	if($Tm > 100){
			#print "$seq_tem\t$len_tem\t$gcn\t$atn\t$Tm\n";
	}
	return($Tm);
			
}
											   
sub Reverse_Complement(){
    my $seqs=shift;
       $seqs=reverse $seqs ;
	   $seqs=~tr/tcgaTCGA/AGCTagct/;
	   $seqs=uc($seqs);
	return($seqs);
}




sub Parse_Primer(){
    my $primer=shift;

	my $id;
	my $number;
	my $max_k=0;
	my $orow;
	my %hash;
	my $nu=0;
	my $add_arm_l;
	my $add_arm_r;
	my $add_varm_l;
	my $add_varm_r;
	my $suffix;
	my $seq;
    
	open IN,$primer;
	open O,">$out/primer.picked.result.xls";
	open F,">$out/primer.final.result.xls";
    print O "Number\tFragmentID\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tProduct_len\n";
    print F "Number\tFragmentID\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tProduct_len\n";
    while(<IN>){
		next if /^$/;
		chomp();

		 my @A=(split /\s+/,$_);
	    if(/PRIMER PICKING RESULTS FOR (.*)/){
		                $id=$1;
						$id=~s///g;
						$nu=1;
						$seq=$Kmer_seq{$id};
						if($id=~/Fragment(\d+)_(\d+):(\d+)/){
						       $number=$1;
						    my $F_len=abs($3 - $2);
							   $suffix=(split /\./,$id)[-1];
							if($suffix=~/Fragment/){
							   $suffix=0;
							}
							my $befor=$number-1;
							my $after=$number+1;
							   $add_arm_l=$arm{"Fragment$befor"}{'right'} if exists $arm{"Fragment$befor"}{'right'};
							   $add_arm_r=$arm{"Fragment$after"}{'left'} if exists $arm{"Fragment$after"}{'left'};
							   $add_arm_r=&Reverse_Complement($add_arm_r);
							   $add_varm_l=&Reverse_Complement($varm);
							   $add_varm_r=$varm;
                           my  $step =$width - $overlap_len;
	                           $max_k=int($F_len/$step);
		                       $max_k=$max_k + 1;
						}
		}

		if(/LEFT PRIMER/){
	        my $flag=$A[-9];
			my $left_primer_add_border="";
			   $left_primer_add_border=substr($seq,0,$A[-7]+$A[-6]) if $seq;
			my $left_Tm=$A[-5];
		    my $left_gc=$A[-4];
			if($A[-1] ne $left_primer_add_border){
			     $left_gc=&GC($left_primer_add_border);
				 $left_Tm=&Tm($left_primer_add_border);
			}
            my $row="$A[-6]\t$A[-5]\t$A[-4]\t$A[-1]\t$left_Tm\t$left_gc\t$left_primer_add_border";
            if($number==1 && ($suffix==0 || $suffix==1)){
			       $row="$row\t$left_primer_add_border\t$add_varm_l$left_primer_add_border";
			}
			if($number==1 && $suffix > 1){
					$row="$row\t$left_primer_add_border\t$left_primer_add_border";
			}
			if($number==$num_chr && ($suffix==0 || $suffix==1)){
					$row="$row\t$add_arm_l$left_primer_add_border\t$add_arm_l$left_primer_add_border";
			}
			if($number==$num_chr && $suffix > 1){
				    $row="$row\t$left_primer_add_border\t$left_primer_add_border";
			}
			if($number > 1 && $number < $num_chr){
				   if($suffix==0 || $suffix==1){
						 $row="$row\t$add_arm_l$left_primer_add_border\t$add_arm_l$left_primer_add_border";
				   }else{
						 $row="$row\t$left_primer_add_border\t$left_primer_add_border";
				   }	
			 }

			#print "$id\tLP:$A[]-1]\tarm:$add_arm_l\tvarm:$add_varm_l\n";
			
			my $right=<IN>;
			   $right=~s/\r//g;
			my @tem=(split /\s+/,$right);
			my $rflag=$tem[-9];
			my $right_primer_add_border="";
			   $right_primer_add_border=substr($seq,$tem[-7]+1,) if $seq;
			   $right_primer_add_border=&Reverse_Complement($right_primer_add_border);
			   $right_primer_add_border="$right_primer_add_border$tem[-1]";
			my $right_Tm=$tem[-5];
		    my $right_gc=$tem[-4];
		   if($tem[-1] ne $right_primer_add_border){
			       $right_gc=&GC($right_primer_add_border);
				   $right_Tm=&Tm($right_primer_add_border);
			}

			my $right_row="$tem[-6]\t$tem[-5]\t$tem[-4]\t$tem[-1]\t$right_Tm\t$right_gc\t$right_primer_add_border";
            if($number==1){
				if($max_k==1 || $suffix==$max_k){
					    $right_row="$right_row\t$add_arm_r$right_primer_add_border\t$add_arm_r$right_primer_add_border";
				}else{
					    $right_row="$right_row\t$right_primer_add_border\t$right_primer_add_border";
				}
			}
			if($number==$num_chr){
				if($max_k==1 || $suffix==$max_k){
					    $right_row="$right_row\t$right_primer_add_border\t$add_varm_r$right_primer_add_border";
				}else{
					    $right_row="$right_row\t$right_primer_add_border\t$right_primer_add_border";
				}
			}
			if($number > 1 && $number < $num_chr){
				if($max_k==1){
						$right_row="$right_row\t$add_arm_r$right_primer_add_border\t$add_arm_r$right_primer_add_border";
				}else{
				    if($suffix==$max_k){
						$right_row="$right_row\t$add_arm_r$right_primer_add_border\t$add_arm_r$right_primer_add_border";
					}else{
						$right_row="$right_row\t$right_primer_add_border\t$right_primer_add_border";
					}
				}
			}
			
			$orow="$id\t$flag\t$row\t$rflag\t$right_row";
			
		}

		if(/PRODUCT SIZE:\s+(\d+),/){
		    my $product_len=$1;
		    print O "$nu\t$orow\t$product_len\n";
		    if($nu==1){
				print F "$nu\t$orow\t$product_len\n";
			}
			$nu++;
		}
	}
 close IN;
 close O;
 close F;
}














