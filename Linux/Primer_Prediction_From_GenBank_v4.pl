################################ 20220401 #####################################################################
#!usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';

my $gb;
my $width;
my $out;
my $insert=-1;
my $overlap_len=100;
my $arm_len=20;
my $left_range=40;
my $right_range=40;
my $vector_gb;
my $varm_left;
my $varm_right;
my $primer_bin;
my $PRIMER_MIN_SIZE=18;
my $PRIMER_MAX_SIZE=30;
my $PRIMER_OPT_SIZE=20;
my $PRIMER_PAIR_MAX_DIFF_TM=5;
my $config;
my $input_conf;


my $system=$^O;
my $start_time=localtime();
my $path=getcwd;
my $pbin="$path/primer3_core.exe";
my $gbbin="$path/output_genbank.pl";
my $mbin="$path/MarkerDown_genbank.pl";

if($system=~/MSWin32/){
	$pbin=~s/\//\\/g;
}

GetOptions(
			"conf:s"	 => \$input_conf,
);


##################################################################################################################################################################
# the help information

if(!$input_conf){
	my $help=<<HH;
		-conf		The input config file
HH
    print "\n $0 \n\n$help\n";

	exit 0;
}

#####################################################################################################################################################################
# read input config
my %config;
my $set="P1";
open(CF,$input_conf) or die "can't open $input_conf:$!";
while(<CF>){
	chomp($_);
	next if $_=~/^$/;
	if($_=~/\[P(\d+)\]/){
		$set="P"."$1";
		#print("$set\n");
	}else{
		my @Row=split("=",$_,2);
		   $config{$set}{$Row[0]}=$Row[1];
		   if ($Row[0] eq "genome"){
				$gb=$Row[1];
			}
		   if($Row[0] eq "vector"){
				$vector_gb=$Row[1];
			}
		   if($Row[0] eq "gene_conf"){
				$config=$Row[1];
			}
		   if($Row[0] eq "outdir"){
				$out=$Row[1];
			}
		   if($Row[0] eq "width"){
				$width=$Row[1];
			}
		   if($Row[0] eq "insert_site"){
				$insert=$Row[1];
			}
		   if($Row[0] eq "overlap"){
				$overlap_len=$Row[1];
			}
		   if($Row[0] eq "arm"){
				$arm_len=$Row[1];
			}
		   if($Row[0] eq "left"){
				$left_range=$Row[1];
			}
		   if($Row[0] eq "right"){
				$right_range=$Row[1];
			}
		   if($Row[0] eq "min_primer"){
				$PRIMER_MIN_SIZE=$Row[1];
			}
		   if($Row[0] eq "max_primer"){
				$PRIMER_MAX_SIZE=$Row[1];
			}
		   if($Row[0] eq "opt_primer"){
				$PRIMER_OPT_SIZE=$Row[1];
			}
		   if($Row[0] eq "diff_tm"){
				$PRIMER_PAIR_MAX_DIFF_TM=$Row[1];
			}
		   if($Row[0] eq "vector_arm_l"){
				$varm_left=$Row[1];
			}
		   if($Row[0] eq "vector_arm_r"){
				$varm_right=$Row[1];
			}


	}

}
close CF;



######################################################################################################################################################################
# define the input and ouput
$gb=abs_path($gb);
$vector_gb=abs_path($vector_gb);


$out=abs_path($out);
$out="$out/Genome_Spliced_Primer_result";
my $primer_dir="$out/Primer_result";
my $frag_dir="$out/Fragments";

if($system=~/MSWin32/){
	$out=~s/\//\\/g;
	$gb=~s/\//\\/g;
	$vector_gb=~s/\//\\/g;
	$primer_dir=~s/\//\\/g;
	$frag_dir=~s/\//\\/g;
}

if(!-d $out){
  `mkdir -p  $out`;
}
if(!-d $primer_dir){
	`mkdir -p $primer_dir`;
}
if(!-d $frag_dir){
	`mkdir -p $frag_dir`;
}


if( !$primer_bin || !-f $primer_bin){
	$primer_bin=$pbin;
}

my $primer_conf="$primer_dir/Primer.input";
my $log="$out/run.log";

if($system=~/MSWin32/){
	$primer_conf=~s/\//\\/g;
	$log=~s/\//\\/g;
}

######################################################################################################################################################
## pick 20bp arm from Vector genome
#
if(-f $vector_gb){
	my $vector = Bio::SeqIO->new(-file=>$vector_gb) or die "couldn't find $vector_gb";
	my $vseq = $vector->next_seq or die "couldn't find a sequence in $vector_gb";
	my $vseqs=$vseq->seq;
	my $vlen=length($vseqs);
	
	if(!$varm_left || !$varm_right){
		if($insert>0 && $insert < $vlen){
			my $start1_v=$insert - $arm_len;
			$varm_left=substr($vseqs,$start1_v,$arm_len);
			$varm_right=substr($vseqs,$insert,$arm_len);
			
		}else{
				$varm_left=substr($vseqs,0,$arm_len);
				$varm_right=substr($vseqs,-$arm_len);
		}
	}

}


#########################################################################################################################################################
# First Insert and Delete the Genome
# Second Splice the Genome with a window default 5000bp and a overlap 100bp
# Third Insert the Fragments into a Vector
#
if(-f $config){
	system("sudo perl $mbin -genome $gb -vector $vector_gb -overlap $overlap_len  -in_vector $insert -win $width -conf $config  -o $frag_dir");
}

my $fragment_f="$frag_dir/Genome.Split.xls";

my $io = Bio::SeqIO->new(-file=>"$gb")    or die "couldn't create Bio::SeqIO";

my $seq_gb = $io->next_seq    or die "couldn't find a sequence in the file";
my $genome_seq=$seq_gb->seq;
my $genome_size=length($genome_seq);


my $num_kmer=0;
my $num_chr=0;
my %Kmer_seq;
my %arm;

open IN,"<$fragment_f";
open O,">$primer_conf";
while(<IN>){
	chomp;
	next if /^$/;
	next if /FragmentID/;
	$num_kmer++;
	my @F=(split /\t/,$_);
	my $F=(split /\./,$F[1])[-1];
	if($F==1){
		$num_chr++;
	}
	my $IDn;
	if($F[1]=~/Fragment(\d+)_(\d+):(\d+)/){
		$IDn="Fragment$1";
	}
	my $tem_left_arm=substr($F[3],0,$arm_len);
	my $tem_right_arm=substr($F[3],-$arm_len);
	   $arm{$IDn}{'left'}=$tem_left_arm;
	   $arm{$IDn}{'right'}=$tem_right_arm;
	my $primer_para=&Primer($F[1],$F[3]);
	print O "$primer_para\n";
	$Kmer_seq{$F[1]}=$F[3];
}
close O;



###########################################################################################################################################################
# Primer Prediction


system("$primer_bin  -format_output  -io_version=4  < $primer_conf >$primer_dir/primer3.out ");

if(-f "$primer_dir/primer3.out"){
	&Parse_Primer("$primer_dir/primer3.out");
}




sub Primer(){
    my @tem=@_;
	my $seq_len=length($tem[1]);
	my $left_arm=substr($tem[1],0,$left_range);
	my $right_arm=substr($tem[1],-$right_range);
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
	#print "$tem[0]\t$seq_len\tSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=0,$start,$end,$start;0,$start2,$end,$start\n";


}


sub GC(){
	my $seq_tag=shift;
	my $len=length($seq_tag);
	my @tag=(split//,$seq_tag);
	my $nu_gc=0;
	foreach my $w(@tag){
		if($w=~/C|G/i){
			$nu_gc++;
		}

	}
	my $gc=sprintf("%0.4f",$nu_gc/$len) if $len;
	   $gc=$gc*100;
	#   $gc=int($gc*100);
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
		print "$seq_tem\t$len_tem\t$gcn\t$atn\t$Tm\n";
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

#my %final_primer;

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
	open O,">$primer_dir/primer.picked.result.xls";
	open F,">$primer_dir/primer.final.result.xls";
	print O "Number\tFragmentID\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tLeft_Final_Length\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tRight_Final_Length\tProduct_len\n";
	print F "Number\tFragmentID\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tLeft_Final_Length\tPrimer_Type\tPrimer_len\tTm1\tGC1\tPrimer_seq\tTm2\tGC2\tPrimer_seq_add_border\tPrimer_seq_after_add_${arm_len}bp_arm\tPrimer_seq_after_add_${arm_len}bp_Vector\tRight_Final_Length\tProduct_len\n";
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
					   $add_varm_r=&Reverse_Complement($varm_right);
					   $add_varm_l=$varm_left;
					
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
				
				#print "$id\tLP:$A[-1]\tarm:$add_arm_l\tvarm:$add_varm_l\n";

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
				  my $left_final=(split /\t/,$row)[-1];
				  my $right_final=(split /\t/,$right_row)[-1];
				  my $left_final_length=length($left_final);
				  my $right_final_length=length($right_final);

				   $orow="$id\t$flag\t$row\t$left_final_length\t$rflag\t$right_row\t$right_final_length";
			       #print "$id\t$flag\t$row\t$rflag\t$right_row\n";
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

my $end_time=localtime();
my $para=<<PP;
	
	Begin at:		$start_time

	The Genome File:		$gb
	The Vector File:		$vector_gb
	The Primer3:			$primer_bin
	The OutputDir:			$out

	Genome_size:					$genome_size
	Kmer_size:					$width
	Kmer From the same Fragment Overlap:		$overlap_len
	Homologous_arm:					$arm_len
	Total_kmers:					$num_kmer
	
	The Left range:		$left_range
	The Right range:	$right_range
	The min primer size:	$PRIMER_MIN_SIZE
	The max primer size:	$PRIMER_MAX_SIZE
	The opt primer size:	$PRIMER_OPT_SIZE
	The diff primer Tm:	$PRIMER_PAIR_MAX_DIFF_TM
    
	Finish at:		$end_time


PP
print "$para\n";
open O,">$log";
print O "$para\n";
close O;
