#!/usr/bin/perl
use strict;
use warnings;
use Tk;
use Cwd;
use Cwd 'abs_path';
use Encode;
use Bio::SeqIO;
use Getopt::Long;


my ($input_mw,$vector,$output_mw,$conf,$insert,$width,$overlap,$varml,$varmr,$left,$right,$min,$max,$opt,$tm,$arm_length,$primer3);
my $path=getcwd;
my $main_primer3="$path/primer3_core.exe";
my $main_bin="$path/Primer_Prediction_From_GenBank_v4.exe";

my $system=$^O;
if($system =~/MSWin32/){
  $main_primer3=~s/\//\\/g;
  $main_bin=~s/\//\\/g;
}


my $mw = MainWindow->new;
   $mw->geometry("1000x920");
   $mw->title("Genome_Spliced_Primer_Prediction");
   #my $frame0 = $mw->Frame->pack(-side=>'top');  
my $code_font=$mw->fontCreate(-family=>'Times',-size=>15,-weight=>'bold');
my $frame0 = $mw->Frame->grid(-row=>1,-column=>5,-rowspan=>2,-columnspan=>30,-pady=>60);  
   $frame0->Button(-text => 'Run',-font=>$code_font,-command =>\&button_para,-height=>3,-width=>4,-background=>'blue')->grid(-row=>2,-column=>6,-sticky => "nsew");
   $frame0->Button(-text => 'Exit',-font=>$code_font,-command =>sub{exit},-height=>3,-width=>4,-background=>'blue')->grid(-row=>2,-column=>7,-sticky => "nsew");


   
###################################################################################################################################################################################
## The input and output file

my $frame1 = $mw->Frame->grid(-row=>3,-column=>3,-rowspan=>5,-columnspan=>20,-pady=>20);
my $button1 = $frame1->Button(-text => "GenBankFile of Genome", -command => \&button1_input, -bg => "grey",-pady=>6)->grid(-row=>3,-column=>3,-sticky => 'w');
my $entry1 = $frame1->Entry(-width=>100)->grid(-row=>3,-column=>4,-pady=>6);
my $button2 = $frame1->Button(-text => "GenBankFile of Vector", -command => \&button2_input, -bg => "grey",-pady=>6)->grid(-row=>4,-column=>3,-sticky => 'w');
my $entry2 = $frame1->Entry(-width=>100)->grid(-row=>4,-column=>4,-pady=>6);
my $button3 = $frame1->Button(-text => "config file", -command => \&button3_input, -bg => "grey",-pady=>6)->grid(-row=>5,-column=>3,-sticky => 'w');
my $entry3 = $frame1->Entry(-width=>100)->grid(-row=>5,-column=>4,-pady=>6);
my $button4 = $frame1->Button(-text => "Outdir", -command => \&button_output, -bg => "grey",-pady=>6)->grid(-row=>6,-column=>3,-sticky => 'w');
my $entry4 = $frame1->Entry(-width=>100)->grid(-row=>6,-column=>4,-pady=>6);

my @help=(" "," Width: The length of genome spliced"," Insert: The site for inserting the genome to a vector,deault add to the end"," "," Overlap: The length of overlap between genome fragment"," "," Length of arm: The length of homologous arm"," "," Length of left range genome: The range length at genome for left primer predicting"," "," Length of right range genome: The range length at genome for right primer predicting"," "," min_primer: The min length of primer"," "," max_primer: The max length of primer"," "," opt_primer: The opted length of primer"," "," Diff_Tm: PRIMER_PAIR_MAX_DIFF_TM"," "," Left arm of Vector: The left arm of Vector,default picked up from Vector genome"," "," Right arm of Vector: The right arm of Vector,default picked up from Vector genome"," "," path of primer3: The path of primer3_core"," ");


#######################################################################################################################################################################################
# The parameters
#
#my $frame2 = $mw->Frame->pack(-side=>'left',-pady=>40,-padx=>10);
my $frame2 = $mw->Frame->grid(-row=>9,-column=>2,-rowspan=>25,-columnspan=>6,-pady=>20);
my $insert_lab = $frame2->Label(-text=>"Insert:",-pady=>6)->grid(-row=>11,-column=>2,-sticky => 'w');   
my $width_lab = $frame2->Label(-text=>"Width:",-pady=>6)->grid(-row=>12,-column=>2,-sticky => 'w');   
my $overlap_lab = $frame2->Label(-text=>"Overlap:",-pady=>6)->grid(-row=>13,-column=>2,-sticky => 'w');   
my $arm_lab = $frame2->Label(-text=>"Length of arm:",-pady=>6)->grid(-row=>14,-column=>2,-sticky => 'w');   
my $left_lab = $frame2->Label(-text=>"Length of left range genome:",-pady=>6)->grid(-row=>15,-column=>2,-sticky => 'w');   
my $right_lab = $frame2->Label(-text=>"Length of right range genome:",-pady=>6)->grid(-row=>16,-column=>2,-sticky => 'w');   
my $min_lab = $frame2->Label(-text=>"min_primer:",-pady=>6)->grid(-row=>17,-column=>2,-sticky => 'w');   
my $max_lab = $frame2->Label(-text=>"max_primer:",-pady=>6)->grid(-row=>18,-column=>2,-sticky => 'w');   
my $opt_lab = $frame2->Label(-text=>"opt_primer:",-pady=>6)->grid(-row=>19,-column=>2,-sticky => 'w');   
my $tm_lab = $frame2->Label(-text=>"Diff_Tm:",-pady=>6)->grid(-row=>20,-column=>2,-sticky => 'w');   
my $varml_lab = $frame2->Label(-text=>"Left arm of Vector:",-pady=>6)->grid(-row=>21,-column=>2,-sticky => 'w');   
my $varmr_lab = $frame2->Label(-text=>"Right arm of Vector:",-pady=>6)->grid(-row=>22,-column=>2,-sticky => 'w');   
my $primer3_lab = $frame2->Label(-text=>"path of primer3:",-pady=>6)->grid(-row=>23,-column=>2,-sticky => 'w');   

my $frame3=$frame2;
#my $frame3 = $mw->Frame->pack(-side=>'left',-pady=>40,-padx=>30);
#my $frame3 = $mw->Frame->grid(-row=>7,-column=>3,-rowspan=>20,-columnspan=>5);
my $insert_ent = $frame2->Entry()->grid(-row=>11,-column=>3,-sticky => 'w');
   $insert_ent->insert('0','END');
   $insert_ent->xview('end');
my $width_ent = $frame3->Entry()->grid(-row=>12,-column=>3,-sticky => 'w');
   $width_ent->insert('0',5000);
   $width_ent->xview('end');
my $overlap_ent = $frame3->Entry()->grid(-row=>13,-column=>3,-sticky => 'w');
   $overlap_ent->insert('0',100);
   $overlap_ent->xview('end');
my $arm_ent = $frame3->Entry()->grid(-row=>14,-column=>3,-sticky => 'w');;
   $arm_ent->insert('0',20);
   $arm_ent->xview('end');
my $left_ent = $frame3->Entry()->grid(-row=>15,-column=>3,-sticky => 'w');
   $left_ent->insert('0',40);
   $left_ent->xview('end');
my $right_ent = $frame3->Entry()->grid(-row=>16,-column=>3,-sticky => 'w');
   $right_ent->insert('0',40);
   $right_ent->xview('end');
my $min_ent = $frame3->Entry()->grid(-row=>17,-column=>3,-sticky => 'w');
   $min_ent->insert('0',18);
   $min_ent->xview('end');
my $max_ent = $frame3->Entry()->grid(-row=>18,-column=>3,-sticky => 'w');
   $max_ent->insert('0',30);
   $max_ent->xview('end');
my $opt_ent = $frame3->Entry()->grid(-row=>19,-column=>3,-sticky => 'w');
   $opt_ent->insert('0',20);
   $opt_ent->xview('end');
my $tm_ent = $frame3->Entry()->grid(-row=>20,-column=>3,-sticky => 'w');
   $tm_ent->insert('0',5);
my $varml_ent = $frame3->Entry()->grid(-row=>21,-column=>3,-sticky => 'w');
my $varmr_ent = $frame3->Entry()->grid(-row=>22,-column=>3,-sticky => 'w');
my $primer3_ent = $frame3->Entry()->grid(-row=>23,-column=>3,-sticky => 'w');
   $primer3_ent->insert('0',$main_primer3);
   $primer3_ent->xview('end');

my $frame4 = $mw->Frame->grid(-row=>12,-column=>20,-pady=>3,-padx=>20,-rowspan=>15);
my $width_h = $frame4->Listbox(-bg=>"white",-height=>26,-width=>85)->pack();
   $width_h->insert('end',@help);


sub button_para	 {
	$width = $width_ent->get();
	$insert = $insert_ent->get();
	$overlap = $overlap_ent->get();
	$arm_length = $arm_ent->get();
	$left = $left_ent->get();
	$right = $right_ent->get();
	$min = $min_ent->get();
	$max = $max_ent->get();
	$opt = $opt_ent->get();
	$tm = $tm_ent->get();
	$varml = $varml_ent->get();
	$varmr = $varmr_ent->get();
	$primer3 = $primer3_ent->get();
	
    if($insert eq 'END'){
		$insert=-1;
	}
	if($min >= $max || $opt >= $max){
		print ("The PRIMER_MAX_SIZE is less than PRIMER_MIN_SIZE or PRIMER_OPT_SIZE\n");
		exit 0;
	}
	if($overlap >= $width){
		print ("The overlap length must be less than the window size\n");
		exit 0;
	}
	$output_mw=$entry4->get();
	$input_mw=$entry1->get();
	$vector=$entry2->get();
	$conf=$entry3->get();
	if (! -d $output_mw){
		`mkdir -p $output_mw`;
	}
	my $input_conf="$output_mw"."/input.conf";
	open(CF,">$input_conf");
	print CF "[P1]\ngenome=$input_mw\nvector=$vector\ngene_conf=$conf\noutdir=$output_mw\n[P2]\nwidth=$width\noverlap=$overlap\narm=$arm_length\nleft=$left\nright=$right\nmin_primer=$min\nmax_primer=$max\nopt_primer=$opt\ndiff_tm=$tm\ninsert_site=$insert\nvector_arm_l=$varml\nvector_arm_r=$varmr\n";
	close CF;

	#system("$main_bin  -genbank $input_mw -vector  $vector -conf $conf -insert_site $insert -width $width -overlap  $overlap  -arm $arm_length -left  $left -right  $right -min_primer $min -max_primer $max  -opt_primer $opt  -diff_tm  $tm  -varmr $varmr  -varml  $varml -primer3 $primer3 -outdir  $output_mw ");
	system("$main_bin  -conf $input_conf ");
	
	print("All Jobs Finished!");
	sleep(5);
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

sub button3_input {
    my $file = $mw->getOpenFile();  #得到输入文件的路径
	   $file = encode("gb2312",$file);
	   if($system=~/MSWin32/){
			$file=~s/\//\\/g;
	   }
	   $conf = $file; 
       $entry3->delete('0','end');
	   $entry3->insert('0',$file);
	   $entry3->xview('end');
	   #print "Your input vector file's "."$file\n";
}	


sub button_output {
	#my $file = $mw->getSaveFile();  #得到输出路径
	my $dir = $mw->chooseDirectory;
	   $dir = encode("gb2312",$dir);
	   $output_mw = abs_path($dir);
	   $output_mw = "$output_mw";
	   if($system=~/MSWin32/){
			$dir=~s/\//\\/g;
			$output_mw=~s/\//\\/g;
	   }
	   if(!-d $output_mw){
		`mkdir  $output_mw`;
	   }
	if (defined $dir and $dir ne '') {
		$entry4->delete('0','end');
		$entry4->insert('0',$dir);
		$entry4->xview('end');
	}
	#print "Your output dir's "."$dir\n";
}


MainLoop;





