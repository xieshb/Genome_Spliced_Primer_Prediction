# Genome_Spliced_Primer_Prediction
This Project is used for generating a PCR primers library from a genome which in genbank format and linked the reconsturcted genome to a Vector.

input:
        Level2_Vector_backbone.gb       #the Vector genome in genbank format
        Model_SVEN.gb                           #the Target genome in genbank format

output:
        Gene_info.xls                                   #Gene information
        Genome_Spliced_Fragment.xls             #Genome Fragment Spliced by DELETE masked in misc_feature
        Genome_Spliced_Fragment.fa              #Genome Fragment Fasta Spliced by DELETE masked in misc_feature
        Genome_Spliced_size_segment.fa  #Genome segment cutted by 5000bp as a Kmer
        Genome_Spliced_Linked.fa                #Genome Fragment Fasta Linked after spliced by DELETE
        genome.log                                              #Genome size statistic after spliced
        Primer.input                                    #input config for primer3
        primer3.out                                             #output of primer3
        primer.picked.result.xls                #Top5 primer prediction result of primer3
        primer.final.result.xls                 #Top1 primer prediction result of primer3

Example:
        perl Primer_Prediction_From_GenBank_v1.pl  -vector Level2_Vector_backbone.gb  -genbank Model_SVEN.gb  -width 5000   -overlap 100  -arm 20  -left 25 -right 25  -outdir output_dir

        or
        ./Primer_Prediction_From_GenBank_v1.exe  -vector Level2_Vector_backbone.gb  -genbank Model_SVEN.gb  -width 5000   -overlap 100  -arm 20  -left 25 -right 25  -outdir output_dir

Parameters:
         -genbank        Model_sevn.gb [*input genome in genbank]
         -width          5000 [*input Kmer length]
         -overlap        100 [input overlap length of kmer]
         -arm            20 [the length of homologous arms]
         -left           100 [the left range of genome]
         -right          100 [the right range of genome ]
         -varm           TCGAGTTCATGTGCAGCTCC [homologous arm for vector]
         -vector         Level2_Vector_backbone.gb [vector genome in genbank]
         -primer3        /home/software/bin/primer3_core [path of primer3_core]
         -outdir         Genome_Spliced [*output dir]
         *                               [Required parameter]
         
         
         
         
GUI version:
        Primer_Prediction_From_GenBank_GUI.v1.exe
        Primer_Prediction_From_GenBank_GUI.v1.pl
command version:
        Primer_Prediction_From_GenBank_v1.exe
        Primer_Prediction_From_GenBank_v1.pl
Test Data:
        Model_SVEN.gb
        Level2_Vector_backbone.gb

