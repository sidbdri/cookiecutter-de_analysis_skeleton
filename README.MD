## rMATS
[rmats](http://rnaseq-mats.sourceforge.net/) ([paper](http://www.pnas.org/content/111/51/E5593)) is a computational tool to detect differential alternative splicing events
 from RNA-Seq data. 
 
 The statistical model of MATS calculates the P-value and false discovery rate that the difference in the isoform ratio 
 of a gene between two conditions exceeds a given user-defined threshold. From the RNA-Seq data, MATS can automatically 
 detect and analyze alternative splicing events corresponding to all major types of alternative splicing patterns.
 MATS handles replicate RNA-Seq data from both paired and unpaired study design.
 
 We implemented rMATS into our analytic pipe line. For each comparison (e.g. mutant vs control), rMATS takes the samples 
 from the two conditions in the comparison and workout the statistically enriched (p<0.05) alternative splicing events.
 
 In the rMATS result folder, you will fine a AS_summary_<species>.csv which contains a brief summary of all the comparisons
 made in a particular species. Each comparison has 5 rows of data, corresponding to the 5 different type of splicing events
 currently supported by rMATS: 
  * skipped exons (SE) 
  * alternative 5′ splice sites (A5SS)
  * alternative 3′ splice sites ((A5SS))
  * mutually exclusive exons (MXE)
  * retained introns (RI) 
 
 Most of the column are self explainable by their names. The Up_regulated, Down_regulated, D.E.total shows the number of
 up/down/total AS event. The the Up_regulated_gene, Down_regulated_gene, D.E.total_gene shows the corresponding number 
 in genes in which the above events are found.


 In the rMATS result folder, you will also find a folder for each comparison with a csv file containing the details of 
 each of the five supported AS events. Some of the columns in this csv file contains intermediate result that were used 
 by the rMATS program for the statistical test, thus are more important to us for debugging prosperous. 
 You might want to focus on the PValue and FDR column.
   

 |Comparison	Condition_tested|Total_number_of_samples_data|Base_level_condition|Number_of_samples_in_base_level_condition|Sample_names_in_base_level_condition|Comparison_level_condition|Number_of_samples_in_comparison_level_condition|Sample_names_in_comparison_level_condition|p.adj.cutoff|Up_regulated|Down_regulated|D.E.total|Up_regulated_gene|Down_regulated_gene|D.E.total_gene|
 |:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|:-------------:|
 |mutant_vs_correction_motor_JC.A3SS|condition_2|12|mutant|6|MC1,MC2,MC3,MH1,MH2,MH3|correction|6|ME1,ME2,ME3,MI1,MI2,MI3|0.05|190|129|320|177|120|287|
 

 
