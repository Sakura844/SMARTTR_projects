# Actual code of network analysis using SMARTTR: simple multi-ensemble atlas registration and statistical testing in R
https://github.com/mjin1812/SMARTTR/tree/9e629dee5f5a49e6c9083bdc30016a4c926e033f

＜data folder＞
 - "cellcounts_per_acronyms.csv"
   Format is specified in following webpage.
   https://mjin1812.github.io/SMARTTR/articles/Part5.ImportingExternalDatasets.html

   Please put csv file of your own experiment data on your environment and change the path in "1_define_experiment_object.R".
   There is a function "formatCsv" which generate formatted csv inside of the "1_define_experiment_object.R".
   
 - "acronyms_to_use.RData"
   This file is an example of how to specify the acronims that you use for this analysis.
   
＜codes folder＞
Change some paths and parameters to excute. 
 - "1_define_experiment_object.R" should be excuted first.
 - "2_within_single_condition.R" is for analysis within each condition such as saline, acute morphine or chronic morphine. This code creates heaetmaps and network plots.
 - "3_corr_between_several_conditions.R" is for analysis between 2 conditions such as saline_vs_acute_morphine. This code creates parallel plots and volcano plots.
