# Actual code of network analysis using SMARTTR: simple multi-ensemble atlas registration and statistical testing in R
https://github.com/mjin1812/SMARTTR/tree/9e629dee5f5a49e6c9083bdc30016a4c926e033f

＜data folder＞
 - Please put csv file of your own experiment data on your environment.
 - Rdata are just samples. These experiment data should be created after you excute "1_define_experiment_object.R"
   
＜codes folder＞
Change some paths and parameters to excute. 
 - "1_define_experiment_object.R" should be excuted first.
 - "2_..." is for analysis within each condition such as saline, acute morphine or chronic morphine. This code creates heaetmaps and network plots.
 - "3_..." is for analysis between 2 conditions such as saline_vs_acute_morphine. This code creates parallel plots and volcano plots.
