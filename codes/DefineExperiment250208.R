#Follow this website:https://mjin1812.github.io/SMARTTR/articles/Part5.ImportingExternalDatasets.html
#For details:https://github.com/mjin1812/SMARTTR/tree/9e629dee5f5a49e6c9083bdc30016a4c926e033f/.github
# Install necessary packages
install.packages("SMARTTR")
install.packages("magrittr")
install.packages("dplyr")
install.packages("/private/var/folders/lx/yxlnvzs96291fgg1fvcnkqf80000gn/T/RtmpPYRBfK/downloaded_packages/SMARTTR_1.0.1.tar.gz", repos = NULL, type = "source")
library(SMARTTR)
library(magrittr)
library(dplyr)

#create experiment object as opioid
opioid <- experiment(experiment_name = "opioid",
                         channels = "cfos",       # If you have more than one channel to import, set this to a character vector, e.g. c("cfos", "eyfp")
                         output_path = "/Users/saccyann/Documents/Sakura_networkanalysis") #Set this to a path location where you want your figures/analysis output to save, e.g. "P:\\DENNYLABV\\Michelle_Jin\\experiment\\folder"
# output_path = tempdir()) #Set this to a path location where you want your figures/analysis output to save, e.g. "P:\\DENNYLABV\\Michelle_Jin\\experiment\\folder"
print(opioid)

#import externally mapped datasets
opioid <- import_mapped_datasets(opioid, 
                                     normalized_count_paths = "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/data/cleaned_countdata_wobg_woco.csv",
                                     show_col_types = FALSE)
print(opioid)
head(SMARTTR::ontology.unified)

ontology <- "unified"   # Set to "allen" if you are using the allen ontology
opioid <- check_ontology_coding(opioid, ontology = ontology)

# Filters to chosen base parent regions and all child subregions
load("/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/data/rejected_acronyms.RData")
opioid <- filter_regions(opioid, 
                             base_regions =  dataset,
                             ontology = ontology)

opioid <- exclude_redundant_regions(opioid, ontology = ontology)

# simplify the acronyms by the keywords and recalculate the normalized counts
simplify_keywords <-c("interfascicular part")
opioid <- simplify_cell_count(opioid, ontology = ontology, simplify_keywords = simplify_keywords, dont_fold = "")
exclusion_acronyms <- c("drt")
opioid <- exclude_by_acronym(opioid, acronyms = exclusion_acronyms, ontology = ontology)
keywords_exclude <- c("nerve", "tract", "pineal gland", "Area subpostrema")
opioid <-  exclude_by_keyword(opioid, keywords = keywords_exclude)

#quality check
opioid$combined_normalized_counts$cfos <- opioid$combined_normalized_counts$cfos %>% dplyr::filter(counts >= thresh)
opioid <- find_outlier_counts(opioid, by = c("group"), n_sd = 2, remove = TRUE, log = TRUE)
opioid <- enough_mice_per_group(opioid, by = c("group"), min_n = 4, remove = TRUE, log = TRUE)

save(opioid, file = "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/opioid_labdata_check.RData")
save_experiment(opioid, timestamp = TRUE)
print("success")
