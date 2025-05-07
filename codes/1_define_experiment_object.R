# Install necessary packages
install.packages("SMARTTR")
install.packages("magrittr")
install.packages("dplyr")
library(SMARTTR)
library(magrittr)
library(dplyr)

formatCsv <- function(original_path = NULL, #Original csv file path
                       set_path = NULL, #New path for cleaned csv file
                       original_col=NULL, #6-14 words for columns to be used in a specified order
                       set_col = c("mouse_ID", "acronym", "name", "counts", "volume.mm3", "normalized.count.by.volume") #you can add ("group", "sex", "genotype", "age", "strain", "drug", "reporter", "cohort") between mouse_ID and acronym
){
  # check if path and column names are defined correctly
  if (is.null(original_path) || is.null(set_path) || is.null(original_col)) {
    stop("original_path, set_path, and original_col must all be specified.")
  }

  if (length(original_col) != length(set_col)) {
    stop("original_col and set_col must be the same length.")
  }
  
  data <- read.csv(original_path, stringsAsFactors = FALSE)
  
  # check if factors of original_col exist in data
  if (!all(original_col %in% names(data))) {
    missing_cols <- original_col[!original_col %in% names(data)]
    stop(paste("These original_col names were not found in the CSV:", paste(missing_cols, collapse = ", ")))
  }
  
  # extract necessary columns
  subset_data <- data[original_col]
  
  # change column names
  colnames(subset_data) <- set_col
  
  # save as new CSV
  write.csv(subset_data, set_path, row.names = FALSE)
  
  message("New CSV saved to: ", set_path)
  return(invisible(subset_data))
}

#excute setcolumns
formatCsv(
  original_path = "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/data/Ex_639_Ch2_stitched_long_merge_Annotated_counts_with_leaf_with_density.csv",
  set_path = "/Users/saccyann/Documents/cleaned_data.csv",　
  original_col = c("ID", "Condition", "Sex", "acronym", "name", "newcounts", "newsize", "density; normalized"),
  set_col = c("mouse_ID", "group", "sex", "acronym", "name", "counts", "volume.mm3", "normalized.count.by.volume")
)
  
#create experiment object as opioid
opioid <- experiment(experiment_name = "opioid",
                         channels = "cfos",       # If you have more than one channel to import, set this to a character vector, e.g. c("cfos", "eyfp")
                         output_path = "/Users/saccyann/Documents/Sakura_networkanalysis") #Set this to a path location where you want your figures/analysis output to save, e.g. "P:\\DENNYLABV\\Michelle_Jin\\experiment\\folder"
print(opioid)

#import externally mapped datasets
opioid <- import_mapped_datasets(opioid, 
                                     normalized_count_paths = "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/data/cleaned_data.csv") #Edit this path",
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

#save
save(opioid, file = "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/opioid_labdata_check.RData")
save_experiment(opioid, timestamp = TRUE)
print("success")
