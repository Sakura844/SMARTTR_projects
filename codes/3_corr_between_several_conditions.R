library(SMARTTR)
library(ggplot2)
library(ggraph)

# Replace the save path with a new save path on your computer
load("/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/opioid_labdata_new.RData") # Edit this path
attr(opioid, "info")$output_path <- "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR" #Edit this path
output_dir <- "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/figures/Parallel/testpdf" #Edit this path

# Create correlation matrices
ontology <- "unified"
groups <- c("Saline", "Acute_Morphine", "Chronic_Morphine", 
            "Chronic_Morphine_21", "Withdrawal_Morphine", "Withdrawal_Morphine_21")
pairs <- combn(groups, 2)

# Create the set of names "pair1_vs_pair2"
results_name <- apply(pairs, 2, function(x) paste(x[1], "vs", x[2], sep = "_"))

# Get regional cross correlations and their p-values in a correlation list object
print(results_name)
for (g in groups) {
  opioid <- get_correlations(opioid,
                                 by = c("group"),
                                 values = g,
                                 channels = "cfos",
                                 p_adjust_method = "fdr",
                                 ontology = ontology,
                                 alpha =  0.05)
}

# Generate all pairs
group_pairs <- combn(groups, 2, simplify = FALSE)
print(group_pairs)

#Permutation analysis to compare the region pairwise correlation coefficients between two different analysis groups.
results <- list()

for (pair in group_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  result_name <- paste0(group1, "_vs_", group2)
  opioid <- correlation_diff_permutation(
    opioid,
    correlation_list_name_1 = group1,
    correlation_list_name_2 = group2,
    channels = c("cfos"),
    seed = 5,
    p_adjust_method = "none",
    alpha = 0.05
  )
  #Export the permutation results as a csv file: not necessary
  opioid <- export_permutation_results(opioid, permutation_groups = "all", filter_significant = TRUE)
}
print(results_name)

# Fine tuning the graph aesthetics
theme.hm <- ggplot2::theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 8),
                           axis.text.y = element_text(vjust = 0.5, size = 8),
                           plot.title = element_text(hjust = 0.5, size = 36),
                           axis.title = element_text(size = 18),
                           legend.text = element_text(size = 22),
                           legend.key.height = unit(100, "points"),
                           legend.title = element_text(size = 22),
                           panel.spacing = unit(0.2, "lines"),
                           strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
                           strip.text.y = element_text(angle = 270,  hjust = 0.5, vjust = 0.5, size = 10),
                           strip.placement = "outside", 
                           strip.background = element_rect(color = "black", fill = "grey"))

# volcano_plot
for (pair in results_name) {
    volcano_plot(
      opioid, 
      permutation_comparison = pair,
      channels = c("cfos"), 
      colors = c("#be0000"), 
      save_plot = TRUE, 
      title = NULL, 
      ylim = c(0, 3), 
      height = 8, 
      width = 10, 
      print_plot = FALSE, 
      image_ext = ".pdf"
    )
}
                      
#If you need pdf of graphs, you can run the following codes.
#extensions <- c(".png", ".pdf")
# for (pair in results_name) {
#   for (ext in extensions){
#   volcano_plot(
#      opioid,
#       permutation_comparison = pair,
#  channels = c("cfos"),
#   colors = c("#be0000"),
#    save_plot = TRUE,
#      title = NULL,
#     ylim = c(0, 3),
#       height = 8,
#  width = 10,
#   print_plot = FALSE,
#    image_ext = ext
#   )
#  }
# }

# parallel coordinate plot
for (pair in results_name) {
  tryCatch({
    file_name <- paste0(output_dir, "/", pair, ".pdf")
    
    # Excute
    p_list <- parallel_coordinate_plot_small(
      opioid, 
      permutation_comparison = pair,
      channels = "cfos", 
      colors = "#be0000",
      save_plot = FALSE #don't save within the function
    )
    
    # Save as pdf
    ggsave(file_name, plot = p_list[[1]], width = 10, height = 8)
    
    # Check
    message(paste("Successfully saved:", file_name))
    
  }, error = function(e) {
    # generate err message and go to the next loop
    message(paste("Error with", pair, ":", e$message))
  })
}
#Save as PNG
for (pair in results_name) {
  tryCatch({
    file_name <- paste0(output_dir, "/", pair, ".png")
    
    # Excute 
    p_list <- parallel_coordinate_plot_small(opioid, 
                                             permutation_comparison = pair,
                                             channels = "cfos", 
                                             colors ="#be0000")
    
    # Save
    ggsave(file_name, plot = p_list[[1]], width = 10, height = 8)
    
    # Check
    message(paste("Successfully saved:", file_name))
    
  }, error = function(e) {
    # generate err message and go to the next loop
    message(paste("Error with", pair, ":", e$message))
  })
}

# Create networks
opioid <- create_joined_networks(opioid,
                                       correlation_list_name = c("Saline", "Acute_Morphine"),
                                       channels = c("cfos"),
                                       ontology = ontology,
                                       alpha = 0.05)

# Plot networks
graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5, size = 50, face = "plain"),
                                             legend.text = element_text(size = 25),
                                             legend.title = element_text(size = 25))

anatomical.colors <- c(Isocortex = "#5571a9", OLF = "#64bdc4", HPF = "#d2875b", 
                       CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af", 
                       HY = "#8e7960",  MB = "#d796c8", HB = "#646464")
p_list <- parallel_coordinate_plot(opioid,
                                   permutation_comparison = "Context_vs_Shock",
                                   channels = "cfos", colors ="#be0000", x_label_group1 = "Context", x_label_group_2 = "Shock")
p_list <- plot_joined_networks(opioid, 
                               correlation_list_names = c("Saline", "Acute_Morphine"),
                               channels = "cfos", 
                               edge_colors = c(male_agg =  "#06537f", female_non = "#C70039"),
                               edge_color_labels = c(male_agg = "Saline_1", female_non = "Acute_Morphine_2"),
                               degree_scale_limit = c(1,45), correlation_edge_width_limit = c(0.8, 1.0),
                               height = 30, width = 30, label_size = 13, label_offset = 0.08, image_ext = ".png")





#Redefined the function to make letters small
parallel_coordinate_plot_small <- function(e,
                                     permutation_comparison = "AD_vs_control",
                                     channels = c("cfos", "eyfp", "colabel"),
                                     colors =  c("#be0000", "#00782e", "#f09b08"),
                                     x_label_group_1 = NULL,
                                     x_label_group_2 = NULL,
                                     height = 10,
                                     width = 10,
                                     print_plot = FALSE,
                                     save_plot = TRUE,
                                     reverse_group_order= FALSE,
                                     force = 1,
                                     plt_theme = NULL,
                                     label_size = 5,
                                     image_ext = ".png",
                                     nudge_x = 2:5
){
  if(is.null(x_label_group_1)){
    group_1 <- stringr::str_split(permutation_comparison, "_vs_", simplify = TRUE)[,1] %>%
      stringr::str_split("_", simplify = TRUE) %>% paste(collapse = " ")
  } else {
    group_1 <- x_label_group_1
  }
  if(is.null(x_label_group_2)){
    group_2 <- stringr::str_split(permutation_comparison, "_vs_", simplify = TRUE)[,2] %>%
      stringr::str_split("_", simplify = TRUE) %>% paste(collapse = " ")
  } else {
    group_2 <- x_label_group_2
  }
  
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels
  
  for (k in 1:length(channels)){
    alpha <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$alpha
    p_vals <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$p_val %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    corr_diffs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$test_statistic %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    sigs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$sig %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    group_1_pearson <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$group_1_pearson %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    group_2_pearson <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$group_2_pearson %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    
    p_vals <- p_vals %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                             values_to = "p_val", names_to = "colreg")
    corr_diffs <- corr_diffs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                                     values_to = "corr_diff", names_to = "colreg")
    sigs <- sigs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                         values_to = "sig", names_to = "colreg")
    group_1_pearson <- group_1_pearson %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                                               values_to = group_1, names_to = "colreg")
    group_2_pearson <- group_2_pearson %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                                               values_to = group_2, names_to = "colreg")
    
    df <- p_vals %>% dplyr::inner_join(corr_diffs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(sigs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_1_pearson, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_2_pearson, by = c("rowreg", "colreg")) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(c(group_1, group_2)), names_to = "group", values_to = "corr")
    
    df <- df %>% dplyr::filter(.data$sig, abs(.data$corr_diff) >= 1) %>%
      dplyr::mutate(group = factor(.data$group, levels = c(group_1, group_2)),
                    nudge = ifelse(.data$group == group_1, -0.1, 0.1)) %>%
      dplyr::arrange(.data$group, .data$corr_diff) %>%
      mutate(text = paste(.data$rowreg, .data$colreg, sep = "."),
             group_plot = paste(.data$rowreg, .data$colreg, sep = "."))
    
    if (isTRUE(reverse_group_order)){
      df <- df %>% dplyr::mutate(group = forcats::fct_rev(.data$group),
                                 nudge = ifelse(.data$group == group_2, -0.1, 0.1))
    }
    
    tryCatch({
      df[seq(2, nrow(df)/2, by = 2),]$text <- ""
      df[seq(nrow(df)/2+1, nrow(df), by = 2),]$text <- ""
      
    }, error = function(e) {
      message("could not divide by 2. Skipping")
    })
    
    if (is.null(plt_theme)){
      plt_theme <- ggplot2::theme_classic() +
        theme(text = element_text(size = 22),
              line = element_line(linewidth = 1),
              plot.title = element_text(hjust = 0.5, size = 36),
              axis.ticks.length = unit(5.5, "points"),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black")
        )
    }
    
    p <- ggplot(df, aes(x = .data$group, y = .data$corr, group = .data$group_plot)) +
      ggplot2::geom_line(alpha = 0.5, color = colors[k], linewidth = 3) +
      ggplot2::geom_point(size = 4, alpha = 0.5, color = colors[k]) +
      ggrepel::geom_text_repel(aes(label = .data$text),
                               size = label_size,
                               color = colors[k], direction = "y",
                               force = force,
                               ylim = c(-1, 1),
                               segment.alpha = 0.3,
                               nudge_x = dplyr::pull(df, .data$nudge)*nudge_x, max.iter = 20000) +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=1.2) +
      xlab("Group") + ylab("Correlation") +
      ggplot2::expand_limits(y=c(-1,1)) + plt_theme
    
    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }
    
    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0(permutation_comparison, "_parallel_coordinate_plot_",
                                                 channels[k], "_", image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

