library(SMARTTR)
library(ggplot2)
library(ggraph)

# Load example dataset
load("/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR/anesthesia_labdata.RData") # Edit this path
# Replace the save path with a new save path on your computer
attr(anesthesia, "info")$output_path <- "/Users/saccyann/Documents/Sakura_networkanalysis/SMARTTR" #Edit this path

# Create correlation matrices
ontology <- "unified"
anesthesia <- get_correlations(anesthesia,
                       by = c("group"),
                       values = "Saline",
                       channels = "cfos",
                       p_adjust_method = "fdr",
                       ontology = ontology,
                       alpha =  0.05)

anesthesia <- get_correlations(anesthesia,
                               by = c("group"),
                               values = "Acute_Morphine",
                               channels = "cfos",
                               p_adjust_method = "fdr",
                               ontology = ontology,
                               alpha =  0.05)

anesthesia <- get_correlations(anesthesia,
                               by = c("group"),
                               values = "Chronic_Morphine",
                               channels = "cfos",
                               p_adjust_method = "fdr",
                               ontology = ontology,
                               alpha =  0.05)

anesthesia <- get_correlations(anesthesia,
                               by = c("group"),
                               values = "Chronic_Morphine_21",
                               channels = "cfos",
                               p_adjust_method = "fdr",
                               ontology = ontology,
                               alpha =  0.05)

anesthesia <- get_correlations(anesthesia,
                               by = c("group"),
                               values = "Withdrawal_Morphine",
                               channels = "cfos",
                               p_adjust_method = "fdr",
                               ontology = ontology,
                               alpha =  0.05)

anesthesia <- get_correlations(anesthesia,
                               by = c("group"),
                               values = "Withdrawal_Morphine_21",
                               channels = "cfos",
                               p_adjust_method = "fdr",
                               ontology = ontology,
                               alpha =  0.05)

print(anesthesia)
summary(anesthesia)

# Plot correlation matrices
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
#plot heatmaps
p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Saline",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5, # helps center and shift the significance asterisks
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Acute_Morphine",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5,
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Chronic_Morphine",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5, # helps center and shift the significance asterisks
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Chronic_Morphine_21",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5,
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = Withdrawal_Morphine",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5, # helps center and shift the significance asterisks
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Withdrawal_Morphine_21",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5,
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "Acute_Cocaine",
                                    correlation_list_name = Withdrawal_Morphine",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5, # helps center and shift the significance asterisks
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)

p_list <- plot_correlation_heatmaps(anesthesia,
                                    channels = "cfos",
                                    correlation_list_name = "Withdrawal_Cocaine",
                                    sig_color = "black",
                                    sig_nudge_y = -0.5,
                                    print_plot = FALSE,
                                    colors = c( "#ff2a04"),
                                    theme.hm = theme.hm,
                                    ontology = ontology,
                                    save_plot = TRUE)
# Create networks
anesthesia <- create_networks(anesthesia,
                             correlation_list_name = "Saline",
                             channels = c("cfos"),
                             alpha = 0.05,
                             ontology = ontology,
                             pearson_thresh = 0)

anesthesia <- create_networks(anesthesia,
                              correlation_list_name = "Acute_Morphine",
                              channels = "cfos",
                              alpha = 0.05,
                              ontology = ontology,
                              pearson_thresh = 0)

anesthesia <- create_networks(anesthesia,
                              correlation_list_name = "Saline",
                              channels = c("cfos"),
                              alpha = 0.05,
                              ontology = ontology,
                              pearson_thresh = 0)

anesthesia <- create_networks(anesthesia,
                              correlation_list_name = "Acute_Morphine",
                              channels = "cfos",
                              alpha = 0.05,
                              ontology = ontology,
                              pearson_thresh = 0)
# Plot networks
graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5, size = 50, face = "plain"),
                                             legend.text = element_text(size = 25),
                                             legend.title = element_text(size = 25))

anatomical.colors <- c(Isocortex = "#5571a9", OLF = "#64bdc4", HPF = "#d2875b", 
                        CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af", 
                        HY = "#8e7960",  MB = "#d796c8", HB = "#646464")

p_list <- plot_networks(anesthesia,
                        network_name = "Saline",
                        channels = "cfos",
                        edge_color = "#ff2a04",
                        label_size = 7,
                        label_offset = 0.13,
                        degree_scale_limit = c(1,20),
                        correlation_edge_width_limit = c(0.8, 1.0),
                        edge_thickness_range = c(1,4),
                        node_size_range = c(1, 6),
                        anatomical.colors = anatomical.colors,
                        graph_theme = graph_theme,
                        height = 20,
                        width = 20,
                        save_plot = TRUE)

p_list <- plot_networks(anesthesia,
                        network_name = "Acute_Morphine",
                        channels = "cfos",
                        edge_color = "#ff2a04",
                        label_size = 7,
                        label_offset = 0.13,
                        degree_scale_limit = c(1,20),
                        correlation_edge_width_limit = c(0.8, 1.0),
                        edge_thickness_range = c(1,4),
                        node_size_range = c(1, 6),
                        anatomical.colors = anatomical.colors,
                        graph_theme = graph_theme,
                        height = 20,
                        width = 20,
                        save_plot = TRUE)




