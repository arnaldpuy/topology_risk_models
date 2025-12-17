## ----setup, include=FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "pdf", cache = TRUE)


## ----warning=FALSE, message=FALSE, results = "hide"----------------------------------------------------

# PRELIMINARY FUNCTIONS #######################################################
################################################################################

sensobol::load_packages(c("data.table", "tidyverse", "openxlsx", "scales", 
                          "cowplot", "readxl", "ggrepel", "tidytext", "here", 
                          "tidygraph", "igraph", "foreach", "parallel", "ggraph", 
                          "tools", "purrr", "sensobol", "benchmarkme"))

# Create custom theme ----------------------------------------------------------

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA), 
          strip.background = element_rect(fill = "white"), 
          legend.text = element_text(size = 7.3), 
          axis.title = element_text(size = 10),
          legend.key.width = unit(0.4, "cm"), 
          legend.key.height = unit(0.4, "cm"), 
          legend.key.spacing.y = unit(0, "lines"),
          legend.box.spacing = unit(0, "pt"),
          legend.title = element_text(size = 7.3), 
          axis.text.x = element_text(size = 7), 
          axis.text.y = element_text(size = 7), 
          axis.title.x = element_text(size = 7.3), 
          axis.title.y = element_text(size = 7.3),
          plot.title = element_text(size = 8),
          strip.text.x = element_text(size = 7.4), 
          strip.text.y = element_text(size = 7.4)) 
}

# Select color palette ---------------------------------------------------------

color_languages <- c("fortran" = "steelblue", "python" = "lightgreen")

# Source all .R files in the "functions" folder --------------------------------

r_functions <- list.files(path = here("functions"), 
                          pattern = "\\.R$", full.names = TRUE)

lapply(r_functions, source)

# Set seed ---------------------------------------------------------------------

seed <- 123


## ----run_analysis, cache.lazy=FALSE--------------------------------------------------------------------

# CREATE DATASET ###############################################################

# Path to folder ---------------------------------------------------------------

path <- "./datasets/call_metrics"

# List CSV files ---------------------------------------------------------------

files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# Split by language ------------------------------------------------------------

python_files  <- grep("python",  files, value = TRUE, ignore.case = TRUE)
fortran_files <- grep("fortran", files, value = TRUE, ignore.case = TRUE)

base_fortran <- file_path_sans_ext(basename(fortran_files))
base_python <- file_path_sans_ext(basename(python_files))

model_names_fortran <- models <- sub(".*_", "", base_fortran)
model_names_python <- models <- sub(".*_", "", base_python)

# Load and name files ----------------------------------------------------------

python_list  <- lapply(python_files, fread)
fortran_list <- lapply(fortran_files, fread)

names(python_list) <- model_names_python
names(fortran_list) <- model_names_fortran

# RBIND ------------------------------------------------------------------------

make_callgraph <- function(lst, lang) {
  rbindlist(lst, idcol = "model") %>%
    .[, language := lang] %>%
    .[, .(model, language, `function`, call)] %>%
    setnames(., c("function", "call"), c("from", "to"))
}

python_callgraphs  <- make_callgraph(python_list,  "python") 
fortran_callgraphs <- make_callgraph(fortran_list, "fortran")

all_callgraphs <- rbind(python_callgraphs, fortran_callgraphs)

# Remove main and module calls -------------------------------------------------

all_callgraphs <- all_callgraphs[!(from %in% c("<module>", "main"))]

# LOAD CYCLOMATIC COMPLEXITY VALUES FOR FUNCTIONS AND SUBROUTINES ##############

cc_unique <- fread("./datasets/cyclomatic_complexity_functions.csv")

# CREATE NETWORK FROM CALL GRAPHS ##############################################

all_graphs <- all_callgraphs[, .(graph = list(as_tbl_graph(.SD, directed = TRUE))), 
                             .(model, language)]

# ADD NODE METRICS #############################################################

# Define the weights to characterize risky nodes -------------------------------

alpha <- 0.6  # Weight to cyclomatic complexity
beta  <- 0.3  # Weight to in-degree (impact of bug upstream)
gamma <- 0.1  # Weight to betweenness (critical bridge)

# Add node metrics -------------------------------------------------------------

all_graphs[, graph:= Map(function(g, m, lang) {
  
  comp_sub <- cc_unique[model == m & language == lang]
  
  # mean cyclomatic complexity for this model & language -----------------------
  
  mean_cyclo <- mean(comp_sub$cyclomatic_complexity, na.rm = TRUE)
  
  g %>%
    activate(nodes) %>%
    
    # Left join with dataset with cyclomatic complexity values -----------------
  
  left_join(comp_sub, by = "name") %>%
    
    # replace NA cyclomatic_complexity with model–language mean ----------------
  
    mutate(cyclomatic_complexity = if (!is.na(mean_cyclo)) {
      
        ifelse(is.na(cyclomatic_complexity), mean_cyclo, cyclomatic_complexity)
      
      } else {
        
        # if even the mean is NA (all NA in comp_sub), leave as-is
        
        cyclomatic_complexity
      }
    ) %>%
    
    # Remove Python MODULE_AGG / CLASS_AGG nodes from this graph 
    # because they are not callable --------------------------------------------
  
  filter(!(language == "python" & type %in% c("MODULE_AGG", "CLASS_AGG"))) %>%
    
    # Calculation of key network metrics ---------------------------------------
  
  mutate(indeg = centrality_degree(mode = "in"),
         outdeg = centrality_degree(mode = "out"),
         btw = centrality_betweenness(directed = TRUE, weights = NULL),
         cyclo_sc = rescale(cyclomatic_complexity),
         indeg_sc = rescale(indeg),
         btw_sc = rescale(btw), 
         risk_score = alpha * cyclo_sc + beta * indeg_sc + gamma * btw_sc)
         }, 
  graph, model, language)]

# EXTRACT NODE DF ##############################################################

all_graphs[, node_df:= lapply(graph, as_tibble, what = "nodes")]

# Export full node df ----------------------------------------------------------

full_node_df <- all_graphs %>%
  mutate(node_df = purrr::map(node_df, ~ select(.x, -model, -language))) %>%
  unnest(node_df) %>%
  select(-graph) %>%
  data.table()

write.xlsx(full_node_df, "full_node_df.xlsx")

# COMPUTE ALL PATHS AND THEIR RISK SCORES ######################################

all_graphs[, paths_tbl:= Map(all_paths_fun, node_df, graph)]

# Export full paths df ---------------------------------------------------------

full_paths_df <- all_graphs %>%
  unnest(paths_tbl) %>%
  select(-c(graph, node_df))

write.xlsx(full_paths_df, "full_paths_df.xlsx")

# CONDUCT UNCERTAINTY AND SENSITIVITY ANALYSIS #################################

# Define sample size and order of effects --------------------------------------

N <- 2^11
order <- "first"

# Run the function -------------------------------------------------------------

all_graphs[, uncertainty_sensitivity:= Map(full_ua_sa_risk_fun, node_df, paths_tbl, N, order)]


# UNNEST APPROPRIATELY #########################################################

unnested_df <- all_graphs %>%
  mutate(us_nodes = map(uncertainty_sensitivity, "nodes"),
         us_paths = map(uncertainty_sensitivity, "paths"))

# Create SA data frame ---------------------------------------------------------

full_sa_df <- unnested_df %>%
  select(us_nodes) %>% 
  unnest(cols = c(us_nodes)) %>%
  select(name, model, language, sensitivity_indices) %>%
  unnest(cols = c(sensitivity_indices))

# Export
fwrite(full_sa_df, "full_sa_df.csv")

# Create UA data frame ---------------------------------------------------------

full_ua_df <- unnested_df %>%
  select(model, language, us_paths) %>%
  unnest(cols = c(us_paths)) %>%
  data.table()

# Export
fwrite(full_ua_df, "full_ua_df.csv")


## ----some_stats, dependson="run_analysis"--------------------------------------------------------------

# CALCULATE SOME DESCRIPTIVE METRICS ###########################################

tmp <- data.table(full_paths_df)[, .(n_paths = .N), .(model, language)] %>%
  .[order(-n_paths)]

tmp2 <- data.table(full_node_df)[, .(n_nodes = .N), .(model, language)] %>%
  .[order(-n_nodes)]

# Path to node ratio: how interconnected the model is.
# Model_cc: Proxy for algorithmic complexity of model.
# Avg_path_length: Proxy for depth of dependency chains (risk-highway potential)
# Model fragility: more (error) propagation routes.
models_metrics <- merge(tmp, tmp2) %>%
  .[, `:=`(path_to_node_ratio = n_paths / n_nodes,
           model_cc = n_paths / log(n_nodes),
           avg_path_length = n_nodes / log(n_paths + 1),
           model_fragility_index = n_paths / (n_nodes * (n_nodes - 1)))] 

models_metrics


# Read descriptive_stats_file --------------------------------------------------

descriptive_stats <- data.table(read_xlsx("./datasets/descriptive_statistics/descriptive_statistics.xlsx"))
all_descriptive_df <- merge(models_metrics, descriptive_stats)

# Sort by model ----------------------------------------------------------------
  model_ordered <- all_descriptive_df[, sum(lines), model] %>%
  .[order(V1)]

# Plot descriptive measures per model ------------------------------------------

plot_descriptive <- melt(all_descriptive_df, measure.vars = c("lines_code", "n_nodes", "n_paths", 
                                          "path_to_node_ratio", "model_cc")) %>% 
  .[, model:= factor(model, levels = model_ordered[, model])] %>%
  ggplot(., aes(model, value, fill = language)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(breaks = breaks_pretty(n = 2)) +
  scale_fill_manual(values = color_languages, name = "") +
  facet_wrap(~ variable, ncol = 7, scales = "free_x") +
  labs(x = "", y = "N") +
  theme_AP() +
  theme(legend.position = c(0.1, 0.3))

plot_descriptive

# METRICS AT THE FILE AND FUNCTION LEVEL #######################################

folder <- "./datasets/results_per_function"

# Get names of files -----------------------------------------------------------

csv_files <- list.files(path = folder, pattern = "\\.csv$", full.names = TRUE)

# Split into file_metrics and func_metrics -------------------------------------

file_metric_files <- grep("file_metrics", csv_files, value = TRUE)
func_metric_files <- grep("func_metrics", csv_files, value = TRUE)

# Build one named list ---------------------------------------------------------

list_metrics <- list(file_metrics = setNames(lapply(file_metric_files, fread), 
                                             basename(file_metric_files)),
                     func_metrics = setNames(lapply(func_metric_files, fread), 
                                             basename(func_metric_files)))

# Create function to combine files ---------------------------------------------

make_combined <- function(subset_list, pattern) {
  rbindlist(subset_list[grep(pattern, names(subset_list))], idcol = "source_file")
}

# Combine files ----------------------------------------------------------------

metrics_combined <- list(file_fortran = make_combined(list_metrics$file_metrics, "fortran"),
                         file_python = make_combined(list_metrics$file_metrics, "python"),
                         func_fortran = make_combined(list_metrics$func_metrics, "fortran"),
                         func_python = make_combined(list_metrics$func_metrics, "python"))


# Functions to extract name of model and language from file --------------------

extract_model <- function(x) 
  sub("^(file|func)_metrics_\\d+_([A-Za-z0-9-]+)_(fortran|python).*", "\\2", x)

extract_lang  <- function(x) 
  sub("^(file|func)_metrics_\\d+_([A-Za-z0-9-]+)_(fortran|python).*", "\\3", x)

# Extract name of model and language -------------------------------------------

metrics_combined <- lapply(metrics_combined, function(dt) {
  dt[, source_file:= sub("\\.csv$", "", basename(source_file))]
  dt[, model:= extract_model(source_file)]
  dt[, language:= extract_lang(source_file)]
  dt
})

# Add column of complexity category --------------------------------------------

metrics_combined <- lapply(names(metrics_combined), function(nm) {
  dt <- as.data.table(metrics_combined[[nm]])
  if (grepl("^func_", nm) && "cyclomatic_complexity" %in% names(dt)) {
    dt[, complexity_category := cut(
      cyclomatic_complexity,
      breaks = c(-Inf, 10, 20, 50, Inf),
      labels = c("b1","b2","b3","b4")
    )]
  }
  dt
}) |> setNames(names(metrics_combined))

# Define labels ----------------------------------------------------------------

lab_expr <- c(b1 = expression(C %in% "(" * 0 * ", 10" * "]"),
              b2 = expression(C %in% "(" * 10 * ", 20" * "]"),
              b3 = expression(C %in% "(" * 20 * ", 50" * "]"),
              b4 = expression(C %in% "(" * 50 * ", " * infinity * ")"))

# Define vector to exclude classes that are not functions ----------------------

excluded_classes_vec <- c("MODULE_AGG", "CLASS_AGG")

# PLOT #########################################################################

## ----plot_c_model, dependson="read_metrics_function_data", fig.height=2.2, fig.width=3.1----

plot_c_model <- metrics_combined[grep("^func_", names(metrics_combined))] %>%
  lapply(., function(x) 
    x[, .(model, language, `function`, cyclomatic_complexity, loc, bugs, type)]) %>%
  rbindlist() %>%
  .[!type %in% excluded_classes_vec] %>%
  .[, model:= factor(model, levels = model_ordered[, model])] %>%
  ggplot(., aes(model, cyclomatic_complexity, fill = language, color = language)) +
  geom_boxplot(outlier.size = 0.7) +
  coord_flip() +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 2)) +
  facet_wrap(~language, scales = "free_x") +
  labs(x = "", y = "C") +
  theme_AP() +
  scale_color_manual(values = color_languages) +
  theme(legend.position = "none", 
        plot.margin = margin(0, 2, 0, 0))

plot_c_model


## ----plot_scatter_and_bar, dependson="read_metrics_function_data", fig.height=2.5, fig.width=3----

# Scatterplot cyclomatic vs lines of code --------------------------------------

plot_c_vs_loc <- metrics_combined[grep("^func_", names(metrics_combined))] %>%
  lapply(., function(x) x[, .(loc, cyclomatic_complexity, language)]) %>%
  rbindlist() %>%
  ggplot(., aes(loc, cyclomatic_complexity, color = language)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_x_continuous(breaks = breaks_pretty(n = 3)) +
  labs(x = "Lines code", y = "C") +
  scale_color_manual(values = color_languages) +
  theme_AP() + 
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  theme(legend.position = "none")

plot_c_vs_loc

# Count & proportion -----------------------------------------------------------

plot_bar_cyclomatic <- metrics_combined[grep("^func_", names(metrics_combined))] %>%
  lapply(., function(x) x[, .(complexity_category, language, type)]) %>%
  rbindlist() %>%
  .[!type %in% excluded_classes_vec] %>%
  .[, .N, .(complexity_category, language)] %>%
  .[, proportion := N / sum(N), language] %>%
  ggplot(., aes(complexity_category, proportion, fill = language)) +
  geom_bar(stat = "identity", position = position_dodge(0.6)) +
  scale_fill_manual(values = color_languages) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  scale_x_discrete(labels = lab_expr) +
  labs(x = "", y = "Proportion") +
  coord_flip() +
  theme_AP() +
  theme(legend.position = "none")

plot_bar_cyclomatic

plot_bar_category <- metrics_combined[grep("^func_", names(metrics_combined))] %>%
  lapply(., function(x) 
    x[, .(model, language, complexity_category, type)]) %>%
  rbindlist() %>%
  .[!type %in% excluded_classes_vec] %>%
  .[, model:= factor(model, levels = model_ordered[, model])] %>%
  .[, .N, .(model, language, complexity_category)] %>%
  .[, proportion := N / sum(N), .(language, model)] %>%
  ggplot(., aes(model, proportion, fill = complexity_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("yellowgreen", "orange", "red", "purple"), 
                    labels = lab_expr, 
                    name = "") +
  facet_wrap(~language) + 
  labs(x = "", y = "Proportion") +
  coord_flip() +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  theme_AP() + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank(), 
        legend.text = element_text(size = 7), 
        plot.margin = margin(0, 0, 0, 2))

plot_bar_category


## ----merge_descriptive_flow, dependson="some_stats", fig.height=3.9, fig.width=6-----------------------

# MERGE FIGURES ################################################################

legend2 <- get_legend_fun(plot_bar_category + theme(legend.position = "top"))
top_plot <- plot_grid(legend2, plot_descriptive, rel_heights = c(0.1, 0.9), ncol = 1, 
                      labels = "a")
bottom <- plot_grid(plot_c_vs_loc, plot_bar_cyclomatic, plot_c_model, 
                    plot_bar_category, ncol = 4, rel_widths = c(0.2, 0.24, 0.34, 0.22), 
                    labels = c("b", "c", "d"))
plot_grid(top_plot, bottom, ncol = 1, rel_heights = c(0.52, 0.48), align = "h",
  axis = "tb")


## ----plot_all_callgraphs, dependson="run_analysis", fig.height=2.5, fig.width=3------------------------

# PLOT FIGURES #################################################################

# Plot graphs ------------------------------------------------------------------

set.seed(seed)

# Thickness of edge:frequency across top-10 riski paths
# Color of edge: mean risk of paths using that edge
all_graphs <- all_graphs[, plot_obj:= mapply(plot_top_paths_fun, call_g = graph, 
                                             paths_tbl = paths_tbl, model.name = model,
                                             language = language, SIMPLIFY = FALSE)]


## ----plot, dependson="plot_all_callgraphs", fig.height=6, fig.width=6----------------------------------

# PLOT OTHER FIGURES ###########################################################

selected_models <- data.table(model = c("CTSM", "PCR-GLOBWB", "DBH", "HYPE", 
                                        "ORCHIDEE", "SWAT", "CWatM", "MHM"),
                              language = c("fortran", "python", "fortran", 
                                           "fortran", "fortran", "fortran", 
                                           "python", "python"))

# Plot call graphs -------------------------------------------------------------

tmp <- all_graphs[selected_models, on = .(model, language)] 
plot_all_risky_paths <- plot_grid(plotlist = tmp$plot_obj, ncol = 2, align = "hv")

# Plot risk_slope --------------------------------------------------------------

a <- full_paths_df %>%
  data.table() %>%
  .[selected_models, on = .(model, language)] %>%
  .[order(-p_path_fail), .SD[1:10], model] %>%
  ggplot(., aes(reorder(model, risk_slope), risk_slope)) +
  geom_boxplot() +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  geom_hline(yintercept = 0, lty = 2, color = "red") +
  coord_flip() +
  labs(x = "", y = expression(theta[1*k])) +
  theme_AP()

# Plot Gini metric -------------------------------------------------------------

b <- full_paths_df %>%
  data.table() %>%
  .[selected_models, on = .(model, language)] %>%
  .[order(-p_path_fail), .SD[1:10], model] %>%
  ggplot(., aes(reorder(model, gini_node_risk), gini_node_risk)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "", y = expression(G[k])) +
  theme_AP()

# Plot Si values ---------------------------------------------------------------

c <- full_sa_df %>%
  data.table() %>%
  .[selected_models, on = .(model, language)] %>%
  .[sensitivity == "Si", .(median = median(original, na.rm = TRUE)), .(model, parameters)] %>%
  ggplot(., aes(x = parameters, y = model, fill = median)) +
  geom_tile() +
  scale_fill_viridis_c(name = expression("Med(" * S[i] * ")"), 
                       limits = c(0, 1), 
                       breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = c(a_raw = expression(alpha),
                              b_raw = expression(beta),
                              c_raw = expression(gamma))) +
  labs(x = NULL, y = NULL) +
  theme_AP() +
  theme(legend.position = "none")

# Plot interaction strength ----------------------------------------------------


tmp <- split(full_sa_df, full_sa_df$language)
tmp2 <- lapply(tmp, data.table) %>%
  lapply(., function(x) {
    dcast(x, name + model + language + parameters ~ sensitivity,
          value.var = "original") %>%
      .[, interaction:= Ti - Si]})

d <- do.call(rbind, tmp2) %>%
  .[selected_models, on = .(model, language)] %>%
  .[, .(median = median(interaction, na.rm = TRUE)), .(parameters, model)] %>%
  ggplot(., aes(x = parameters, y = model, fill = median)) +
  geom_tile() +
  scale_x_discrete(labels = c(a_raw = expression(alpha),
                              b_raw = expression(beta),
                              c_raw = expression(gamma))) +
  scale_fill_viridis_c(name   = expression("Med(" * T[i] - S[i] * ")"), 
                       limits = c(0, 0.06), 
                       breaks = c(0, 0.03, 0.06), 
                       option = "magma") +
  labs(x = NULL, y = NULL) +
  theme_AP() +
  theme(legend.position = "none")
       

# MERGE #######################################################################

p_for_fill_legend <- all_graphs$plot_obj[[6]] +
  guides(size = "none", fill = guide_legend(title = ""))
fill_legend <- get_legend_fun(p_for_fill_legend + theme(legend.position = "top"))
plot_top_paths <- plot_grid(fill_legend, plot_all_risky_paths, ncol = 1, rel_heights = c(0.05 ,0.9), 
                            labels = "a")
heatmap_legend <- get_legend(c + theme(legend.position = "top"))
ti_legend <- get_legend(d + theme(legend.position = "top"))
dada <- plot_grid(a, b, c, d, ncol = 1, labels = c("b", "c", "d", "e"))
all_legends <- plot_grid(heatmap_legend, ti_legend, ncol = 1)
right_plot <- plot_grid(all_legends, dada, ncol = 1, rel_heights = c(0.1, 0.9))
plot_grid(plot_top_paths, right_plot, ncol = 2, rel_widths = c(0.7, 0.3))



## ----nodes_proportion, dependson="run_analysis", fig.height=3, fig.width=3.5---------------------------

# PATH-LEVEL RISK ACCOUNTED FOR THE TOP 5% NODES ###############################

setDT(full_paths_df)

# To long format ---------------------------------------------------------------

paths_long <- full_paths_df[, .(node = unlist(path_nodes),
                                p_path_fail = p_path_fail,
                                gini_node_risk = gini_node_risk,
                                risk_slope = risk_slope,
                                risk_mean = risk_mean,
                                risk_sum  = risk_sum), 
                            .(model, language, path_id)]


# Aggregate at function level --------------------------------------------------

node_from_paths <- paths_long[, .(n_paths = .N,
                                  mean_p_path = mean(p_path_fail, na.rm = TRUE),
                                  max_p_path = max(p_path_fail,  na.rm = TRUE),
                                  sum_p_path = sum(p_path_fail,  na.rm = TRUE),
                                  mean_gini = mean(gini_node_risk, na.rm = TRUE),
                                  mean_slope = mean(risk_slope, na.rm = TRUE),
                                  mean_risksum = mean(risk_sum, na.rm = TRUE)),
                              .(model, language, node)]

# Join with nodes --------------------------------------------------------------

node_summary <- merge(node_from_paths, full_node_df, by.x = c("model", "language", "node"),
                      by.y = c("model", "language", "name"), all.x = TRUE)

# Calculate risk mass ----------------------------------------------------------

node_summary[, risk_mass:= mean_p_path * n_paths]

# share of risk mass in top X% nodes, per model .-------------------------------


top_share <- function(X = 0.05) {
  
  node_summary[!is.na(risk_mass) & risk_mass >= 0, {
    dt <- .SD[order(-risk_mass)]
    n_top <- max(1L, ceiling(.N * X))
    .(X = X, n_nodes = .N, n_top = n_top, 
      share_risk_mass_topX = sum(dt$risk_mass[1:n_top]) / sum(dt$risk_mass))
  },
  .(model, language)
  ]
}

# Run function -----------------------------------------------------------------

tmp <- top_share(0.05) %>%
  .[order(-share_risk_mass_topX)]

tmp

# Plot--------------------------------------------------------------------------

merge(all_descriptive_df, tmp, by = c("model", "language")) %>%
  ggplot(., aes(share_risk_mass_topX, path_to_node_ratio, color = language)) +
  geom_point() +
  scale_color_manual(values = color_languages, name = "") +
  geom_text_repel(aes(label = model), size = 2, max.overlaps = Inf, show.legend = FALSE) +
  scale_y_log10() +
  labs(x = "Proportion path-level risk top 5% nodes", y = "path_to_node_ratio") +
  theme_AP() +
  theme(legend.position = c(0.2, 0.8))


## ----plot_paths, dependson="run_analysis"--------------------------------------------------------------

# PLOT THE TOP 50 PATHS PER MODEL ##############################################

tmp <- full_ua_df %>%
  .[order(-P_k_mean), .SD[1:50], .(model, language)] %>%
  split(., list(.$model, .$language)) %>%
  lapply(., na.omit)

# Remove empty slots ----------------------------------------------------------

tmp2 <- tmp[sapply(tmp, function(x) nrow(x) > 0)]
  
# Plot in a for loop ----------------------------------------------------------

out <- list()

for ( i in 1:length(tmp2)) {
  
  out[[i]] <- ggplot(tmp2[[i]], aes(P_k_mean, reorder(path_str, P_k_mean), color = risk_slope))  +
    geom_point(size = 1) +
    geom_errorbar(aes(xmin = P_k_q025, xmax = P_k_q975), height = 0.2) +
    scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0,             
                          name = expression(beta[k])) +
    labs(y = "Path ID", x = expression(P[k])) +
    theme_AP() +
    scale_x_continuous(breaks = breaks_pretty(n = 3), 
                       limits = c(0, 1)) +
    theme(axis.text.y = element_text(size = 4), 
          legend.position = "top") +
    ggtitle(names(tmp[i]))
  
}

out


## ----session_information-------------------------------------------------------------------------------

# SESSION INFORMATION ##########################################################

sessionInfo()

## Return the machine CPU ------------------------------------------------------

cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores -------------------------------------------------

cat("Num cores:   "); print(detectCores(logical = FALSE))

## Return number of threads ---------------------------------------------------

cat("Num threads: "); print(detectCores(logical = FALSE))






full_paths_df[order(-p_path_fail), .SD[1:10], .(model, language)] %>%
  .[, .(model, language, path_str)]

