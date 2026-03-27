## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "pdf", cache = TRUE)


## ----warning=FALSE, message=FALSE, results = "hide"-------------------------------------------

# PRELIMINARY FUNCTIONS ########################################################
###############################################################################

sensobol::load_packages(c("data.table", "tidyverse", "openxlsx", "scales", 
                          "cowplot", "readxl", "ggrepel", "tidytext", "here", 
                          "tidygraph", "igraph", "foreach", "parallel", "ggraph", 
                          "tools", "purrr", "sensobol", "benchmarkme", "softwareRisk"))

# Select color palette ---------------------------------------------------------

color_languages <- c("fortran" = "steelblue", "python" = "lightgreen")

# Color to functional forms ----------------------------------------------------

color_functional_forms <- c("additive" = "#C65D09", "compensatory" = "#3B5B92")

# Set seed ---------------------------------------------------------------------

seed <- 123

# Source all .R files in the "functions" folder --------------------------------

r_functions <- list.files(path = here("functions"), pattern = "\\.R$", full.names = TRUE)
lapply(r_functions, source)


## ----run_analysis, cache.lazy=FALSE-----------------------------------------------------------

# CREATE DATASET ##############################################################

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
    .[, .(file, model, language, `function`, call)] %>%
    setnames(., c("function", "call"), c("from", "to"))
}

python_callgraphs  <- make_callgraph(python_list,  "python") 
fortran_callgraphs <- make_callgraph(fortran_list, "fortran")

all_callgraphs <- rbind(python_callgraphs, fortran_callgraphs)

# SOURCE CODE CLASSIFICATION BY FUNCTIONAL ROLE ################################

# Strip leading "./models/"-----------------------------------------------------
all_callgraphs[, file_clean:= sub("^\\./models/", "", file)]

# Provenance: file must be inside ./models to be eligible at all----------------
all_callgraphs[, in_model_tree:= !is.na(file) & nchar(file) > 0L & 
                 grepl("^\\./models/", file)]

# Rest = everything after first segment-----------------------------------------
all_callgraphs[, rest:= sub("^[^/]+/", "", file_clean)]

# If rest starts with model_id/ then drop it (this is the duplicate)------------
all_callgraphs[, rest:= fifelse(startsWith(rest, paste0(tolower(model), "/")),
                                sub("^[^/]+/", "", rest), rest)]

# In case of triple nesting like vic/vic/vic/...--------------------------------
all_callgraphs[, rest:= fifelse(startsWith(rest, paste0(tolower(model), "/")),
                                sub("^[^/]+/", "", rest), rest)]

# Top_level---------------------------------------------------------------------
all_callgraphs[, top_level:= tstrsplit(rest, "/", fixed = TRUE, keep = 1L)]

# Useful generic external patterns (works across models)------------------------
all_callgraphs[, is_generic_external:= in_model_tree & (
  grepl("(^|/)(extern|external|externals|third[_-]?party|vendor|vendored)(/|$)", 
        rest, ignore.case = TRUE) |
    grepl("(^|/)\\.lib(/|$)", rest, ignore.case = TRUE))]

# CTSM-specific external (framework + FoX parser under cdeps + CESM shared infra)
all_callgraphs[, is_ctsm_external:= (model == "CTSM") & in_model_tree & (
  grepl("^cime(/|$)", rest) |
    grepl("^cime_config(/|$)", rest) |
    grepl("^components/cdeps(/|$)", rest) |      # FoX XML/SAX parser
    grepl("(^|/)share_esmf(/|$)", rest) |        # shared ESMF infrastructure
    grepl("^src/unit_test_shr(/|$)", rest)       # unit tests 
)]

# Extra CTSM allowlist (tightens "model core" to land-model code)---------------
CTSM_STRICT <- TRUE

all_callgraphs[, is_ctsm_allowed:= TRUE]

if (CTSM_STRICT) {
  all_callgraphs[model == "CTSM" & in_model_tree, is_ctsm_allowed :=
                   grepl("^src(/|$)", rest) |
                   grepl("^lilac(/|$)", rest) |            
                   grepl("^components/cism(/|$)", rest)]
}

## Component classification (order matters)-------------------------------------
all_callgraphs[nchar(file) == 0L | is.na(file), component:= NA_character_]

all_callgraphs[!is.na(file) & nchar(file) > 0L, component:= fcase(
  
  # Anything not in ./models is external immediately
  !in_model_tree, "external_lib",
  
  # Generic external dirs (vendored/third_party/etc.)
  is_generic_external, "external_lib",
  
  # CI and vendored libraries 
  top_level == ".github", "ci_cd",
  top_level == ".lib" | grepl("/\\.lib/", rest), "vendored_lib",
  
  # CIME / CESM / CDEPS infrastructure (framework) 
  grepl("^cime/CIME/", rest), "framework",
  grepl("^cime_config/", rest), "framework",
  grepl("^cime/doc/", rest), "framework",
  grepl("^components/cdeps/cime_config/", rest), "framework",
  
  # CTSM-specific: treat cime + cdeps (incl FoX) as external/framework
  is_ctsm_external, "framework",
  
  # CTSM shared "shr_*" routines (framework by symbol)
  (model == "CTSM") & (grepl("^shr_", from) | grepl("^shr_", to)), "framework",
  
  # Tests (incl. SystemTests, case-insensitive) 
  grepl("SystemTests/", rest, ignore.case = TRUE), "tests",
  grepl("/tests?/|^tests?/", rest, ignore.case = TRUE), "tests",
  grepl("(^|/)unit_test", rest, ignore.case = TRUE), "tests",
  
  # Drivers: Fortran code clearly in drivers directories 
  grepl("(^|/)drivers(/|$)", rest) & language == "fortran", "driver",
  
  # Couplers: NUOPC / LILAC / CESM / cpl 
  grepl("cpl_|/cesm/|/cpl/|/cpl_", rest), "coupler",
  
  # CLI / tools: scripts, setup, CTSM python utilities 
  grepl("tools/|scripts/|cli\\.py$|setup\\.py$", rest), "cli_or_tool",
  grepl("^python/ctsm/", rest), "cli_or_tool",
  
  # Everything else: model core 
  default = "model_core"
)]

# include flag for risk calculations ----------------------------------------
# Key changes:
# - require in_model_tree
# - exclude external_lib/framework/tests/etc (already via component filter)
# - CTSM strict allowlist (optional)
all_callgraphs[, include_in_risk := (
  in_model_tree &
    language %chin% c("fortran", "python") &
    component %chin% c("model_core", "coupler") &
    (model != "CTSM" | !CTSM_STRICT | is_ctsm_allowed)
)]

# --- Sanity check: should now be FALSE for the SAX/FoX functions in CTSM
all_callgraphs[
  model == "CTSM" & include_in_risk &
    (grepl("sax|parseDTD|parsefile|parsestring|runParser", from, ignore.case = TRUE) |
       grepl("sax|parseDTD|parsefile|parsestring|runParser", to,   ignore.case = TRUE)),
  .(from, to, file, rest, component)
]

# SUmmarize --------------------------------------------------------------------

all_callgraphs[, .N, component]
all_callgraphs[, .N, include_in_risk]

# Remove module calls -------------------------------------------------

all_callgraphs <- all_callgraphs[!(from %in% "<module>")] %>%
  .[include_in_risk == TRUE]

# LOAD CYCLOMATIC COMPLEXITY VALUES FOR FUNCTIONS AND SUBROUTINES ##############

cc_unique <- fread("./datasets/cyclomatic_complexity_functions.csv")

# CREATE NETWORK FROM CALL GRAPHS ##############################################

all_graphs <- all_callgraphs[, .(graph = list(as_tbl_graph(.SD, directed = TRUE))), 
                             model]

# ADD NODE METRICS #############################################################

# Define the weights to characterize risky nodes -------------------------------

alpha <- 0.6  # Weight to cyclomatic complexity
beta  <- 0.3  # Weight to in-degree (impact of bug upstream)
gamma <- 0.1  # Weight to betweenness (critical bridge)

# Add node metrics -------------------------------------------------------------

all_graphs[, graph:= Map(function(g, m) {
  
  comp_sub <- cc_unique[model == m]
  
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
  
  mutate(type = type, 
         indeg = centrality_degree(mode = "in"),
         outdeg = centrality_degree(mode = "out"),
         btw = centrality_betweenness(directed = TRUE, weights = NULL),
         cyclo_sc = rescale(cyclomatic_complexity),
         indeg_sc = rescale(indeg),
         btw_sc = rescale(btw), 
         risk_score = alpha * cyclo_sc * beta * indeg_sc * gamma * btw_sc)
}, 
graph, model)]

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

# list of risk-form specifications ---------------------------------------------
# We run an additive and a power-mean version with p = -1

risk_specs <- list(additive = list(risk_form = "additive", p = NA),
                   power_mean = list(risk_form = "power_mean", p = -1))

# Run programatically ----------------------------------------------------------

for (nm in names(risk_specs)) {
  spec <- risk_specs[[nm]]
  
  all_graphs[, (paste0("paths_tbl_", nm)):= Map(all_paths_fun, graph = graph,
                                                risk_form = spec$risk_form,
                                                alpha = alpha,
                                                beta  = beta,
                                                gamma = gamma,
                                                p  = spec$p,                
                                                complexity_col = "cyclomatic_complexity")]}

# Unnest tables ----------------------------------------------------------------

# For additive form -----------------

nodes_additive <- unnest_paths_tbl_fun(all_graphs, "paths_tbl_additive", "nodes") %>%
  .[, risk_form := "additive"]
paths_additive <- unnest_paths_tbl_fun(all_graphs, "paths_tbl_additive", "paths") %>%
  .[, risk_form := "additive"]

# For power mean form ---------------

nodes_power_mean <- unnest_paths_tbl_fun(all_graphs, "paths_tbl_power_mean", "nodes") %>%
  .[, risk_form := "power_mean"]
paths_power_mean <- unnest_paths_tbl_fun(all_graphs, "paths_tbl_power_mean", "paths") %>%
  .[, risk_form := "power_mean"]

# Rbind ------------------------------------------------------------------------

nodes_all <- rbind(nodes_additive, nodes_power_mean, fill = TRUE)
paths_all <- rbind(paths_additive, paths_power_mean, fill = TRUE)[, path_str:= NULL]

# Export -----------------------------------------------------------------------

write.xlsx(nodes_all, "full_nodes_df.xlsx")
write.xlsx(paths_all, "full_paths_df.xlsx")


## ----ua_sa, dependson="run_analysis"----------------------------------------------------------

# CONDUCT UNCERTAINTY AND SENSITIVITY ANALYSIS #################################

# Define sample size and order of effects --------------------------------------

N <- 2^10
order <- "first"

# Run the function (we remove the vic and python model implementation because 
# there are not paths) ---------------------------------------------------------

idx <- all_graphs[model != "VIC", which = TRUE]

# Define variants (input column -> output column + extra args) -----------------

ua_specs <- list(additive = list(paths_col = "paths_tbl_additive",
                                 out_col = "uncertainty_sensitivity_additive",
                                 extra = list()),
                 power_mean = list(paths_col = "paths_tbl_power_mean",
                                   out_col = "uncertainty_sensitivity_power_mean",
                                   extra = list(risk_form = "power_mean")))

# Define parallel settings -----------------------------------------------------

ncores <- max(1L, parallel::detectCores() - 1L)

# Run UA/SA --------------------------------------------------------------------

for (nm in names(ua_specs)) {
  spec <- ua_specs[[nm]]
  
  all_graphs[idx, (spec$out_col):= parallel::mclapply(X = get(spec$paths_col),
                                            FUN = function(x) 
                                              do.call(uncertainty_fun,
                                                      c(list(all_paths_out = x, 
                                                             N = N, 
                                                             order = order), 
                                                        spec$extra)),
                                            mc.cores = ncores
  )]
}

# Extract data -----------------------------------------------------------------

# Additive -----------------------

nodes_us_additive <- unnest_us_fun(all_graphs, 
                                   us_col = "uncertainty_sensitivity_additive",
                                   slot   = "nodes") %>%
  .[, risk_form:= "additive"]

paths_us_additive <- unnest_us_fun(all_graphs, 
                                   us_col = "uncertainty_sensitivity_additive",
                                   slot = "paths") %>%
  .[, risk_form:= "additive"]

# Power-mean --------------------

nodes_us_power_mean <- unnest_us_fun(all_graphs,
                                     us_col = "uncertainty_sensitivity_power_mean",
                                     slot   = "nodes") %>%
  .[, risk_form:= "power_mean"]

paths_us_power_mean <- unnest_us_fun(all_graphs,
                                     us_col = "uncertainty_sensitivity_power_mean",
                                     slot   = "paths") %>%
  .[, risk_form:= "power_mean"]

# Rbind ------------------------------------------------------------------------

nodes_us_all <- rbind(nodes_us_additive, nodes_us_power_mean, fill = TRUE)
paths_us_all <- rbind(paths_us_additive, paths_us_power_mean, fill = TRUE)

# Calculate median, q5 and q95 -------------------------------------------------
    
cols <- c("uncertainty_analysis", "risk_trend", "gini_index")
    
for (nm in cols) {
  
  prefix <- switch(nm, uncertainty_analysis = "P", risk_trend = "risk",
                   gini_index = "gini")
  paths_us_all[, paste0(prefix, "_median"):= vapply(get(nm), median, numeric(1), na.rm = TRUE)]
  paths_us_all[, paste0(prefix, "_q5"):= vapply(get(nm), quantile, numeric(1), probs = 0.05, na.rm = TRUE)]
  paths_us_all[, paste0(prefix, "_q95"):= vapply(get(nm), quantile, numeric(1), probs = 0.95, na.rm = TRUE)]

}

# Extract sensitivity indices --------------------------------------------------

indices_all <- extract_sa_fun(nodes_us_all)


## ----some_stats, dependson="run_analysis"-----------------------------------------------------

# CALCULATE SOME DESCRIPTIVE METRICS ###########################################

tmp <- data.table(paths_all[risk_form == "additive"])[, .(n_paths = .N), model] %>%
  .[order(-n_paths)]

tmp2 <- data.table(nodes_all[risk_form == "additive"])[, .(n_nodes = .N), model] %>%
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

num_cols <- c("files","functions","modules","lines","lines_code","lines_comments")

descriptive_stats <- data.table(read_xlsx("./datasets/descriptive_statistics/descriptive_statistics.xlsx")) 

descriptive_stats <- dcast(melt(descriptive_stats, id.vars="model", measure.vars=num_cols),
                           model ~ variable, value.var ="value", fun.aggregate = function(z) sum(z, na.rm = TRUE)) %>%
  .[, lines_function := lines_code / functions]

all_descriptive_df <- merge(models_metrics, descriptive_stats)

# Sort by model ----------------------------------------------------------------
model_ordered <- all_descriptive_df[, sum(lines), model] %>%
  .[order(V1)]

# Plot descriptive measures per model ------------------------------------------

plot_descriptive <- melt(all_descriptive_df, measure.vars = c("lines_code", "n_nodes", "n_paths", 
                                                              "path_to_node_ratio", "model_cc")) %>% 
  .[, model:= factor(model, levels = model_ordered[, model])] %>%
  ggplot(., aes(model, value)) +
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

plot_c_model <- metrics_combined[grep("^func_", names(metrics_combined))] %>%
  lapply(., function(x) 
    x[, .(model, language, `function`, cyclomatic_complexity, loc, bugs, type)]) %>%
  rbindlist() %>%
  .[!type %in% excluded_classes_vec] %>%
  .[, model:= factor(model, levels = model_ordered[, model])] %>%
  na.omit() %>%
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


## ----merge_descriptive_flow, dependson="some_stats", fig.height=3.9, fig.width=6--------------

# MERGE FIGURES ################################################################

legend2 <- get_legend_fun(plot_bar_category + theme(legend.position = "top"))
top_plot <- plot_grid(legend2, plot_descriptive, rel_heights = c(0.1, 0.9), ncol = 1, 
                      labels = "a")
bottom <- plot_grid(plot_c_vs_loc, plot_bar_cyclomatic, plot_c_model, 
                    plot_bar_category, ncol = 4, rel_widths = c(0.2, 0.24, 0.34, 0.22), 
                    labels = c("b", "c", "d"))
plot_grid(top_plot, bottom, ncol = 1, rel_heights = c(0.52, 0.48), align = "h",
          axis = "tb")


## ----plot_all_callgraphs, dependson="run_analysis", fig.height=2.5, fig.width=3---------------

# PLOT FIGURES #################################################################

# Define language of models ----------------------------------------------------

python_models <- c("CWatM", "GR4J", "HBV", "PCR-GLOBWB")
fortran_models <- c("ORCHIDEE", "H08", "HYPE", "DBH", "SACRAMENTO")
python_and_fortran <- c("CTSM", "SWAT", "MHM", "HydroPy", "VIC")

all_graphs[, language := fcase(model %chin% python_models, "python",
                               model %chin% fortran_models, "fortran",
                               model %chin% python_and_fortran, "python+fortran",
                               default = NA_character_)]


# Define specifications --------------------------------------------------------

plot_specs <- list(additive = list(paths_col = "paths_tbl_additive", 
                                   out_col = "plot_obj_additive"),
                   power_mean = list(paths_col = "paths_tbl_power_mean",  
                                     out_col = "plot_obj_power"))

# Plot -------------------------------------------------------------------------

# Thickness of edge:frequency across top-10 risky paths
# Color of edge: mean risk of paths using that edge

# Plot graphs ------------------------------------------------------------------

set.seed(seed)

for (nm in names(plot_specs)) {
  spec <- plot_specs[[nm]]
  
  all_graphs[, (spec$out_col):= mapply(FUN = plot_top_paths_fun, graph = graph,
                                       all_paths_out= get(spec$paths_col),
                                       model.name = model, language = language,
                                       SIMPLIFY = FALSE)]
}

# Remove the fill legend -------------------------------------------------------

remove_fill_legend <- function(p) {
  p + guides(fill = "none")
}

plot_cols <- grep("^plot_obj_", names(all_graphs), value = TRUE)

for (col in plot_cols) {
  all_graphs[, (col) := lapply(get(col), remove_fill_legend)]}

# Define models to plot --------------------------------------------------------

selected_models <- data.table(model = c("CTSM", "PCR-GLOBWB", "DBH", "HYPE", 
                                        "ORCHIDEE", "SWAT", "CWatM", "MHM"))

# Plot call graphs -------------------------------------------------------------

tmp <- all_graphs[selected_models, on = .(model)] 
plot_all_risky_paths <- plot_grid(plotlist = tmp$plot_obj_additive, ncol = 2, align = "hv")
plot_all_risky_paths <- plot_grid(legend2, plot_all_risky_paths, ncol = 1, rel_heights = c(0.1, 0.9))


## ----same_additive_power_paths, dependson="ua_sa"---------------------------------------------

# HOW MANY OF THE TOP N PATHS ARE THE SAME IN THE ADDITIVE AND POWER-MEAN? ###

top_n <- 10
top10 <- paths_us_all[, .SD[order(-P_median)][1:top_n], .(model, risk_form)]

# Separate the two risk definitions --------------------------------------------

top_add <- unique(top10[risk_form == "additive",  .(model, path_id)]) 
top_pow <- unique(top10[risk_form == "power_mean", .(model, path_id)])

# Count overlaps per model -----------------------------------------------------

top_add[top_pow, on = .(model, path_id), .N, model] %>%
  na.omit() %>%
  .[selected_models, on = .(model)] %>%
  .[order(-N)] 



## ----risk_slope, dependson="run_analysis", fig.height=2.2, fig.width=3------------------------

# PLOT DISTRIBUTION OF RISK SLOPE AND GINI FOR THE TOP TEN PATHS ###############

plot_spread_risk <- paths_all[order(-path_risk_score), .SD[1:10], .(model, risk_form)] %>%
  .[selected_models, on = .(model)] %>%
  .[, risk_form:= fifelse(risk_form == "power_mean", "compensatory", risk_form)] %>%
  ggplot(., aes(reorder(model, risk_slope), risk_slope, fill = risk_form)) +
  scale_fill_manual(values = color_functional_forms, name = "") +
  geom_boxplot() +
  labs(x = "", y = expression(theta[1*k])) +
  coord_flip() +
  theme_AP() + 
  theme(legend.position = c(0.4, 0.7))

plot_spread_risk


## ----path_level_risk, dependson = "run_analysis", fig.height=2.5, fig.width=3-----------------

# PATH-LEVEL RISK ACCOUNTED FOR THE TOP 5% NODES ###############################

# To long format ---------------------------------------------------------------

paths_long <- paths_all[, .(node = unlist(path_nodes),
                            path_risk_score = path_risk_score,
                            gini_node_risk = gini_node_risk,
                            risk_slope = risk_slope,
                            risk_mean = risk_mean,
                            risk_sum  = risk_sum), 
                        .(model, path_id, risk_form)]

# Aggregate at function level --------------------------------------------------

node_from_paths <- paths_long[, .(n_paths = .N,
                                  mean_p_path = mean(path_risk_score, na.rm = TRUE),
                                  max_p_path = max(path_risk_score,  na.rm = TRUE),
                                  sum_p_path = sum(path_risk_score,  na.rm = TRUE),
                                  mean_gini = mean(gini_node_risk, na.rm = TRUE),
                                  mean_slope = mean(risk_slope, na.rm = TRUE),
                                  mean_risksum = mean(path_risk_score, na.rm = TRUE)),
                              .(model, node, risk_form)]

# Join with nodes --------------------------------------------------------------

node_summary <- merge(node_from_paths, nodes_all[, .(name, model, risk_score, risk_form)], 
                      by.x = c("model", "node", "risk_form"),
                      by.y = c("model", "name", "risk_form"), all.x = TRUE)

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
  .(model, risk_form)
  ]
}

# Run function -----------------------------------------------------------------

tmp <- top_share(0.05) 

table_top5 <- tmp %>%
  dcast(., model ~ risk_form, value.var = "share_risk_mass_topX") %>%
  .[order(-additive)] 

table_top5

table_top5 %>%
  .[selected_models, on = .(model)] %>%
  melt(., measure.vars = c("additive", "power_mean")) %>%
  ggplot(., aes(reorder(model, value), value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(0.5)) +
  coord_flip() +
  labs(y = "Path-level risk mass", x = "") +
  scale_fill_manual(values = color_functional_forms, name = "") +
  theme_AP() + 
  theme(legend.position = "top")

# Rerun name of the functions that dominate the top 5% of risk mass -----------

top_share_functions <- function(X = 0.05) {
  
  node_summary[
    !is.na(risk_mass) & risk_mass >= 0,
    {
      dt <- .SD[order(-risk_mass)]
      n_top <- max(1L, ceiling(.N * X))
      
      .(
        X = X,
        n_nodes = .N,
        n_top = n_top,
        share_risk_mass_topX = sum(dt$risk_mass[1:n_top]) / sum(dt$risk_mass),
        top_nodes = list(dt$node[1:n_top])
      )
    },
    by = .(model, risk_form)
  ]
}

tmp <- top_share_functions(0.05) 
tmp


## ----sensitivity_functions, dependson="ua_sa", fig.height=4.5, fig.width=2.5------------------

# SENSITIVITY OF SELECTED FUNCTIONS ############################################

functions_to_check <- c("main_grassland_management", "stomate_main", "model", 
                        "gwflow_simulate", "time_control")

plot_sa_functions <- indices_all[name %in% functions_to_check] %>%
  .[, name:= factor(name, levels = functions_to_check)] %>%
  .[sensitivity == "Si"] %>%
  ggplot(., aes(parameters, original, fill = model)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = expression(S[p])) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"), name = "") +
  scale_x_discrete(labels = c(a_raw = expression(alpha),
                              b_raw = expression(beta),
                              c_raw = expression(gamma), 
                              p = expression(p))) +
  facet_wrap(~name, ncol = 1) +
  theme_AP() + 
  scale_y_continuous(breaks = breaks_pretty(n = 2)) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 5),
        strip.text.y.right = element_text(size = 5), 
        strip.text.x.top = element_text(size = 4.5),
        strip.text.y = element_text(size = 7.4),
        strip.text = element_text(margin = margin(t = 1.5, b = 1.5)))

plot_sa_functions


## ----merged_figure, dependson=c("sensitivity_functions", "plot_all_callgraphs", "risk_slope"), fig.height=5.5, fig.width=6----

# MERGE FIGURES ################################################################

top <- plot_grid(plot_spread_risk, plot_sa_functions, ncol = 1, 
                 rel_heights = c(0.3, 0.7), labels = c("b", "c"))

plot_grid(plot_all_risky_paths, top, rel_widths = c(0.62, 0.38), 
          labels = c("a", ""))


## ----plot_tile_top_n, dependson="run_analysis", fig.height=3, fig.width=4---------------------

# WHICH PATHS IMPROVE IF A NODE'S RISK IS FIXED TO 0? ##########################

number_of_interest <- 30

all_graphs[!model == "VIC", 
           path_fix_heatmap:= lapply(paths_tbl_additive,
                                     function(x) {
                                       path_fix_heatmap(all_paths_out = x, 
                                                        n_nodes = number_of_interest, 
                                                        k_paths = number_of_interest)
  }
)]

# Add name of the model by reference -------------------------------------------

all_graphs[, path_fix_heatmap:= Map(
  function(res, m) {if (is.null(res) || is.null(res$plot)) return(res)
    res$plot <- res$plot + ggtitle(m)
    res
  },
  path_fix_heatmap,
  model
)]

# Plot -------------------------------------------------------------------------

all_graphs$path_fix_heatmap


## ----git_logs, fig.height=3, fig.width=4------------------------------------------------------

# METRICS AT THE FILE AND FUNCTION LEVEL #######################################

folder <- "./datasets/git_logs"

# Get names of files -----------------------------------------------------------

csv_files <- list.files(path = folder, pattern = "_2\\.csv$", full.names = TRUE)
csv_files

# Merge and flatten ------------------------------------------------------------

logs_dt <- lapply(csv_files, fread) %>%
  lapply(., function(x) x[, .(model, language, `function`, cyclomatic_complexity, number_changes, 
                              change_per_days, age_days)]) %>%
  rbindlist() %>%
  na.omit() %>%
  setnames(., "function", "node")

# Create and run function ------------------------------------------------------

spearman_fun <- function(dt) {
  
  cor.test(dt[, cyclomatic_complexity],
           dt[, number_changes], method = "spearman")$estimate[[1]]
}

rho_values <- logs_dt[, .(rho = spearman_fun(.SD)), model]

# Create dt to position labels ------------------------------------------------

pos_df <- logs_dt %>%
  filter(is.finite(cyclomatic_complexity), is.finite(number_changes)) %>%
  group_by(model) %>%
  summarise(x = min(cyclomatic_complexity, na.rm = TRUE),
            y = max(number_changes, na.rm = TRUE), .groups = "drop")

# join rho and make label text -------------------------------------------------

lab_df <- rho_values %>%
  left_join(pos_df, by = "model") %>%
  mutate(label = ifelse(is.na(rho), "rho==NA", paste0("rho==", round(rho, 2))))

# Plot -------------------------------------------------------------------------

plot_logs <- logs_dt %>%
  ggplot(aes(cyclomatic_complexity, number_changes, color = language)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = color_languages, name = "") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = expression(C), y = "number_changes") +
  facet_wrap(~model, scales = "free_y") +
  geom_text(data = lab_df %>% na.omit(),
            aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.2, size = 2, parse = TRUE) +
  coord_cartesian(clip = "off") +
  theme_AP() +
  theme(legend.position = "top")

plot_logs

# ARE REFACTORING ACTIONS MORE COMMON IN RISKY PATHS? ##########################

models_to_check <- unique(logs_dt$model)

# Extract name of functions in top ten risky paths per model -------------------

list_nodes <- paths_all %>%
  data.table() %>%
  .[order(-path_risk_score), .SD[1:10], model] %>%
  .[model %in% models_to_check] %>%
  na.omit() %>%
  .[, .(node = unlist(path_nodes)), model]

# keep only unique risk nodes
risk_nodes <- unique(list_nodes[, .(model, node)])
risk_nodes[, high_risk_path := TRUE]

# Filter out and keep only models that have log data ---------------------------

logs_dt <- logs_dt[model %in% models_to_check]

# left join into logs_dt -------------------------------------------------------

logs_dt[risk_nodes, high_risk_path:= i.high_risk_path,
        on = .(model, node)]

# Turn NAs into FALSE ----------------------------------------------------------

logs_dt[is.na(high_risk_path), high_risk_path := FALSE]

# We now quantify effect size using the common-language statistic A,
# defined as the probability that a randomly selected high-risk function 
# has higher churn than a randomly selected non-risk function. -----------------

wilcox_dt <- logs_dt[is.finite(number_changes),
                     {
                       x <- number_changes[high_risk_path]
                       y <- number_changes[!high_risk_path]
                       
                       n_x <- length(x)
                       n_y <- length(y)
                       
                       if (n_x >= 3L && n_y >= 3L) {
                         
                         wt <- wilcox.test(x, y, alternative = "greater", exact = FALSE)
                         
                         # Wilcoxon rank-sum statistic (sum of ranks for x)------------------------
                         W <- as.numeric(wt$statistic)
                         
                         # Convert to Mann–Whitney U ----------------------------------------------
                         U <- W - n_x * (n_x + 1) / 2
                         
                         # Common-language effect size (A statistic) ------------------------------
                         A <- U / (n_x * n_y)
                         
                         .(
                           p_value = wt$p.value,
                           A_common_language = A,
                           n_risk = n_x,
                           n_nonrisk = n_y,
                           median_risk = median(x),
                           median_nonrisk = median(y)
                         )
                         
                       } else {
                         
                         .(
                           p_value = NA_real_,
                           A_common_language = NA_real_,
                           n_risk = n_x,
                           n_nonrisk = n_y,
                           median_risk = if (n_x) median(x) else NA_real_,
                           median_nonrisk = if (n_y) median(y) else NA_real_
                         )
                       }
                     },
                     model
]


wilcox_dt



## ----print_top_paths, dependson="run_analysis", eval=FALSE, results="hide"--------------------
# 
# # FUNCTIONS TO SELECT THE TOP TEN RISKY PATHS PER MODEL AND
# # PRINT THEM OUT FOR LATEX #####################################################
# 
# # Top ten per model ------------------------------------------------------------
# 
# tmp <- paths_all[risk_function == "additive"] %>%
#   .[order(-path_risk_score), .SD[1:10], model] %>%
#   .[, .(model, path_str)] %>%
#   na.omit() %>%
#   split(., .$model)
# 
# 
# to_tex_list_fun(tmp)
# 
# # Random sample from bottom 20% per model --------------------------------------
# 
# set.seed(seed)
# 
# tmp2 <- paths_additive[, {n_bot <- ceiling(.N * 0.2)
# .SD[order(path_risk_score)] %>%
#     .[sample(seq_len(n_bot), size = min(5, n_bot))]
#   }, model] %>%
#   .[, .(model, path_str)] %>%
#   na.omit() %>%
#   split(., .$model)
# 
# to_tex_list_fun(tmp2)


## ----session_information----------------------------------------------------------------------

# SESSION INFORMATION ##########################################################

sessionInfo()

## Return the machine CPU ------------------------------------------------------

cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores -------------------------------------------------

cat("Num cores:   "); print(parallel::detectCores(logical = FALSE))

## Return number of threads ---------------------------------------------------

cat("Num threads: "); print(parallel::detectCores(logical = FALSE))

