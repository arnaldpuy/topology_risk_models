# 00b_extract_ua_summary.R
# One-off extractor: reduces the ~1.7 GB full_ua_df.csv to a slim
# full_ua_summary.csv (~tens of MB) that the dynamic-validation pipeline
# reads. Drops the per-path P_k_vec column after using it to compute the
# exact P_k_sd (and CV_k) for each path.
#
# Run from repository root *once* after the static uncertainty analysis has
# been re-run, or whenever full_ua_df.csv changes. The slim file is small
# enough to commit to the repo, so the validation pipeline is reproducible
# without the full 1.7 GB raw file.

library(data.table)

in_file  <- "full_ua_df.csv"
out_file <- "full_ua_summary.csv"

if (!file.exists(in_file)) stop("Missing ", in_file, " at ", getwd())

cat("[extract_ua] reading ", in_file, " (",
    format(file.size(in_file) / 1e9, digits = 3), " GB) ...\n", sep = "")

t0 <- Sys.time()
ua <- fread(in_file)
cat("[extract_ua] read ", nrow(ua), " rows x ", ncol(ua), " cols in ",
    format(round(difftime(Sys.time(), t0, units = "secs"), 1)), "\n", sep = "")

# Exact sd from the ensemble vector
cat("[extract_ua] parsing P_k_vec and computing P_k_sd ...\n")
t1 <- Sys.time()
ua[, P_k_sd := vapply(
  strsplit(P_k_vec, "|", fixed = TRUE),
  function(x) sd(as.numeric(x), na.rm = TRUE),
  numeric(1)
)]
ua[, CV_k := P_k_sd / pmax(P_k_mean, 1e-12)]
cat("[extract_ua] done in ",
    format(round(difftime(Sys.time(), t1, units = "secs"), 1)), "\n", sep = "")

slim <- ua[, .(model, path_id, path_nodes, path_str, hops,
               P_k_min, P_k_mean, P_k_max,
               P_k_q025, P_k_q50, P_k_q975,
               P_k_sd, CV_k)]

fwrite(slim, out_file)
cat("[extract_ua] wrote ", out_file, " (",
    format(file.size(out_file) / 1e6, digits = 3), " MB).\n", sep = "")
