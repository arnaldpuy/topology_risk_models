
# Load the static edge list and the runtime edge counts, both normalized -------
# format "two_file": static callgraph_edges.csv + runtime_callgraph_edges.csv
# format "combined": one folder with runtime_edge_overlap_diagnostics.csv
#                    (matched edges carry runtime_calls) and
#                    static_edges_not_hit_by_runtime.csv (runtime_calls = 0)

load_edges_fun <- function(static_path, runtime_path, format) {

  if (format == "two_file") {

    d  <- fread(static_path)
    el <- unique(data.table(from = norm_name_fun(d$caller),
                            to   = norm_name_fun(d$callee)))
    r  <- fread(runtime_path)
    r[, `:=`(cn = norm_name_fun(caller), ce = norm_name_fun(callee))]
    rt <- r[, .(runtime_calls = sum(runtime_calls, na.rm = TRUE)), by = .(cn, ce)]

  } else {

    diag <- fread(file.path(static_path, "runtime_edge_overlap_diagnostics.csv"))
    mt   <- diag[edge_gap_cause == "matched",
                 .(cn = norm_name_fun(caller_norm), ce = norm_name_fun(callee_norm),
                   runtime_calls)]
    nh   <- fread(file.path(static_path, "static_edges_not_hit_by_runtime.csv"))
    el   <- unique(rbind(mt[, .(from = cn, to = ce)],
                         nh[, .(from = norm_name_fun(caller_norm),
                                to   = norm_name_fun(callee_norm))]))
    rt   <- mt[, .(runtime_calls = sum(runtime_calls, na.rm = TRUE)), by = .(cn, ce)]
  }

  el <- el[from != "" & to != "" & !is.na(from) & !is.na(to)]
  setkey(rt, cn, ce)

  list(el = el, rt = rt)
}
