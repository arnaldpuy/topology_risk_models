
# FUNCTIONS TO SELECT THE TOP TEN RISKY PATHS PER MODEL AND
# PRINT THEM OUT FOR LATEX #####################################################

escape_latex <- function(x) {
  
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%#$&_{}])", "\\\\\\1", x)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

to_tex_list_fun <- function(lst) {
  cat("\\subsection{Top risky paths}\n\n")
  for (nm in names(lst)) {
    dt <- lst[[nm]]
    if (!is.data.table(dt) || nrow(dt) == 0) next
    
    # split "CTSM.fortran" -> "CTSM (fortran)"
    title <- sub("\\.(.+)$", " (\\1)", nm)
    
    paths <- dt[["path_str"]]
    paths <- escape_latex(paths)
    paths <- gsub("→", "\\\\textrightarrow{}", paths)  
    
    cat(sprintf("\\subsection{%s}\n", title))
    cat("\\begin{enumerate}[leftmargin=*, itemsep=0.25em]\n")
    cat(paste0("\\item {\\small\\ttfamily ", paths, "}\n", collapse = ""))
    cat("\\end{enumerate}\n\n")
  }
}
