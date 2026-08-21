
# Normalize function names: lowercase, trim, strip leading class qualifier -----

norm_name_fun <- function(x) sub(".*\\.", "", tolower(trimws(as.character(x))))
