args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(arg_name, default = NULL) {
  match <- grep(paste0("^--", arg_name, "="), args, value = TRUE)
  if (length(match) == 0) return(default)
  sub(paste0("^--", arg_name, "="), "", match)
}