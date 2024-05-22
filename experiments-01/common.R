# Check if lgpmodel has correct format
check_model_compatibility <- function(model) {
  stopifnot(is(model, "lgpmodel"))
  ver <- model@info$lgpr_version
  if (ver < "1.1.4" || ver >= "1.2.0") {
    stop(
      "model was created with lgpr ", ver, ", but should be created with ",
      "lgpr version number at least 1.1.4 and less than 1.2.0!"
    )
  }
  TRUE
}

# Check lgpr package version
check_lgpr_version <- function() {
  ver <- packageVersion("lgpr")
  if (ver < "1.1.4" || ver >= "1.2.0") {
    stop(
      "using lgpr version ", ver, ", but ",
      "lgpr version number at least 1.1.4 and less than 1.2.0 is needed!"
    )
  }
  TRUE
}

# Startup for all experiments
startup <- function(array_idx = NULL, create_dir = TRUE) {
  library(lgpr)
  check_lgpr_version()
  library(ggplot2)
  library(ggpubr)
  library(posterior)
  library(rstan)
  library(RColorBrewer)
  library(methods)
  rstan::rstan_options(javascript = FALSE)
  rstan::rstan_options(auto_write = TRUE)
  library(cmdstanr)
  if (is.null(array_idx)) {
    outdir <- "results"
  } else {
    outdir <- file.path("results", paste0("repl_", array_idx))
    if (create_dir) {
      if (!dir.exists("results")) dir.create("results")
    }
  }
  if (create_dir) {
    if (!dir.exists(outdir)) dir.create(outdir)
    cat(" * results will be saved to:", outdir, "\n")
  } else {
    outdir <- NULL
  }
  # Source R files and set stan dir
  r_dir <- normalizePath("R")
  stan_dir <- normalizePath("stan")
  options(stan_dir = stan_dir)
  cat("* set option stan_dir =", stan_dir, "\n")
  cat("* sourcing all R files in", r_dir, "\n")
  for (f in dir(r_dir)) {
    source(file.path(r_dir, f))
  }
  return(outdir)
}



# Expose all Stan functions without creating a complete Stan model
expose_stanfuns <- function() {
  stan_dir <- getOption("stan_dir")
  funs_dir <- file.path(stan_dir, "chunks")
  FILES <- list.files(path = funs_dir, pattern = "functions")
  FILES <- file.path(funs_dir, FILES)
  cat(" * exposing functions in following files:\n")
  print(FILES)

  # Create Stan model containing only a functions block with all the functions
  two_spaces <- "  "
  f_list <- lapply(FILES, FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
  header <- "// Automatically generated\n\n"
  model_code <- paste0(header, model_code, "\n")

  # Write Stan code to file and expose
  tfile <- tempfile()
  fn <- paste0(tfile, ".stan")
  cat(model_code, file = fn)
  rstan::expose_stan_functions(fn, verbose = TRUE)
  file.remove(fn)
}

# Matrix rows to a list
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- as.numeric(x[i, ])
  }
  return(L)
}

# Matrix rows to a list of lists
matrix_to_list_of_lists <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- as.list(x[i, ])
  }
  return(L)
}

# Standardize variable to zero mean and unit variance
normalize_var <- function(x) (x - mean(x)) / stats::sd(x)
