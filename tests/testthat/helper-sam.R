################################################################################
### Shared test fixtures
##
## Building the GPS runs `nnet::multinom` over the bundled data, which is slow
## enough that repeating it in every test file dominates the suite. Fixtures
## are therefore memoised for the lifetime of the session.
################################################################################

sam_covariates <- function(data) {
  setdiff(names(data), c("synthetic_id", "treatment", "mortality_28d"))
}

sam_anchor <- function(data) {
  unique(as.character(data$treatment))[1]
}

#' Standard pipeline fixture
#'
#' Returns a list with `data`, `X_vars`, `anchor`, `groups`, `gps`, `search`
#' and `match` for one bundled dataset.
sam_fixture <- local({
  cache <- new.env(parent = emptyenv())

  function(dataset = "sample_4group", gps_space = "logit", top_m = 10L) {
    key <- paste(dataset, gps_space, top_m, sep = "/")
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }

    data <- get(utils::data(list = dataset, package = "SAMatch",
                            envir = environment()))
    X_vars <- sam_covariates(data)
    anchor <- sam_anchor(data)

    fit <- estimate_gps_multinom(
      data,
      X_vars = X_vars,
      treatment_var = "treatment",
      anchor_level = anchor
    )
    search <- gps_candidate_search(
      data, fit$gps,
      treatment_var = "treatment",
      anchor_level = anchor,
      top_m = top_m,
      gps_space = gps_space
    )
    match_result <- sam_match(
      data, search,
      X_vars = X_vars,
      treatment_var = "treatment"
    )

    out <- list(
      data = data, X_vars = X_vars, anchor = anchor,
      groups = search$groups, gps = fit$gps,
      search = search, match = match_result
    )
    cache[[key]] <- out
    out
  }
})

#' A small, fully controlled dataset
#'
#' Well-conditioned covariates and a balanced treatment vector, for tests that
#' need a pipeline to run but do not care about the bundled data.
sam_small_data <- function(n = 120L, n_groups = 3L, seed = 20260814L) {
  set.seed(seed)
  out <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    x3 = stats::rnorm(n)
  )
  out$T <- rep(LETTERS[seq_len(n_groups)], length.out = n)
  out$y <- stats::rbinom(n, 1, 0.3)
  out
}
