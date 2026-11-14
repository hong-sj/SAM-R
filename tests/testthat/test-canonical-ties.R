# Tie-breaking that does not depend on rounding noise or index internals.
#
# The fixtures under golden/ties_* are shared verbatim with the Python
# package's test suite: both implementations must reproduce ties_expected.csv
# exactly. The data has only 72 distinct covariate profiles, so mathematically
# equal distances, losses and perimeters are everywhere -- before the canonical
# tie-break, R and Python (and two releases of either) could pick different
# matches from the same input. Expected rows are 0-based (Python) indices.

ties_fixture <- function() {
  list(
    data = utils::read.csv(testthat::test_path("golden", "ties_data.csv")),
    gps = as.matrix(utils::read.csv(testthat::test_path("golden", "ties_gps.csv"))),
    expected = utils::read.csv(testthat::test_path("golden", "ties_expected.csv"))
  )
}

for (method in c("sam", "3way")) {
  for (gps_space in c("raw", "logit")) {
    local({
      method <- method
      gps_space <- gps_space
      test_that(paste("tied data matches the shared fixture:", method, gps_space), {
        fx <- ties_fixture()
        search <- gps_candidate_search(fx$data, fx$gps, treatment_var = "T",
                                       anchor_level = "A", top_m = 5,
                                       gps_space = gps_space)
        result <- if (method == "sam") {
          sam_match(fx$data, search, X_vars = paste0("x", 0:4),
                    treatment_var = "T")
        } else {
          match_3way(fx$data, search, fx$gps, treatment_var = "T",
                     gps_space = gps_space, top_n = 5)
        }
        want <- fx$expected[fx$expected$method == method &
                              fx$expected$gps_space == gps_space, ]
        for (column in c("anchor", "B", "C")) {
          expect_identical(as.integer(result$matched[[column]]),
                           as.integer(want[[column]]) + 1L, info = column)
        }
      })
    })
  }
}

test_that("order_key absorbs rounding noise but not real differences", {
  expect_identical(order_key(1), order_key(1 + 4e-16))
  expect_identical(order_key(0.1 + 0.2), order_key(0.3))
  expect_lt(order_key(1), order_key(1 + 1e-9))
  # Round-half-even, as numpy.rint and Python's round().
  expect_identical(order_key(0.25e-10), 0)
  expect_identical(order_key(0.75e-10), 1)
})

test_that("sam_match is unaffected by a large covariate offset", {
  # Adding a constant to a covariate changes no Mahalanobis distance. The
  # expanded form |a|^2 + |b|^2 - 2ab loses about eps * offset^2 / var of
  # accuracy to cancellation -- here around 1e-4 -- while differencing
  # whitened rows is exact to rounding.
  set.seed(3L)
  n <- 900L
  data <- data.frame(x0 = stats::rnorm(n), x1 = stats::rnorm(n),
                     x2 = stats::rnorm(n))
  data$T <- sample(c("A", "B", "C"), n, replace = TRUE)
  raw <- matrix(stats::rgamma(3L * n, shape = 5), n, 3L)
  gps <- raw / rowSums(raw)
  colnames(gps) <- c("A", "B", "C")
  X_vars <- c("x0", "x1", "x2")

  search <- gps_candidate_search(data, gps, treatment_var = "T",
                                 anchor_level = "A")
  base <- sam_match(data, search, X_vars = X_vars, treatment_var = "T")$matched

  shifted_data <- data
  shifted_data$x0 <- shifted_data$x0 + 1e6
  shifted_search <- gps_candidate_search(shifted_data, gps, treatment_var = "T",
                                         anchor_level = "A")
  shifted <- sam_match(shifted_data, shifted_search, X_vars = X_vars,
                       treatment_var = "T")$matched

  expect_identical(shifted[c("anchor", "B", "C")], base[c("anchor", "B", "C")])
  expect_lt(max(abs(shifted$loss - base$loss)), 1e-8)
})

test_that("nearest prefers the lowest row among equidistant subjects", {
  coords <- rbind(c(0, 0), c(1, 1), c(1, 1), c(2, 2))
  group <- .active_group(coords)

  expect_identical(.ag_nearest(group, c(1, 1)), 2L)
  # Equidistant from (1.5, 1.5): rows 2, 3 and 4. The lowest wins.
  expect_identical(.ag_nearest(group, c(1.5, 1.5)), 2L)
  expect_identical(.nearest_untouched(group, coords, coords), c(1L, 2L, 2L, 4L))

  .ag_deactivate(group, 2L)
  expect_identical(.ag_nearest(group, c(1, 1)), 3L)

  .ag_deactivate(group, 3L)
  expect_identical(.ag_nearest(group, c(1.5, 1.5)), 4L)
})

test_that("within is strict and sorted after consumption", {
  coords <- rbind(c(0, 0), c(0.5, 0), c(1, 0), c(0, 0.5))
  group <- .active_group(coords)

  # The boundary point at distance exactly 1 is excluded.
  expect_identical(.ag_within(group, c(0, 0), 1), c(1L, 2L, 4L))

  .ag_deactivate(group, 1L)
  .ag_deactivate(group, 4L)  # triggers a rebuild

  expect_identical(.ag_within(group, c(0, 0), 1), 2L)
  expect_identical(.ag_within(group, c(0, 0), 1 + 1e-12), c(2L, 3L))
})
