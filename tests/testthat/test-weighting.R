################################################################################
### Multi-arm balancing weights, trimming and effective sample size (#9)
################################################################################

separated_data <- function() {
  # A covariate that predicts treatment almost perfectly drives GPS columns
  # toward zero, which is exactly the near-separation an unregularised
  # multinomial fit cannot protect against. The four reassigned subjects then
  # hold a treatment their own covariates make nearly impossible, so it is
  # their own-group score -- the divisor -- that approaches zero.
  set.seed(20260814L)
  n <- 150L
  data <- data.frame(x1 = c(rep(-8, 50), rep(0, 50), rep(8, 50)))
  data$x2 <- stats::rnorm(n)
  data$T <- rep(c("A", "B", "C"), each = 50)
  data$T[c(1, 2)] <- "C"
  data$T[c(149, 150)] <- "A"
  data
}

test_that("the bundled data trims nothing, so results are unchanged", {
  fixture <- sam_fixture()

  for (method in c("iptw", "overlap", "matching")) {
    weighting <- compute_balancing_weights(
      fixture$data, method = method, X_vars = fixture$X_vars,
      treatment_var = "treatment", anchor_level = fixture$anchor
    )
    expect_equal(weighting$n_trimmed, 0L)

    untrimmed <- compute_balancing_weights(
      fixture$data, method = method, X_vars = fixture$X_vars,
      treatment_var = "treatment", anchor_level = fixture$anchor, trim = 0
    )
    expect_equal(weighting$weights, untrimmed$weights)
  }
})

test_that("trimming bounds an IPTW weight that would otherwise run away", {
  data <- separated_data()
  X_vars <- c("x1", "x2")

  untrimmed <- compute_balancing_weights(
    data, method = "iptw", X_vars = X_vars, treatment_var = "T",
    anchor_level = "A", trim = 0
  )
  trimmed <- compute_balancing_weights(
    data, method = "iptw", X_vars = X_vars, treatment_var = "T",
    anchor_level = "A", trim = 1e-3
  )

  expect_gt(trimmed$n_trimmed, 0)
  expect_equal(untrimmed$n_trimmed, 0L)
  expect_true(all(is.finite(trimmed$weights)))

  # A single subject holding a treatment their covariates make near-impossible
  # produced a weight an order of magnitude beyond anything else in the sample.
  expect_gt(max(untrimmed$weights), 1000)
  expect_lt(max(trimmed$weights), max(untrimmed$weights) / 10)
})

test_that("overlap and matching weights were already bounded", {
  # The file header says these two do not need stabilization because they are
  # bounded by construction; trimming should therefore barely move them, which
  # is why the runaway above is specific to IPTW.
  data <- separated_data()
  X_vars <- c("x1", "x2")

  for (method in c("overlap", "matching")) {
    untrimmed <- compute_balancing_weights(
      data, method = method, X_vars = X_vars, treatment_var = "T",
      anchor_level = "A", trim = 0
    )
    expect_lte(max(untrimmed$weights), 1 + 1e-8)
  }
})

test_that("n_trimmed is a positivity diagnostic, not just a counter", {
  data <- separated_data()
  X_vars <- c("x1", "x2")

  gps <- estimate_gps_multinom(data, X_vars = X_vars, treatment_var = "T",
                               anchor_level = "A")$gps
  expected <- sum(apply(gps < 1e-3, 1, any))

  weighting <- compute_balancing_weights(
    data, method = "overlap", gps = gps, treatment_var = "T",
    anchor_level = "A"
  )
  expect_equal(weighting$n_trimmed, as.integer(expected))
})

test_that("the returned gps is the untrimmed one", {
  data <- separated_data()
  gps <- estimate_gps_multinom(data, X_vars = c("x1", "x2"),
                               treatment_var = "T", anchor_level = "A")$gps

  weighting <- compute_balancing_weights(
    data, method = "matching", gps = gps, treatment_var = "T",
    anchor_level = "A", trim = 0.05
  )
  expect_equal(weighting$gps, gps)
  expect_gt(weighting$n_trimmed, 0)
})

test_that("trim is validated", {
  fixture <- sam_fixture()

  weights_with <- function(trim) {
    compute_balancing_weights(
      fixture$data, method = "iptw", X_vars = fixture$X_vars,
      treatment_var = "treatment", anchor_level = fixture$anchor, trim = trim
    )
  }

  expect_error(weights_with(-0.1), "\\[0, 1\\)")
  expect_error(weights_with(1), "\\[0, 1\\)")
  expect_error(weights_with(c(0.1, 0.2)), "\\[0, 1\\)")
  expect_error(weights_with(NA_real_), "\\[0, 1\\)")
  expect_silent(weights_with(0))
})

test_that("each weighting method produces the shape it claims", {
  fixture <- sam_fixture()

  for (method in c("iptw", "overlap", "matching")) {
    weighting <- compute_balancing_weights(
      fixture$data, method = method, X_vars = fixture$X_vars,
      treatment_var = "treatment", anchor_level = fixture$anchor
    )
    expect_length(weighting$weights, nrow(fixture$data))
    expect_length(weighting$h, nrow(fixture$data))
    expect_identical(weighting$method, method)
    expect_true(all(weighting$weights > 0))
    expect_true(all(is.finite(weighting$weights)))
  }
})

test_that("the tilting functions match their definitions", {
  fixture <- sam_fixture()
  gps <- fixture$gps

  overlap <- compute_balancing_weights(
    fixture$data, method = "overlap", gps = gps, treatment_var = "treatment",
    anchor_level = fixture$anchor, trim = 0
  )
  expect_equal(overlap$h, as.numeric(1 / rowSums(1 / gps)))

  matching <- compute_balancing_weights(
    fixture$data, method = "matching", gps = gps, treatment_var = "treatment",
    anchor_level = fixture$anchor, trim = 0
  )
  expect_equal(matching$h, as.numeric(apply(gps, 1, min)))

  iptw <- compute_balancing_weights(
    fixture$data, method = "iptw", gps = gps, treatment_var = "treatment",
    anchor_level = fixture$anchor, trim = 0, stabilize = FALSE
  )
  own <- gps[cbind(seq_len(nrow(gps)),
                   match(treatment_labels(fixture$data, "treatment"),
                         colnames(gps)))]
  expect_equal(iptw$weights, as.numeric(1 / own))
})

test_that("stabilization multiplies by the marginal prevalence", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")

  plain <- compute_balancing_weights(
    fixture$data, method = "iptw", gps = fixture$gps,
    treatment_var = "treatment", anchor_level = fixture$anchor,
    stabilize = FALSE, trim = 0
  )
  stabilized <- compute_balancing_weights(
    fixture$data, method = "iptw", gps = fixture$gps,
    treatment_var = "treatment", anchor_level = fixture$anchor,
    stabilize = TRUE, trim = 0
  )

  prevalence <- as.numeric(prop.table(table(labels))[labels])
  expect_equal(stabilized$weights, plain$weights * prevalence)

  # Ignored for the bounded methods.
  for (method in c("overlap", "matching")) {
    with_stabilize <- compute_balancing_weights(
      fixture$data, method = method, gps = fixture$gps,
      treatment_var = "treatment", anchor_level = fixture$anchor,
      stabilize = TRUE
    )
    without <- compute_balancing_weights(
      fixture$data, method = method, gps = fixture$gps,
      treatment_var = "treatment", anchor_level = fixture$anchor,
      stabilize = FALSE
    )
    expect_equal(with_stabilize$weights, without$weights)
  }
})

test_that("compute_effective_sample_size matches Kish's formula", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")
  weights <- compute_balancing_weights(
    fixture$data, method = "overlap", gps = fixture$gps,
    treatment_var = "treatment", anchor_level = fixture$anchor
  )$weights

  ess <- compute_effective_sample_size(fixture$data, weights,
                                       treatment_var = "treatment")

  expect_setequal(ess$group, unique(labels))
  for (i in seq_len(nrow(ess))) {
    w <- weights[labels == ess$group[i]]
    expect_equal(ess$n[i], length(w))
    expect_equal(ess$ess[i], sum(w)^2 / sum(w^2))
    expect_equal(ess$max_over_mean[i], max(w) / mean(w))
  }
  # The effective sample size can never exceed the actual one.
  expect_true(all(ess$ess <= ess$n + 1e-8))
})

test_that("evaluate_comparator_weighting surfaces n_trimmed", {
  data <- separated_data()

  report <- evaluate_comparator_weighting(
    data, method = "iptw", X_vars = c("x1", "x2"),
    treatment_var = "T", anchor_level = "A"
  )

  expect_equal(report$n_trimmed, report$weights$n_trimmed)
  expect_gt(report$n_trimmed, 0)
  expect_true(all(is.finite(report$ess$ess)))
  expect_true(all(is.finite(report$balance$summary$mean_abs_smd)))
})

test_that("compute_weighted_balance validates its weights", {
  fixture <- sam_fixture()

  balance_with <- function(weights) {
    compute_weighted_balance(fixture$data, weights, X_vars = fixture$X_vars,
                             treatment_var = "treatment",
                             anchor_level = fixture$anchor)
  }

  expect_error(balance_with(rep(1, 3)), "same number of observations")
  expect_error(balance_with(c(NA_real_, rep(1, nrow(fixture$data) - 1))),
               "must all be finite")
  expect_error(balance_with(c(Inf, rep(1, nrow(fixture$data) - 1))),
               "must all be finite")
  expect_error(balance_with(c(-1, rep(1, nrow(fixture$data) - 1))),
               "non-negative")
})

test_that("a treatment group with zero total weight is named, not NaN", {
  fixture <- sam_fixture()
  labels <- treatment_labels(fixture$data, "treatment")
  starved <- fixture$groups[1]

  weights <- rep(1, nrow(fixture$data))
  weights[labels == starved] <- 0

  # Previously: every weighted mean for that group became 0/0 and the table
  # came back full of NaN.
  expect_error(
    compute_weighted_balance(fixture$data, weights, X_vars = fixture$X_vars,
                             treatment_var = "treatment",
                             anchor_level = fixture$anchor),
    starved,
    fixed = TRUE
  )
})
