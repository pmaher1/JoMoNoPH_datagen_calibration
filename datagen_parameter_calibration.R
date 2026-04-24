library(tidyverse)

# ============================================================
# 1. Extract posterior summary statistics from Stan fit
# ============================================================

# previous inappropriate version using b2 = exp(b2); B2 = exp(B2) .* tau_aft;
# BP_objects <- readRDS("~/Documents/github/JoMoNoPH_datagen_calibration/BP_objects.rds")

# # current amended version using b2 = exp(b2) ./ tau_aft; B2 = exp(B2);
BP_objects <- readRDS("~/Documents/github/JoMoNoPH_datagen_calibration/BP_objects_DM.rds")

# Pull only the parameters needed for baseline hazard + RE
sum_re <- BP_objects$fit$summary(
  variables = c(
    "beta_long[1]", "beta_long[2]", "beta_long[3]",
    "sigma_long",
    "sd_1_long", "L_1_long",
    "beta_surv[1]", "alpha",
    "gamma[1]", "gamma[2]", "gamma[3]", "gamma[4]", "gamma[5]"
  )
)
sum_re

# ============================================================
# 2. Compute random-effects correlation (intercept–slope)
# ============================================================

# Extract Cholesky elements for L_1_long (2×2 lower-triangular)
L_mean <- sum_re$mean[grepl("^L_1_long", sum_re$variable)]

# Rebuild the Cholesky matrix
L <- matrix(0, 2, 2)
# L[lower.tri(L, diag = TRUE)] <- L_mean
L <- matrix(L_mean, nrow = 2, ncol = 2, byrow = FALSE) # this fix correct or not?

# Correlation matrix Ω = L Lᵀ
Omega <- L %*% t(L)

# Extract correlation between random intercept & slope
corr_intercept_slope <- Omega[1, 2]
corr_intercept_slope


# ============================================================
# 3. Extract Bernstein–polynomial γ parameters (baseline hazard)
# ============================================================

gamma <- sum_re %>%
  filter(str_detect(variable, "^gamma")) %>%
  arrange(variable) %>%
  pull(mean)

m <- length(gamma)   # number of Bernstein basis terms = degree


# Max observed event time (defines Bernstein domain scaling)
# t_max <- max(cPFS$time)
t_max <- 93.33881 # from accessible data set
tau   <- t_max        # rescale t → u ∈ (0,1)


# ============================================================
# 4. Define Bernstein basis functions (CDF + PDF components)
# ============================================================

# CDF basis used for cumulative hazard
bernstein_cdf <- function(u, m) {
  k <- 1:m
  pbeta(u, shape1 = k, shape2 = m - k + 1)
}

# PDF basis used for instantaneous hazard
bernstein_pdf <- function(u, m) {
  k <- 1:m
  dbeta(u, shape1 = k, shape2 = m - k + 1)
}


# ============================================================
# 5. Compute cumulative baseline hazard H0(t)
# ============================================================

grid <- tibble(
  t = seq(0, t_max, length.out = 500)
) %>%
  mutate(
    # Normalize time → u ∈ (0,1)
    u  = pmin(pmax(t / tau, 1e-6), 1 - 1e-6),
    # Each row is a vector of Bernstein basis CDF values
    B  = map(u, ~ bernstein_cdf(.x, m)),
    # Baseline cumulative hazard under the amended DM scaling:
    # H0(t) = Σ γ_k F_k(u), where u = t / τ
    H0 = map_dbl(B, ~ sum(gamma * .x))
    # Old scaling for BP_objects.rds:
    # H0(t) = τ Σ γ_k F_k(u)
    # H0 = map_dbl(B, ~ tau * sum(gamma * .x))
  )

# ---- Plot: Bernstein cumulative baseline hazard ----
ggplot(grid, aes(t, H0)) +
  geom_line(size = 1.2, color = "steelblue") +
  labs(
    x = "Time",
    y = "Cumulative baseline hazard H0(t)",
    title = "Estimated Bernstein Cumulative Baseline Hazard",
    subtitle = paste0("t ∈ [0, ", round(t_max, 1), "]")
  ) +
  theme_minimal(base_size = 14)


# ============================================================
# 6. Compute instantaneous baseline hazard h0(t)
# ============================================================

grid_h <- tibble(
  t = seq(0, t_max, length.out = 500)
) %>%
  mutate(
    u  = pmin(pmax(t / tau, 1e-6), 1 - 1e-6),
    f  = map(u, ~ bernstein_pdf(.x, m)),
    # Instantaneous hazard: h0(t) = dH0(t)/dt = (1 / τ) Σ γ_k f_k(u)
    h0 = map_dbl(f, ~ sum(gamma * .x) / tau)
    # Old scaling for BP_objects.rds:
    # h0(t) = Σ γ_k f_k(u)
    # h0 = map_dbl(f, ~ sum(gamma * .x))
  )

# ---- Plot: Bernstein baseline hazard ----
ggplot(grid_h, aes(t, h0)) +
  geom_line(size = 1.2, color = "firebrick") +
  labs(
    x = "Time",
    y = "Baseline hazard h0(t)",
    title = "Estimated Bernstein Instantaneous Baseline Hazard"
  ) +
  theme_minimal(base_size = 14)


# ============================================================
# 7. Fit log-logistic parametric baseline to Bernstein H0(t)
# ============================================================

# Log-logistic cumulative hazard:
#     H(t) = log(1 + (t/scale)^shape)
H0_loglogis <- function(t, shape, scale) {
  log(1 + (t / scale)^shape)
}

# Use only positive times (avoid log(0))
grid_fit <- grid %>% filter(t > 0)

# Objective: minimise squared error between Bernstein H0 and log-logistic H0
obj_llog <- function(par, t, H_target) {
  shape <- exp(par[1])        # enforce positivity
  scale <- exp(par[2])
  H_ll  <- H0_loglogis(t, shape, scale)
  sum((H_target - H_ll)^2)
}

# Heuristic starting values:
# - shape ≈ 1.5 gives unimodal hazard
# - scale ≈ t where H0 = log(2) (median-like)
H_half <- log(2)
t_half <- grid_fit$t[which.min(abs(grid_fit$H0 - H_half))]

start_shape <- log(1.5)
start_scale <- log(t_half)

# Fit parameters
fit_par <- optim(
  par      = c(start_shape, start_scale),
  fn       = obj_llog,
  t        = grid_fit$t,
  H_target = grid_fit$H0,
  method   = "BFGS"
)

shape_hat <- exp(fit_par$par[1])
scale_hat <- exp(fit_par$par[2])

shape_hat
scale_hat

# ============================================================
# 8. Overlay: Bernstein vs fitted log-logistic cumulative hazard
# ============================================================

grid <- grid %>%
  mutate(
    H_ll = H0_loglogis(t, shape_hat, scale_hat)
  )

ggplot(grid, aes(t)) +
  geom_line(aes(y = H0),   size = 1.2, color = "steelblue") +
  geom_line(aes(y = H_ll), size = 1.0, color = "darkred", linetype = "dashed") +
  labs(
    x = "Time",
    y = "Cumulative baseline hazard",
    title = "Bernstein vs Log-logistic Approximation of Baseline Hazard"
  ) +
  theme_minimal(base_size = 14)


# ============================================================
# 9. Fit Weibull parametric baseline to Bernstein H0(t)
#    Same raw squared-error approach as the log-logistic fit
# ============================================================

# Weibull cumulative hazard:
#     H(t) = (t/scale)^shape
H0_weibull <- function(t, shape, scale) {
  (t / scale)^shape
}

# Use only positive times (avoid log(0))
grid_fit <- grid %>% filter(t > 0)

# Objective: minimise squared error between Bernstein H0 and Weibull H0
obj_weib <- function(par, t, H_target) {
  shape <- exp(par[1])        # enforce positivity
  scale <- exp(par[2])
  H_weib <- H0_weibull(t, shape, scale)
  if (any(!is.finite(H_weib))) return(.Machine$double.xmax / 10)
  err <- sum((H_target - H_weib)^2)
  if (!is.finite(err)) return(.Machine$double.xmax / 10)
  err
}

# Heuristic starting values:
# - shape ≈ 1.5 gives a flexible increasing/decreasing hazard
# - scale ≈ t where H0 = 1, since Weibull H(scale) = 1
H_one <- 1
t_one <- grid_fit$t[which.min(abs(grid_fit$H0 - H_one))]

start_shape_weib <- log(1.5)
start_scale_weib <- log(t_one)

# Fit parameters
fit_par_weib <- optim(
  par      = c(start_shape_weib, start_scale_weib),
  fn       = obj_weib,
  t        = grid_fit$t,
  H_target = grid_fit$H0,
  method   = "L-BFGS-B",
  lower    = c(log(1e-3), log(min(grid_fit$t) * 1e-3)),
  upper    = c(log(50),   log(max(grid_fit$t) * 1e3))
)

shape_weib_plain_hat <- exp(fit_par_weib$par[1])
scale_weib_plain_hat <- exp(fit_par_weib$par[2])

shape_weib_plain_hat
scale_weib_plain_hat

# ---- Overlay: Bernstein vs fitted Weibull cumulative hazard ----
grid <- grid %>%
  mutate(
    H_weib_plain = H0_weibull(t, shape_weib_plain_hat, scale_weib_plain_hat)
  )

ggplot(grid, aes(t)) +
  geom_line(aes(y = H0),           size = 1.2, color = "steelblue") +
  geom_line(aes(y = H_weib_plain), size = 1.0, color = "darkred", linetype = "dashed") +
  labs(
    x = "Time",
    y = "Cumulative baseline hazard",
    title = "Bernstein vs Weibull Approximation of Baseline Hazard"
  ) +
  theme_minimal(base_size = 14)


# ============================================================
# 10. Fit Weibull parametric baseline to Bernstein H0(t) (stable)
#     Tail-mass rule: choose p_tail so a target fraction of weight
#     lies in the tail (t >= quantile(t, tail_q)), then weighted log-fit
# ============================================================

if (exists("p_tail")) rm(p_tail)

# Use only positive times and strictly positive hazards
grid_fit <- grid %>% dplyr::filter(t > 0, H0 > 0) %>% dplyr::arrange(t)

# ---- Starting values from Weibull log-linear identity ----
init_from_lm <- function(dat) {
  lm_init <- lm(log(H0) ~ log(t), data = dat)
  shape0  <- max(1e-3, unname(coef(lm_init)[2]))
  scale0  <- exp(-unname(coef(lm_init)[1]) / shape0)
  c(shape0 = shape0, scale0 = scale0)
}

# ---- Tail-mass rule to pick p_tail ----
# Choose p_tail so that:
#   sum_{t >= t_cut} w(t) / sum w(t) = target_mass
choose_p_tail_by_mass <- function(t, tail_q = 0.7, target_mass = 0.6) {
  t <- sort(t)
  tmax <- max(t)
  tcut <- as.numeric(quantile(t, tail_q))
  f <- function(p) {
    w <- (t / tmax)^p
    mass <- sum(w[t >= tcut]) / sum(w)
    mass - target_mass
  }
  # p=0 => mass ~ 1 - tail_q ; increasing p increases tail mass
  uniroot(f, interval = c(0, 50))$root
}

# ---- Set your tail definition + desired emphasis ----
tail_q      <- 0.70   # tail starts at 70th percentile of time
target_mass <- 0.60   # e.g., 60% of total weight should lie in the tail

p_tail <- choose_p_tail_by_mass(grid_fit$t, tail_q = tail_q, target_mass = target_mass)
p_tail

# ---- Fit weighted Weibull on log cumulative hazard ----
fit_weib_weighted <- function(dat, p_tail) {
  init <- init_from_lm(dat)
  shape0 <- init["shape0"]; scale0 <- init["scale0"]
  w <- (dat$t / max(dat$t))^p_tail
  obj_weib_log_w <- function(par, t, H_target, w) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    eps <- 1e-12
    logH_target <- log(pmax(H_target, eps))
    logH_weib   <- log(pmax(H0_weibull(t, shape, scale), eps))
    sum(w * (logH_target - logH_weib)^2)
  }
  fit <- optim(
    par      = c(log(shape0), log(scale0)),
    fn       = obj_weib_log_w,
    t        = dat$t,
    H_target = dat$H0,
    w        = w,
    method   = "L-BFGS-B",
    lower    = c(log(1e-3), log(min(dat$t) * 1e-3)),
    upper    = c(log(50),   log(max(dat$t) * 1e3))
  )
  list(shape = exp(fit$par[1]), scale = exp(fit$par[2]))
}

# ---- Fit on all data using chosen p_tail ----
pars_final <- fit_weib_weighted(grid_fit, p_tail)
shape_weib_hat <- pars_final$shape
scale_weib_hat <- pars_final$scale

shape_weib_hat
scale_weib_hat

# ---- Overlay: Bernstein vs tuned weighted-Weibull cumulative hazard ----
grid <- grid %>%
  mutate(H_weib = H0_weibull(pmax(t, 1e-12), shape_weib_hat, scale_weib_hat))

ggplot(grid, aes(t)) +
  geom_line(aes(y = H0),     size = 1.2, color = "steelblue") +
  geom_line(aes(y = H_weib), size = 1.0, color = "darkred", linetype = "dashed") +
  labs(
    x = "Time",
    y = "Cumulative baseline hazard",
    title = "Bernstein vs Weighted-Weibull Approximation of Baseline Hazard",
    subtitle = paste0(
      "Tail-mass rule: tail_q=", tail_q,
      ", target_mass=", target_mass,
      " ⇒ p_tail=", signif(p_tail, 4),
      "; shape=", signif(shape_weib_hat, 4),
      ", scale=", signif(scale_weib_hat, 4)
    )
  ) +
  theme_minimal(base_size = 14) 


# ============================================================
# 11. Overlay all cumulative baseline hazard estimates
# ============================================================

grid_overlay <- grid %>%
  mutate(
    `Bernstein H0` = H0,
    `Log-logistic` = H0_loglogis(t, shape_hat, scale_hat),
    `Weibull, raw SSE` = H0_weibull(t, shape_weib_plain_hat, scale_weib_plain_hat),
    `Weibull, weighted log` = H0_weibull(pmax(t, 1e-12), shape_weib_hat, scale_weib_hat)
  ) %>%
  select(
    t,
    `Bernstein H0`,
    `Log-logistic`,
    `Weibull, raw SSE`,
    `Weibull, weighted log`
  ) %>%
  pivot_longer(
    cols = -t,
    names_to = "curve",
    values_to = "H"
  )

ggplot(grid_overlay, aes(x = t, y = H, color = curve, linetype = curve)) +
  geom_line(aes(linewidth = curve)) +
  scale_color_manual(
    values = c(
      `Bernstein H0` = "black",
      `Log-logistic` = "steelblue",
      `Weibull, raw SSE` = "darkred",
      `Weibull, weighted log` = "darkgreen"
    )
  ) +
  scale_linetype_manual(
    values = c(
      `Bernstein H0` = "solid",
      `Log-logistic` = "dashed",
      `Weibull, raw SSE` = "dotdash",
      `Weibull, weighted log` = "twodash"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      `Bernstein H0` = 1.2,
      `Log-logistic` = 1.0,
      `Weibull, raw SSE` = 1.0,
      `Weibull, weighted log` = 1.0
    )
  ) +
  labs(
    x = "Time",
    y = "Cumulative baseline hazard",
    color = NULL,
    linetype = NULL,
    linewidth = NULL,
    title = "Bernstein H0 vs Parametric Baseline Approximations",
    subtitle = "All fitted curves use the same Bernstein H0 target from grid$H0"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")


# ============================================================
# 12. Overlay all baseline hazard estimates
# ============================================================

# Log-logistic baseline hazard:
#     h(t) = d/dt log(1 + (t/scale)^shape)
h0_loglogis <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1) /
    (1 + (t / scale)^shape)
}

# Weibull baseline hazard:
#     h(t) = d/dt (t/scale)^shape
h0_weibull <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1)
}

grid_h_overlay <- grid_h %>%
  filter(t > 0) %>%
  mutate(
    `Bernstein h0` = h0,
    `Log-logistic` = h0_loglogis(t, shape_hat, scale_hat),
    `Weibull, raw SSE` = h0_weibull(t, shape_weib_plain_hat, scale_weib_plain_hat),
    `Weibull, weighted log` = h0_weibull(t, shape_weib_hat, scale_weib_hat)
  ) %>%
  select(
    t,
    `Bernstein h0`,
    `Log-logistic`,
    `Weibull, raw SSE`,
    `Weibull, weighted log`
  ) %>%
  pivot_longer(
    cols = -t,
    names_to = "curve",
    values_to = "h"
  )

ggplot(grid_h_overlay, aes(x = t, y = h, color = curve, linetype = curve)) +
  geom_line(aes(linewidth = curve)) +
  scale_color_manual(
    values = c(
      `Bernstein h0` = "black",
      `Log-logistic` = "steelblue",
      `Weibull, raw SSE` = "darkred",
      `Weibull, weighted log` = "darkgreen"
    )
  ) +
  scale_linetype_manual(
    values = c(
      `Bernstein h0` = "solid",
      `Log-logistic` = "dashed",
      `Weibull, raw SSE` = "dotdash",
      `Weibull, weighted log` = "twodash"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      `Bernstein h0` = 1.2,
      `Log-logistic` = 1.0,
      `Weibull, raw SSE` = 1.0,
      `Weibull, weighted log` = 1.0
    )
  ) +
  labs(
    x = "Time",
    y = "Baseline hazard",
    color = NULL,
    linetype = NULL,
    linewidth = NULL,
    title = "Bernstein h0 vs Parametric Baseline Hazard Approximations",
    subtitle = "Parametric h0 curves are analytic derivatives of the fitted cumulative hazards"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
