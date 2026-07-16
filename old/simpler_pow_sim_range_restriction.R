library(ggplot2)
set.seed(1)

# ---- Function to simulate unrestricted data ----
sim_full <- function(n = 5000, beta = 0.4, sd_x = 1, sd_e = 1) {
  x <- rnorm(n, 0, sd_x)
  y <- beta * x + rnorm(n, 0, sd_e)
  data.frame(x, y)
}

# ---- Function to apply range restriction to X ----
restrict_x <- function(dat, proportion_kept = 0.5) {
  # keep middle proportion_kept of the distribution
  lower <- quantile(dat$x, (1 - proportion_kept) / 2)
  upper <- quantile(dat$x, 1 - (1 - proportion_kept) / 2)
  dat[dat$x >= lower & dat$x <= upper, ]
}

# ---- Run simulation across several restriction levels ----
proportions <- c(1.0, 0.75, 0.50, 0.30, 0.10)  # percent of X retained

results <- lapply(proportions, function(p) {
  full <- sim_full()
  
  restricted <- restrict_x(full, proportion_kept = p)
  
  fit <- lm(y ~ x, data = restricted)
  
  data.frame(
    proportion = p,
    sd_x = sd(restricted$x),
    cor_xy = cor(restricted$x, restricted$y),
    slope = coef(fit)[2],
    n = nrow(restricted)
  )
})

results <- do.call(rbind, results)
print(results)

# ---- Plot attenuation of cor and slope ----
library(reshape2)

df_long <- melt(results,
                id.vars = "proportion",
                measure.vars = c("cor_xy", "slope"))

ggplot(df_long, aes(proportion, value)) +
  geom_line(size = 1) +
  geom_point() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Proportion of X retained",
       y = "Value",
       title = "Effect of range restriction of X on correlation and regression slope")
