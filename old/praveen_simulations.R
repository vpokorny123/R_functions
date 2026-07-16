library(ggplot2)
set.seed(6)
timepoints = 160 
subjects = 100 #assuming 300 subjects per group
mean = .5
sd = .1
n = timepoints * subjects

#next create two random matrices (mean of .5 and sd of .25) with dimensions of subject by timepoint
group_1_time_series <- matrix(rnorm(n = n, mean = mean, sd = sd),subjects, timepoints)
group_2_time_series <- matrix(rnorm(n = n, mean = mean, sd = sd),subjects, timepoints)

matplot(t(group_1_time_series), type = "l", lty = 1,
        xlab = "Time", ylab = "Value", main = "Group 1 Timeseries")

matplot(t(group_2_time_series), type = "l", lty = 1,
        xlab = "Time", ylab = "Value", main = "Group 2 Timeseries")

#randomly select two belief parameter values for each group
belief_1 <- rnorm(1, mean = mean, sd = sd)
belief_2 <- rnorm(1, mean = mean, sd = sd)

#average across subjects and scale by belief parameter
group_1_avg <- colMeans(group_1_time_series) * belief_1
group_2_avg <- colMeans(group_2_time_series) * belief_2

#plot timeseries
plot(x = seq(timepoints),y = group_1_avg, 'l', ylim = c(0,1), ylab = '')
lines(x = seq(timepoints),y = group_2_avg, 'l')


# Combine into a data frame
df <- data.frame(
  value = c(group_1_avg, group_2_avg),
  group = factor(rep(c("Group 1", "Group 2"), each = length(group_1_avg)))
)

# Plot
ggplot(df, aes(x = group, y = value)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.15)) +
  # or geom_point(), geom_violin(), etc.
  theme_minimal()
#now test for group difference
pub_ready_stats(t.test(group_1_avg, group_2_avg ,var.equal = TRUE))
