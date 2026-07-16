library(boot)
f_boot <- function(dat, idx) coef(lm(rating ~ network + ill_ness_group, data = dat[idx, ]))
b <- boot(df, f_boot, R = 2000)
plot(b)
boot.ci(b, type = "perc", index = 2)  # CI for x1
