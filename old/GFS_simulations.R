set.seed(1)
library(MASS)
#correlation cbetween some criterion variable and GFS is .3
sd = 3 # make a little bigger due to range restriction
cor = .2 
iterations = 10000

reg_p_all <- NULL
max10_p_all <- NULL
reg_r_all <- NULL
max10_r_all <- NULL
for (j in seq(iterations)){
sigma<-matrix(c(sd^2,cor*(sd^2),cor*(sd^2),sd^2),2,2)
mu = c(0,7.8)
df<-data.frame(mvrnorm(mu = mu,Sigma = sigma, n = 100))
df$criterion<-df$X1
df$gfs<-df$X2

#now truncate to 10
df$gfs_max_10<-df$gfs
df$gfs_max_10[df$gfs>10] = 10

reg_p_all[j]<-cor.test(df$criterion, df$gfs)$p.value
max10_p_all[j]<-cor.test(df$criterion, df$gfs_max_10)$p.value
reg_r_all[j]<-cor.test(df$criterion, df$gfs)$estimate
max10_r_all[j]<-cor.test(df$criterion, df$gfs_max_10)$estimate

}

sum(reg_p_all<.05)/iterations
sum(max10_p_all<.05)/iterations

mean(reg_r_all)
mean(max10_r_all)
 
sd(reg_r_all)
sd(max10_r_all)

t.test(reg_r_all,max10_r_all)
