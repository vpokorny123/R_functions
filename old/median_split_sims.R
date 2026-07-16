library(pwr)
#ok the question I'm really interested in is this:
# if there is no significant association between continuous variables
#is there any chance a binary split group difference is not type-1 error
iterations = 10000
n = 100
sig.level = .05
power = .5
# what about power
pwr_res<-pwr.r.test(n = n, sig.level = sig.level, power = power,
           alternative = 'two.sided')
pwr_res
r = pwr_res$r

all_cont_p = NULL
all_bin_p = NULL
set.seed(1)
for (i in seq(iterations)){
  sigma = matrix(c(1,r,
                   r,1),2,2, byrow = TRUE)
sim_data<-data.frame(MASS::mvrnorm(n=n, c(0,0),Sigma = sigma))
cont_res<-summary(lm(X1 ~X2, sim_data))
all_cont_p[i] = cont_res$coefficients[2,'Pr(>|t|)']

sim_data$X2_split <- as.numeric(sim_data$X2 > median(sim_data$X2))
bin_res<-t.test(X1 ~ X2_split, sim_data)
all_bin_p[i] <-bin_res$p.value
if (bin_res$p.value<.05 & cont_res$coefficients[2,'Pr(>|t|)'] >.05){
  data_of_interest <- sim_data
}
}

#how do we confirm power? power is the number of correct rejections
#complement of type-II


1-sum(all_cont_p>.05)/iterations
1-sum(all_bin_p>.05)/iterations

sum(all_bin_p <.05 & all_cont_p > .05)/iterations
sum(all_bin_p >.05 & all_cont_p < .05)/iterations
#sum(all_bin_p < all_cont_p)/iterations
#sum(all_bin_p > all_cont_p)/iterations

#scatterplot(data_of_interest$X2, data_of_interest$X1)
#scatterplot(data_of_interest$X2_split, data_of_interest$X1)





