data_generating_process<-rbinom(n = 1000, size = 1, prob = .5) 
plot(data_generating_process)
random_noise<- rnorm(n = 1000, sd = .5 )
plot(random_noise)
observed_signal<- data_generating_process + random_noise
plot(observed_signal)

#now split by .5

recovered_vals<-as.numeric(observed_signal<=.5)

scatterplot(data_generating_process, recovered_vals, )
