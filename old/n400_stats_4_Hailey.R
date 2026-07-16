dat = read.csv("/Users/victorpokorny/Downloads/N400_Stim_MeanAmps.csv", header = TRUE)
hist(dat$Unrelated_m_Related)
t.test(Unrelated_m_Related ~ Group, 
       data = dat)

