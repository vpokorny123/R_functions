df<-read.csv('/Users/victorpokorny/Downloads/all_tasks - Sheet1.csv')

#you will likely have to run 'install.packages('psych')'
library(psych)

#select only the variables we want to include in the factor analysis
data_matrix <- df[,c('CC.mean','SGI.B.Total','VA','CS','HA','VVIQ')]

#drop percentage from CS
data_matrix$CS<- substr(data_matrix$CS,1,nchar(data_matrix$CS)-1)

#convert all variables to numeric (instead of character)
data_matrix <- sapply(data_matrix,as.numeric)

#compute correlation matrix
r_matrix = cor(data_matrix, use = 'pairwise')

# run factor analysis
results<-psych::fa(r = r_matrix, nfactors = 3)
print(results)
