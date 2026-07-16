df<-data.table::fread('/Users/victorpokorny/Downloads/threshold_data.csv')

#we can do tests between any two conditions like this
t.test(df$grating,df$grating_para, paired = TRUE)
t.test(df$grating_para,df$grating_orth, paired = TRUE)
t.test(df$noise,df$noise_noise, paired = TRUE)
t.test(df$stratified,df$stratified_para, paired = TRUE)
t.test(df$stratified_para,df$stratified_orth, paired = TRUE)
t.test(df$fibrous_para,df$fibrous_orth, paired = TRUE)

library(afex)
library(dplyr)
df_long<-tidyr::pivot_longer(df,
                             cols = c('grating','grating_para','grating_orth',
                                      'noise','noise_noise',
                                      'stratified','stratified_para','stratified_orth',
                                      'stratified_noise'),
                             names_to = 'condition',
                             values_to = 'thresholds')
df_long$texture = ifelse(grepl('grating',df_long$condition),'grating',
                         ifelse(grepl('noise',df_long$condition),'noise','stratified'))

df_long$surround = ifelse(grepl('_para',df_long$condition),'parallel',
                         ifelse(grepl('_orth',df_long$condition),'orthogonal','no_surround'))

res<-aov_ez(dv = 'thresholds',
            within = c('texture','surround'),
            id = 'V1' ,
            data = df_long)



df_long<-tidyr::pivot_longer(df,
                             cols = c('braided','braided_para',
                                      'fibrous','fibrous_para'),
                             names_to = 'condition',
                             values_to = 'thresholds')

df_long$surround = ifelse(grepl('_para',df_long$condition),'parallel',
                          ifelse(grepl('_orth',df_long$condition),'orthogonal','no_surround'))

res<-aov_ez(dv = 'thresholds',
            within = c('surround'),
            id = 'V1' ,
            data = df_long)
res

df_long<-tidyr::pivot_longer(df,
                             cols = c('braided_para','braided_orth',
                                      'fibrous_para','fibrous_orth'),
                             names_to = 'condition',
                             values_to = 'thresholds')

df_long$surround = ifelse(grepl('_para',df_long$condition),'parallel',
                          ifelse(grepl('_orth',df_long$condition),'orthogonal','no_surround'))

res<-aov_ez(dv = 'thresholds',
            within = c('surround'),
            id = 'V1' ,
            data = df_long)
res

library(lmerTest)
df_long<-tidyr::pivot_longer(df,
                             cols = -c('V1'),
                             names_to = 'condition',
                             values_to = 'thresholds')

df_long$texture = ifelse(grepl('grating',df_long$condition),'grating',
                         ifelse(grepl('noise',df_long$condition),'noise',
                                ifelse(grepl('fibrous',df_long$condition),'fibrous',
                                         ifelse(grepl('stratified',df_long$condition),'stratified','braided'))))

df_long$surround = ifelse(grepl('para',df_long$condition),'para',
                         ifelse(grepl('orth',df_long$condition),'orth','no_surround'))
df_long <- df_long[!(grepl('noise',df_long$condition) | grepl('grating',df_long$condition)),]


model1<-lmer(thresholds ~ surround + (1|V1) , data =df_long)
summary(model1)
model2<-lmer(thresholds ~ (1|V1) , data =df_long)
summary(model2)
anova(model1,model2)



