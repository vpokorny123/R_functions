# Install and load pscl for more pseudo-R² measures
if (!require(pscl)) install.packages("pscl")
library(pscl)
library(pROC)

source("~/Desktop/R_functions/funcs.R")
local_dir <- '/Users/victorpokorny/Library/CloudStorage/GoogleDrive-vpokorny123@gmail.com/My Drive/'
load_dir <- paste0(local_dir,'CAPR Screener Project/RData/')
load(paste0(load_dir ,'merged_df.RData'))

#for now we only care about baseline
merged_df<-merged_df[merged_df$visit == 'bl',]

#let's drop folks that are NA for group
merged_df<-merged_df[!is.na(merged_df$chr_or_not),]

#let's drop folks that don't have pqb total
merged_df<-merged_df[!is.na(merged_df$pqb_total),]

#let's drop folks that don't have prime total
merged_df<-merged_df[!is.na(merged_df$prime_total),]



#main goal is to assess whether capr tasks provide extra accuracy in determining
#who will meet criteria for CHR

#first let's look at how PQB and PRIME do in predicting CHR status
#for now logistic regression is fine because we onlt want to predict CHR vs non-CHR
model1 <- glm(chr_or_not ~ scale(pqb_total), data = merged_df, family = binomial)
pub_ready_stats(model1)
get_auc(model1)

# Optional: plot ROC curve

model2 <- glm(chr_or_not ~ scale(pqb_total) + scale(pqb_distress), data = merged_df, family = binomial)
pub_ready_stats(model2)
get_auc(model2)

anova(model1, model2)
#hmm so distress doesn't seem to matter much. Let's just focus on total for now

#now let's do prime
model <- glm(chr_or_not ~ scale(prime_total), data = merged_df, family = binomial)
pub_ready_stats(model)
get_auc(model)


#how about together?
model <- glm(chr_or_not ~ scale(prime_total) + scale(pqb_total), data = merged_df, family = binomial)
pub_ready_stats(model)
get_auc(model)

#interesting so they both predict unique variance and improve AUC a bit

#now we need to look at task stuff
library(caret)





