set.seed(1)
library(caret)
local_dir <- '/Users/victorpokorny/Library/CloudStorage/GoogleDrive-vpokorny123@gmail.com/My Drive/'
load(paste0(local_dir,'CAPR Screener Project/RData/merged_df.RData'))

dependent_variable = 'chr_or_not_fac'
predictors = c('prime_total','pqb_total')

split_num = .8
inTrain <- createDataPartition(
  y = merged_df[,dependent_variable],
  ## the outcome data are needed
  p = split_num,
  ## The percentage of data in the
  ## training set
  list = FALSE
)

training <- merged_df[ inTrain,]
testing  <- merged_df[-inTrain,]

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 20,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary,
  savePredictions = "final")  # enables AUC/ROC as a metric

gbmGrid <-  expand.grid(interaction.depth = c(1,2), 
                        n.trees = c(seq(5,50,2)),
                        shrinkage = c(.05,.1,.15),
                        n.minobsinnode = c(20))

lm_formula<-build_lm_formula(c(dependent_variable,predictors))

gbmFit1 <- train(form = lm_formula, 
                 data = training, 
                 method = "gbm", 
                 trControl = fitControl,
                 #na.action = 'na.omit',
                 preProcess = c('center','scale',"knnImpute"),
                 metric = 'ROC',
                 tuneGrid = gbmGrid,
                 verbose = FALSE)
plot(gbmFit1)






