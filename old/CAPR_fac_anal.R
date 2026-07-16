local_dir <- '/Users/victorpokorny/Library/CloudStorage/GoogleDrive-vpokorny123@gmail.com/My Drive/'
load(paste0(local_dir,'CAPR Screener Project/RData/merged_df.RData'))

#try looking at chr only at some point


library(fungible)
library(psych)
task_predictors = c(
  'kamin_blocking_score',
  'kamin_blocking_score_raw',
  'prl_total_win_switches',
  'total_yes_responses_all',
  'ch_nu',
  'glat_GAINvLA', 
  'eefrt_Highest_Mag_HardProp', 
  'LrgK', 
  'iaps_PosPleasant', 
  'meantaps_dominant',
  'context_sensitivity_all_trials',
  'inverted_faces_reported'
  #'bacs_total',
  #'sips_n_total'
  #'sips_p_total'
  #'sips_d_total'
  #'sips_g_total'
  )

task_df <- merged_df[,task_predictors]
completed_df <- task_df[complete.cases(task_df),]

fa_res<-faMain(X = completed_df, numFactors = 3, facMethod = "pca", rotate = 'varimax')
summary(fa_res)
res<-faSort(fa_res$loadings)
task_predictors[res$sortOrder]
fa_res$loadings[res$sortOrder,]

#I think the idea is to correlate with disorganized, 




#----- I also wanna try looking at just SIPS items
library(fungible)
task_predictors = c(
  paste0('n',seq(6),'_sev'),
  paste0('sips_p',seq(5),'_sev'),
  paste0('d',seq(3),'_sev')
  )


task_df <- merged_df[,task_predictors]
completed_df <- task_df[complete.cases(task_df),]

fa_res<-faMain(X = completed_df, numFactors = 3, facMethod = "fals")
summary(fa_res)

#----- I also wanna try looking at SIPS items and tasks
library(fungible)
task_predictors = c(
  paste0('n',seq(6),'_sev'),
  paste0('sips_p',seq(5),'_sev'),
  paste0('d',seq(3),'_sev'),
  #'kamin_blocking_score_raw',
  'prl_total_win_switches',
  'total_yes_responses_all',
  'ch_nu',
  'glat_GAINvLA', 
  'eefrt_Highest_Mag_HardProp', 
  'LrgK', 
  'iaps_PosPleasant', 
  'meantaps_dominant',
  'context_sensitivity_all_trials',
  'inverted_faces_reported'
)


task_df <- merged_df[,task_predictors]
completed_df <- task_df[complete.cases(task_df),]

fa_res<-faMain(X = completed_df, numFactors = 3, facMethod = "fals")
summary(fa_res)

#what if we do bassackward
bassAckward(r = cor(completed_df), nfactors = 3, fm = "minres", rotate = "oblimin", cut = .2)
