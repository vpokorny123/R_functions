ch_subject<- ch[ch$src_subject_id=='10005',]

#we want to look at whether strength of previous trial predicts hallucination
#let's grab idxs of all zero trials
new_df <- NULL
new_df$zero_trial<-which(ch_subject$condition==0)
#now idxs of trial before
new_df$zero_trial_minus_1<-which(ch_subject$condition==0)-1

new_df$responses<- as.numeric(ch_subject$response[new_df$zero_trial])
new_df$preceding_stim <-as.numeric(ch_subject$condition[new_df$zero_trial_minus_1])/100
res<-glm(responses ~ preceding_stim, data = new_df, family = 'binomial')
summary(res)


library(dplyr)
library(purrr)
library(broom)

# get vector of subject IDs
subjects <- unique(ch$src_subject_id)

run_model_for_subj <- function(subj_id){
  
  dat <- ch[ch$src_subject_id == subj_id, ]
  
  # Identify zero trials + preceding trials
  zero_idx  <- which(dat$condition == 0)
  prev_idx  <- zero_idx - 1
  
  # Exclude invalid (index 0 or < 1)
  keep <- prev_idx >= 1
  zero_idx  <- zero_idx[keep]
  prev_idx  <- prev_idx[keep]
  
  # Build data frame
  df <- data.frame(
    responses      = as.numeric(dat$response[zero_idx]),
    preceding_stim = as.numeric(dat$condition[prev_idx]) / 100
  )
  
  # Need variation in response for logistic regression
  if(length(unique(df$responses)) < 2){
    return(data.frame(
      src_subject_id = subj_id,
      beta_intercept = NA, se_intercept = NA,
      beta_preceding_stim = NA, se_preceding_stim = NA,
      n_trials_used = nrow(df)
    ))
  }
  
  # Fit model
  fit <- glm(responses ~ preceding_stim, data = df, family = binomial())
  
  # Extract coefficients + SEs
  out <- broom::tidy(fit)
  
  # Reshape to wide format
  tibble(
    src_subject_id      = subj_id,
    beta_intercept      = out$estimate[out$term == "(Intercept)"],
    se_intercept        = out$std.error[out$term == "(Intercept)"],
    beta_preceding_stim = out$estimate[out$term == "preceding_stim"],
    se_preceding_stim   = out$std.error[out$term == "preceding_stim"],
    n_trials_used       = nrow(df)
  )
}

# Apply to all subjects
results <- map_df(subjects, run_model_for_subj)

results$beta_preceding_stim[abs(results$beta_preceding_stim)>50] <- NA

t.test(results$beta_preceding_stim)

ch_merge<-merge(results, wide_df, by = 'src_subject_id')

anova(lm(beta_preceding_stim ~phenotype_four_groups_bl, data = ch_merge))

