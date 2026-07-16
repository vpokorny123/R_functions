library(data.table)
library(dplyr)
baseline_only<-fread('/Users/victorpokorny/Downloads/baseline_only_ID.csv')
final_sample_only<-fread('/Users/victorpokorny/Downloads/final_sample_ID.csv')
big_list <- nda_study_intake_clean %>% group_by(visit,src_subject_id) %>%
  summarize(n = n())

#subset big list to only folks in the final sample
big_list_subset<-big_list[ big_list$src_subject_id %in% as.character(final_sample_only$chr_id), ]

big_list[!(big_list$src_subject_id %in% as.character(final_sample_only$chr_id)),]
table(big_list_subset$visit)

437-344

#why is there only 437? not 438? someone must have gotten dropped?
final_sample_only[!(as.character(final_sample_only$chr_id) %in% big_list$src_subject_id)]

final_sample_only[!(as.character(final_sample_only$chr_id) %in% gfr_clean$src_subject_id)]
final_sample_only[!(as.character(final_sample_only$chr_id) %in% sips_p_clean$src_subject_id)]

nda_study_intake_clean %>% group_by(visit,nda_phenotype) %>%
  summarize(n = n())
