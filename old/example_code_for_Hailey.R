library(readxl)

df<-read_excel('/Users/victorpokorny/Downloads/EEG_DATABASE_ADAPT.xlsx')
clinical_df<- df[!df$ID_unblinded %in% c('3043','3057','3060','4003','4033'),]
clinical_df$SubID <- clinical_df$ID_blinded
clinical_df$GF_S_current_baseline <- as.numeric(clinical_df$GF_S_current_baseline)

eeg_df <- read.csv("/Users/victorpokorny/Downloads/N400_Stim_MeanAmps.csv", header = TRUE)

merged_df<-base::merge(clinical_df, eeg_df, by = 'SubID')

cor.test(merged_df$Unrelated_m_Related, merged_df$GF_S_current_baseline)
plot(merged_df$Unrelated_m_Related, merged_df$GF_S_current_baseline)

t.test(Unrelated_m_Related ~ sex, data = merged_df)
cor.test(merged_df$Unrelated_m_Related, merged_df$age)
