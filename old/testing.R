library(wizaRdry)
clean('demo')
clean('sips_p')
clean('dd')
redcap.index()
redcap('bacs',exclude_pii = FALSE)

remove.packages("wizaRdry")
rstudioapi::restartSession()
devtools::install_github("belieflab/wizaRdry")

dict<-redcap.dict('nda_study_intake')
