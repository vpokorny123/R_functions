praveen_subs<-data.table::fread("~/Downloads/capr_id_praveen.csv")
dataRequest('demo')

praveen_demo <- base::merge(praveen_subs,demo_clean, by = 'src_subject_id', all.x = TRUE)
write.csv(praveen_demo,"~/Downloads/praveen_demo.csv" )
