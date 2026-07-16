ddk<-fread('/Users/victorpokorny/Downloads/preprocessed-analysis-input-ddk.csv')
afex::aov_ez(id = 'participant',
             within = 'subtask',
             dv = 'avg_ddk_rate',
             between = 'part_grp',
             data = ddk)

rainbow<-fread('/Users/victorpokorny/Downloads/preprocessed-analysis-input-rp.csv')
afex::aov_ez(id = 'participant',
             dv = 'avg_local_speech_rate',
             between = 'part_grp',
             data = rainbow)

pbj<-fread('/Users/victorpokorny/Downloads/preprocessed-analysis-input-pbj.csv')
afex::aov_ez(id = 'participant',
             dv = 'avg_local_speech_rate',
             between = 'part_grp',
             data = pbj)
