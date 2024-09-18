#checking 
vjp_hist <- function(x, label){
  hist(x, main = NULL)
  mu<-round(mean(x, na.rm = TRUE),3)
  sigma<-round(sd(x, na.rm = TRUE),3)
  mini<- round(min(x, na.rm = TRUE),3)
  maxi<- round(max(x, na.rm = TRUE),3)
  mtext(paste0('m = ', mu,', sd = ',sigma,', min = ',mini,', max = ',maxi))
  title(label)
}

run.cocor.dep.groups.overlap<- function(j,k,h){
  library(cocor)
  r.jk = cor.test(j,k)$estimate
  r.jh = cor.test(j,h)$estimate
  r.hk = cor.test(h,k)$estimate
  #check n is equal across vectors
  if (length(unique(sapply(list(j,k,h),length))) == 1){
    n = length(j)
  }
  
  cocor.dep.groups.overlap(r.jk,r.jh,r.hk, n)
  
}
unifactorScores <- function(X){
  library(fungible)
  complete_cases <- complete.cases(X)
  na_idx <- which(!complete.cases(X))
  X_complete = X[complete_cases]
  fac_res<-faMain(X = X_complete, numFactors = 1)
  scores<-faScores(X =X_complete, faMainObject = fac_res, Method = 'Thurstone')
  # FSI can be computed from fac res as S'RS where S is the structure (i.e. pattern or loadings %*% phi)
  FSI <- fac_res$facIndeterminacy
  #put NAs back in
  if (sum(complete_cases) != nrow(X)){
    final_scores <- insertValuesAtIdx(scores$fscores,rep(NA,sum(!complete_cases)),na_idx)
  } else {
    final_scores <- list(scores$fscores, FSI)
    names(final_scores) <- c('scores','FSI')
  }
  return(final_scores)
}

insertValuesAtIdx <- function(vec,values,idx)
{
  res<-vector(mode=mode(vec),length = length(vec)+length(idx))
  res[-idx]<-vec
  res[idx]<-values
  return(res)
}

#this takes multiple X variables and a single y vector, converts the raw X variables
#to a single vector of factor scores and then correlates fscores with y.
fscore_cor <- function(X, y){ 
  complete_cases <- complete.cases(X)
  X_complete = X[complete_cases]
  fac_res<-faMain(X = X_complete, numFactors = 1)
  scores<-faScores(X =X_complete, faMainObject = fac_res, Method = 'Thurstone')
  result<-cor.test(scores$fscores,y[complete_cases])
  return(result)
}
scatterplot<- function(x,y){
  xlab = deparse(substitute(x))
  ylab = deparse(substitute(y))
  r2 = cor.test(x, y)$estimate^2
  plot(x, y, xlab = xlab, ylab = ylab)
  mtext(round(r2,2))
}

plot_groups <- function(df,group_var){
  summary_stats <- df %>% 
    group_by(!!group_var,condition) %>%
    summarise(means = mean(responses,na.rm = TRUE),
              sd = sd(responses,na.rm = TRUE),
              n = n(),
              se = sd / sqrt(n),
              conf_int = se*2 )
  ggplot(data = summary_stats, aes(x = condition, y = means, group = !!group_var, color = !!group_var)) +
    geom_line(position = position_dodge(0.1)) + 
    geom_pointrange(aes(ymin = means-conf_int, ymax = means+conf_int), position = position_dodge(0.1))+
    theme_bw()
}

mediating_func <- function(x,m,y,control.value, treat.value, data){
  data<-subset(data, select = c(x, m, y))
  data<-data[complete.cases(data)]
  mediator_model = lm(data[[m]] ~ data[[x]])
  summary(mediator_model)
  outcome_model = lm(data[[y]] ~ data[[m]] + data[[x]])
  summary(outcome_model)
  if (missing(control.value)) {
    result<-mediate(mediator_model, outcome_model, 
                    treat = 'data[[x]]',
                    mediator = 'data[[m]]',
                    sims = 5000)
  } else {
    result<-mediate(mediator_model, outcome_model, 
                  treat = 'data[[x]]',
                  mediator = 'data[[m]]',
                  control.value = control.value,
                  treat.value = treat.value,
                  sims = 5000)
  }
  return(result) 
}


pub_ready_stats<-function(x) {
  #this is for afex aov_ez
  #correction should either be none, hf or gg
  output = NULL
  if (class(x)[1]=="afex_aov") {
    anov_table<-x[["anova_table"]]
    pest <- afex::nice(x, 'pes')
    rownames(pest)<- pest$Effect
    for (j in rownames(anov_table)) {
      fstat = round(anov_table[j,"F"], digits = 2)
      pval = round(anov_table[j,"Pr(>F)"], digits = 2)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      df1 = round(anov_table[j,"num Df"],digits =2)
      df2 = round(anov_table[j,"den Df"],digits =2)
      p_eta <- pest[j,'pes']
      pub_ready = paste0('F(',df1,',',df2,')=',fstat,', p',pval,', ', greekLetters::greeks("eta^2"), '=', p_eta)
      pub_ready = unname(cbind(j,pub_ready))
      output = unname(rbind(output,pub_ready))
    }
  }
  if (all(class(x)==c("psych","fa"))){
    #x$STATISTIC is what we want to report unless something funky's goin on
    chi<-round(x$STATISTIC,digits = 2)
    ifelse(round(x$PVAL,digits = 2)==0, pval<-"<.001",pval<-paste0('=',round(x$PVAL,digits = 2)))
    chi_rep<- paste0(greekLetters::greeks("chi^2"),'(',x$dof,')=',chi,', p',pval)
    x$CFI<-((x$null.chisq-x$null.dof)-
              (x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
    other_fits<-paste0(chi_rep,', TLI=',round(x$TLI,digits=2),
                       ', CFI=', round(x$CFI,digits=2),
                       ', RMSEA=', round(x$RMSEA,digits=2)[1])
    output = other_fits
  }
  #browser()
  if ("method" %in% names(x)){
    if (x$method == "Pearson's product-moment correlation") {
      library(MBESS)
      r = round(x$estimate,digits = 2)
      df = x$parameter
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      ci_res <- ci.cc(r, df+2)
      cis <- paste0(round(ci_res$Lower.Limit,2),',',round(ci_res$Upper.Limit,2))
      output<-paste0("r(",df,")=",r,", p",pval,', 95% CI [',cis,']')}
    if (grepl('t-test',x$method )) {
      t = unname(round(x$statistic,digits = 2))
      df = unname(round(x$parameter,digits = 2))
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      if (x$method == " Two Sample t-test"){
        d = round(t/ sqrt((df+2)),2)
        output<-paste0("t(",df,")=",t,", p",pval,", Cohen's d=",d)
      } else {
        output<-paste0("t(",df,")=",t,", p",pval)
      }
      
      
    }
    if (x$method == "Pearson's Chi-squared test"){
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      output<-paste0('X^2(',round(x$parameter,2),')=',round(x$statistic,2),
                     ', p', pval)
      
    }
    if (x$method == "Kruskal-Wallis rank sum test"){
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      output<-paste0('X^2(',round(x$parameter,2),')=',round(x$statistic,2),
                     ', p', pval)
      
    }
  }
  if (all(grepl("mediate",class(x)))) {
    library(greekLetters)
    digits = 2
    acme_b<-round(x$d0,digits = digits)
    while (acme_b == 0){
      digits = digits+1
      acme_b = round(x$d0,digits = digits)
    }
    pval<-round(x$d0.p,digits = 3)
    digits = 2
    acme_ci<-round(x$d0.ci,digits = 2)
    while (acme_ci[1] == 0 & acme_ci[2] == 0){
      digits = digits+1
      acme_ci<-round(x$d0.ci,digits = digits)
    }
    ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
    #browser()
    acme<-paste0("ACME ", greeks("beta"),"=",acme_b,", 95% CI [",acme_ci[1],", ",
                   acme_ci[2],"], p",pval)
    acme_output<-cbind('mediation',acme)
    model_m<-x$model.m
    x_to_m <- pub_ready_stats(model_m)
    x_to_m_output<- cbind('x to m',x_to_m[2,2])
    model_y<-x$model.y['model']
    model_data<-model_y$model
    m_to_y <- pub_ready_stats(lm(`data[[y]]` ~ `data[[m]]`, model_data))
    m_to_y_output <- cbind('m to y',m_to_y[2,2])
    output<-rbind(acme_output,x_to_m_output,m_to_y_output)
  }
  if (length(class(x))==2){
    if (all(class(x) %in% c('aov','lm'))){
      res_sum<-summary(x)[[1]]
      pval = round(res_sum$`Pr(>F)`[1],3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      output <- paste0('F(',res_sum$Df[1],',',res_sum$Df[2],')=',
                       round(res_sum$`F value`[1],2), 
                       ', p',pval)
    }
  }
  if (all(class(x) %in% c('anova','data.frame'))){
    res_sum<-x
    pval = round(res_sum$`Pr(>F)`[1],3)
    ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
    output <- paste0('F(',res_sum$Df[1],',',res_sum$Df[2],')=',
                     round(res_sum$`F value`[1],2), 
                     ', p',pval)
  }
  if (length(class(x))==1){
    if (class(x)=='lm'){
      #browser
      res_sum<-summary(x)$coefficients
      df<-summary(x)$df[2]
      if (nrow(res_sum)==2){
        r2_coefs = cor(unlist(x$model[1]),as.numeric(unlist(x$model[2])))^2
        r2_coefs = c(NA,r2_coefs)
        r2_text <- ', r^2 = '
      }
      if (nrow(res_sum)>2){
        library(sensemakr)
        r2_coefs<-partial_r2(x)
        r2_text <- ', partial r^2 = '
      }
      idx = 1
      for (j in rownames(res_sum)) {
        tstat = round(res_sum[j,"t value"], digits = 2)
        pval = round(res_sum[j,"Pr(>|t|)"], digits = 2)
        ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
        b = round(res_sum[j,"Estimate"],3)
        r2_coef = round(r2_coefs[idx],3)
        pub_ready = paste0('b=',b,', t(',df,')=',tstat,', p',pval,
                           r2_text,r2_coef)
        pub_ready = unname(cbind(j,pub_ready))
        output = unname(rbind(output,pub_ready))
        idx = idx+1
      }
    }
  }
  
  return(output)
}
