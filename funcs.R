library(dplyr)
library(rstatix)
library(tidyr)
library(ggplot2)

library(lme4)
library(ggplot2)
library(dplyr)

library(lmerTest)

library(lme4)
library(lmerTest)
library(ggplot2)

plot_partial_regression_lmer <- function(model, save_dir = NULL) {
  
  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  if (!inherits(model, "lmerModLmerTest")) {
    model <- lmerTest::as_lmerModLmerTest(model)
  }
  
  model_coefs <- summary(model)$coefficients
  data        <- model@frame
  response    <- names(data)[1]
  re_vars     <- unique(unlist(lapply(findbars(formula(model)), all.vars)))
  all_vars    <- setdiff(names(data), c(response, re_vars))
  
  re_fitted <- fitted(model) - predict(model, re.form = NA)
  y_adj     <- data[[response]] - re_fitted
  
  plot_list <- list()
  
  for (pred in all_vars) {
    
    others <- all_vars[all_vars != pred]
    
    if (length(others) > 0) {
      resid_y <- residuals(lm(as.formula(paste("y_adj ~", paste(others, collapse = " + "))),
                              data = cbind(data, y_adj = y_adj)))
    } else {
      resid_y <- y_adj - mean(y_adj)
    }
    
    pred_vals <- data[[pred]]
    if (is.numeric(pred_vals)) {
      if (length(others) > 0) {
        resid_x <- residuals(lm(as.formula(paste(pred, "~", paste(others, collapse = " + "))),
                                data = data))
      } else {
        resid_x <- pred_vals - mean(pred_vals)
      }
    } else {
      pred_numeric <- as.numeric(as.factor(pred_vals)) - 1
      tmp_df <- cbind(data, .pred_tmp = pred_numeric)
      if (length(others) > 0) {
        resid_x <- residuals(lm(.pred_tmp ~ ., data = tmp_df[, c(".pred_tmp", others)]))
      } else {
        resid_x <- pred_numeric - mean(pred_numeric)
      }
    }
    
    coef_rows <- model_coefs[grepl(pred, rownames(model_coefs)), , drop = FALSE]
    
    if (nrow(coef_rows) >= 1) {
      slope   <- coef_rows[1, "Estimate"]
      p_val   <- min(coef_rows[, "Pr(>|t|)"])
      p_label <- if (p_val < 0.001) "p < .001" else sprintf("p = %.3f", p_val)
      if (nrow(coef_rows) > 1) p_label <- paste0(p_label, "*")
    } else {
      slope   <- NA
      p_label <- "p = N/A"
    }
    
    label   <- sprintf("β = %.4f  |  %s", slope, p_label)
    plot_df <- data.frame(resid_x = resid_x, resid_y = resid_y)
    
    p <- ggplot(plot_df, aes(x = resid_x, y = resid_y)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
      geom_point(alpha = 0.5, color = "black", size = 2) +
      geom_smooth(method = "lm", se = TRUE,
                  color = "black", fill = "black", alpha = 0.15, linewidth = 1) +
      annotate("text", x = Inf, y = Inf, label = label,
               hjust = 1.05, vjust = 1.8, size = 5, color = "grey30") +
      labs(
        title = paste("Partial Regression:", pred),
        x     = paste0("e(", pred, " | others)"),
        y     = paste0("e(", response, " | others)")
      ) +
      theme_bw(base_size = 18) +
      theme(
        plot.title  = element_text(face = "bold", size = 12),
        panel.grid  = element_line(color = "grey92"),
        axis.title  = element_text(size = 10),
        plot.margin = margin(10, 15, 10, 10)
      )
    
    plot_list[[pred]] <- p
    print(p)
    
    if (!is.null(save_dir)) {
      filename <- file.path(save_dir, paste0("partial_regression_", pred, ".svg"))
      ggsave(filename, plot = p, device = "svg", width = 6, height = 5)
      message("Saved: ", filename)
    }
  }
  
  invisible(plot_list)
}
pub_ready_reg_table <- function(x, save_file = NULL) {
  
  is_lmer <- inherits(x, c("lmerMod", "lmerModLmerTest"))
  
  res_sum      <- summary(x)
  coefficients <- round(res_sum$coefficients, 3)
  
  # --- P-values ---
  if ("Pr(>|t|)" %in% colnames(coefficients)) {
    ps <- unname(coefficients[, 'Pr(>|t|)'])
    ps <- ifelse(ps == 0, "<.001", as.character(ps))
    coefficients[, 'Pr(>|t|)'] <- ps
  } else if (is_lmer) {
    coefficients <- cbind(coefficients, `Pr(>|t|)` = NA_character_)
    message("Note: p-values unavailable for lmer models without lmerTest.")
  }
  
  # --- Estimates ---
  betas <- unname(coefficients[, 'Estimate'])
  betas <- ifelse(betas == 0, "<.001", as.character(betas))
  coefficients[, 'Estimate'] <- betas
  
  # --- Clean rownames ---
  rownames(coefficients) <- gsub('`', '', rownames(coefficients))
  
  # --- Partial R² / eta² ---
  if (is_lmer) {
    
    library(effectsize)
    coef_mat <- res_sum$coefficients
    t_vals   <- coef_mat[, 't value']
    
    if ("df" %in% colnames(coef_mat)) {
      dfs <- coef_mat[, 'df']
    } else {
      dfs <- rep(df.residual(x), length(t_vals))
    }
    
    eta2     <- t_to_eta2(t = t_vals, df_error = dfs)$Eta2_partial
    r2_coefs <- round(eta2, 3)
    r2_coefs <- ifelse(r2_coefs == 0, "<.001", as.character(r2_coefs))
    
  } else {
    
    library(sensemakr)
    r2_coefs <- round(partial_r2(x), 3)
    r2_coefs <- ifelse(r2_coefs == 0, "<.001", as.character(r2_coefs))
  }
  
  reg_table <- cbind(coefficients, r2_coefs)
  n_cols    <- ncol(reg_table)
  
  # --- Model-level fit statistics ---
  if (is_lmer) {
    
    library(MuMIn)
    r2_vals           <- r.squaredGLMM(x)
    rsquared_marginal <- round(r2_vals[1, "R2m"], 3)
    rsquared_cond     <- round(r2_vals[1, "R2c"], 3)
    
    fit_row <- c(
      paste0('Marginal R^2 (fixed) = ', rsquared_marginal, '; ',
             'Conditional R^2 (fixed + random) = ', rsquared_cond),
      rep("", n_cols - 1)
    )
    
  } else {
    
    rsquared      <- round(res_sum$r.squared, 3)
    rsquared_adj  <- round(res_sum$adj.r.squared, 3)
    fstat         <- res_sum$fstatistic['value']
    df1           <- res_sum$fstatistic['numdf']
    df2           <- res_sum$fstatistic['dendf']
    rsquared_pval <- round(pf(fstat, df1, df2, lower.tail = FALSE), 3)
    rsquared_pval <- ifelse(rsquared_pval == 0, "<.001", paste0('=', rsquared_pval))
    
    fit_row <- c(
      paste0('Full Model R^2 = ', rsquared, '; ',
             'Adjusted R^2 = ', rsquared_adj, '; ',
             'F(', df1, ',', df2, ')=', round(fstat, 2),
             ', p', rsquared_pval),
      rep("", n_cols - 1)
    )
  }
  
  pub_ready_reg_table <- rbind(reg_table, ` ` = fit_row)
  
  # --- Random effects rows (lmer only) ---
  if (is_lmer) {
    
    vc    <- lme4::VarCorr(x)
    vc_df <- as.data.frame(vc)
    
    blank_row <- matrix("", nrow = 1, ncol = n_cols,
                        dimnames = list("---Random Effects---", NULL))
    pub_ready_reg_table <- rbind(pub_ready_reg_table, blank_row)
    
    for (i in seq_len(nrow(vc_df))) {
      grp      <- vc_df$grp[i]
      var_name <- ifelse(is.na(vc_df$var2[i]),
                         vc_df$var1[i],
                         paste0(vc_df$var1[i], " | ", vc_df$var2[i]))
      variance <- round(vc_df$vcov[i], 3)
      sd_val   <- round(vc_df$sdcor[i], 3)
      
      re_row        <- matrix("", nrow = 1, ncol = n_cols,
                              dimnames = list(paste0(grp, ": ", var_name), NULL))
      re_row[1, 1]  <- paste0('Variance = ', variance)
      pub_ready_reg_table <- rbind(pub_ready_reg_table, re_row)
    }
  }
  
  # --- Save ---
  if (!is.null(save_file)) {
    write.csv(pub_ready_reg_table, file = save_file)
    message("Table saved to: ", save_file)
  }
  
  return(pub_ready_reg_table)
}

vjp_rectangle_corrplot <- function(row_variables, column_variables, adjust = 'fdr'){
  #NOTE if you're having trouble running this. It may be because you have
  #a single variables in row_variables. this function CAN handle this, but
  #the input needs to be df[rows_variables], NOT df[,row_variables]. The
  #latter way removes the name from the vector which causes problems later
  #on
  #browser()
  library(corrplot)
  cor_res <- psych::corr.test(x = row_variables, y =column_variables, 
                              use = 'pairwise', adjust = adjust)
  df = cbind(row_variables,column_variables)
  #partial.r(df, x = colnames(row_variables), y =colnames(column_variables), 
  #          use = 'pairwise')
  cor_mat <- cor_res$r
  p_mat <- cor_res$p.adj
  p_mat_uncorrected <- cor_res$p
  #corrplot gets mad if none of the p-values are significant AND we are doing   
  #a rectangle corr plot. so let's avoid this with an if statement
  if (any(p_mat<.05)){
    corrplot(cor_mat, p.mat = p_mat,
             method = c('circle'),
             tl.cex = .8, number.cex = .65,
             sig.level = .05, diag = TRUE,
             insig = 'label_sig', pch.cex = 1.3,
             pch.col = 'black', tl.col = 'black',
             col = rev(COL2('RdBu', 100))) 
  } else {
    corrplot(cor_mat,
             method = c('circle'),
             tl.cex = .6, number.cex = .65,
             sig.level = .05, diag = TRUE,
             insig = 'label_sig', pch.cex = 1,
             pch.col = 'black',tl.col = 'black',
             col = rev(COL2('RdBu', 100))) 
  }
  output <- list(round(p_mat,3),round(p_mat_uncorrected,3),cor_mat)
  names(output)<-c('p vals adjusted','p vals not adjusted','cors')
  return(output)
  
}

vjp_corrplot <- function(df,p){
  library(corrplot)
  library(psych)
  #cor_mat <- cor(df, use = "pairwise.complete.obs")
  cor_res <- psych::corr.test(df, adjust = 'fdr')
  cor_mat <- cor_res$r
  p_mat <- cor_res$p
  corrplot(cor_mat, p.mat = p_mat,
           method = c('circle'),type = "upper",
           tl.cex = .6, number.cex = .65,
           sig.level = .05, diag = FALSE,
           insig = 'label_sig', pch.cex = 1,
           pch.col = 'black',
           col = rev(COL2('RdBu', 100)))
  
}
#checking 
build_lm_formula <-function(x, interactions = FALSE){
  for (j in seq(x)){
    variable = x[j]
    if (j == 1){
      lm_formula = paste0('`',variable,'` ~ `')
    } else if (j == 2) {
      lm_formula = paste0(lm_formula, variable, '`')
    }else if (j>2 & interactions == TRUE){
      lm_formula = paste0(lm_formula,' * `',variable, '`')
    }else if (j>2 & interactions ==FALSE){
      lm_formula = paste0(lm_formula,' + `',variable, '`')
    }
  }
  return(formula(lm_formula))
}
impute_glm <- function(y,x, data){
  model_variables <- c(y, x)
  
  #just constructing the lm_formula
  for (j in seq(model_variables)){
    variable = model_variables[j]
    if (j == 1){
      lm_formula = paste0(variable,' ~ ')
    } else if (j == 2) {
      lm_formula = paste0(lm_formula,variable)
    }else {
      lm_formula = paste0(lm_formula,' + ',variable)
    }
  }
  # for (j in seq(model_variables)){
  #   variable = model_variables[j]
  #   if (j == 1){
  #     lm_formula = paste0('scale(',variable,') ~ ')
  #   } else if (j == 2) {
  #     lm_formula = paste0(lm_formula,'scale(',variable,')')
  #   }else {
  #     lm_formula = paste0(lm_formula,' + scale(',variable,')')
  #   }
  # }
  #lm_formula<-build_lm_formula(model_variables)
  
  #now we create an impute object
  imputed_data <- data %>% 
    select(all_of(model_variables)) %>%
    mice(pred = quickpred(.), seed = 123, m = 10) 
  
  #now we can fit 5 linear regressions and pool across them 
  fit<-with(imputed_data, lm(formula = formula(lm_formula)))
  results<-summary(pool(fit))
  rsquared_adj<-round(pool.r.squared(fit, adjusted = TRUE)[1,1],3)
  rsquared<-round(pool.r.squared(fit, adjusted = FALSE)[1,1],3)
  round_ps<-round(results$p.value,3)
  #browser()
  cleaner_results<-cbind(results, round_ps)
  final_results <- list(cleaner_results, 
                        rsquared_adj, 
                        rsquared)
  names(final_results) <- c("main_results","rsquared_adj","rsquared")
  return(final_results)
  
}

impute_logistic <- function(y,x, data){
  model_variables <- c(y, x)
  lm_formula<-build_lm_formula(model_variables)
  
  #browser()
  #now we create an impute object
  imputed_data <- data %>% 
    select(all_of(model_variables)) %>%
    mice(., m = 5, method = 'pmm', seed = 123) 
  
  # Create a wrapper function that accepts data and uses the outer formula
  my_glm <- function(data) {
    glm(formula = lm_formula, family = binomial, data = data)
  }
  #now we can fit 5 linear regressions and pool across them 
  fit <- with(imputed_data, my_glm(data = data))
  results<-summary(pool(fit))
  
  # Step 2: Extract completed datasets
  completed_datasets <- lapply(1:imputed_data$m, function(i) complete(imputed_data, i))
  
  # Step 3: For each completed dataset, fit model and compute AUC
  auc_list <- sapply(completed_datasets, function(d) {
    model <- glm(formula = lm_formula, data = d, family = binomial)
    preds <- predict(model, type = "response")
    roc_obj <- roc(d[[y]], preds)
    auc(roc_obj)
  })
  
  # Step 4: Average AUC
  mean_auc <- round(mean(auc_list),3)
  print(paste0('Average AUC =',mean_auc))
  return(list(fit = fit, results = results))
  
}


vjp_hist <- function(x, label = NULL, gaussian = TRUE){
  #browser()
  mu    <-round(mean(x, na.rm = TRUE),3)
  sigma <-round(sd(x, na.rm = TRUE),3)
  mini  <- round(min(x, na.rm = TRUE),3)
  maxi  <- round(max(x, na.rm = TRUE),3)
  
  if (gaussian==TRUE) {
    hist(x, main = NULL, breaks = 10, freq = FALSE)
    curve(dnorm(x, mean = mu, sd = sigma), 
          col = "blue", lwd = 2, add = TRUE)
  } else {
    hist(x, main = NULL, breaks = 'Scott', freq = TRUE)
  }
  mtext(paste0('m = ', mu,', sd = ',sigma,', min = ',mini,', max = ',maxi))
  if (!is.null(label)){
    title(label)
  } 
  
}

run.cocor.dep.groups.overlap<- function(j,k,h){
  
  library(cocor)
  #testing whether jk differs from jh
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
  X_complete = X[complete_cases,]
  #browser()
  fac_res<-faMain(X = X_complete, numFactors = 1)
  scores<-faScores(X =X_complete, faMainObject = fac_res, Method = 'Thurstone')
  # FSI can be computed from fac res as S'RS where S is the structure (i.e. pattern or loadings %*% phi)
  FSI <- fac_res$facIndeterminacy
  #put NAs back in
  if (sum(complete_cases) != nrow(X)){
    final_scores <- insertValuesAtIdx(scores$fscores,rep(NA,sum(!complete_cases)),na_idx)
    final_scores <- list(final_scores, fac_res)
  } else {
    final_scores <- list(scores$fscores, fac_res)
  } 
  names(final_scores) <- c('scores','fa.res')
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
scatterplot<- function(x,y,xlab=NULL, ylab= NULL, tail = NULL, jitter_x = NULL,
                       jitter_y = NULL, grouping_var = NULL,
                       save_name = NULL){
  library(scales)
  if (!is.null(save_name)){
    svg(save_name)
  }
  #browser()
  if (!is.null(grouping_var)){
    group_levels <- unique(grouping_var)
    group_colors <- setNames(hue_pal()(length(group_levels)), group_levels)
    alpha <- 0.4
    group_colors <- sapply(group_colors, function(clr) adjustcolor(clr, alpha.f = alpha))
    col <- group_colors[grouping_var]
    
  } else {
    col = 'black'
  }
  
  if (is.null(xlab) & is.null(ylab)){
    xlab = deparse(substitute(x))
    ylab = deparse(substitute(y))
  }
  
  if (!is.null(tail)){
    res = cor.test(x, y, alternative = tail)
  } else {
    res = cor.test(x, y)
  }
  model<-lm(y ~ x)
  #browser()
  if (!is.null(jitter_x)){
    plot(jitter(x, factor = jitter_x), y, xlab = xlab, ylab = ylab, pch = 20, 
         col = col)
  } else {
    plot(x, y, xlab = xlab, ylab = ylab, pch = 20, col = col)
  }
  
  if (!is.null(grouping_var)){
    #par(mar = c(5, 4, 4, 10))           # enlarge right margin
    legend("topright", inset = c(-0.25, 0), xpd = TRUE,
           legend = group_levels, col = group_colors, pch = 16, title = "Group")}
  abline(model,lwd = 2)
  mtext(pub_ready_stats(res))
  print(pub_ready_stats(res))
  if (!is.null(save_name)){
    dev.off()
  }
  
  
}

ggscatterplot<- function(df, x, y, xlab=NULL, ylab= NULL, tail = NULL, jitter_x = NULL,
                         jitter_y = NULL, grouping_var = NULL, group_specific_lines = FALSE,
                         save_name = NULL, save_width = 5, save_height = 4.5){
  font_size = 15
  #browser()
  library(ggplot2)
  if (is.null(xlab)){
    xlab = x
  }
  if (is.null(ylab)){
    ylab = y
  }
  
  #browser()
  final_plot<-ggplot(df, aes(x = !!as.name(x), y = !!as.name(y), 
                             color = !!as.name(grouping_var),
                             group = !!as.name(grouping_var))) +
    geom_point(alpha = .7) 
  #browser()
  if (group_specific_lines == FALSE){
    stats<-pub_ready_stats(cor.test(df[[x]],df[[y]]))
    print(stats)
    final_plot <- final_plot + geom_smooth(method = "lm", se = TRUE, color = "black", group = 1) +
      xlab(xlab) + ylab(ylab)+
      theme_bw(base_size = font_size) + 
      annotate('text',  x=Inf, y = Inf, label = stats, vjust=1.6, hjust=1.05) #,size = 4 )
  } else {
    #browser()
    formula<-build_lm_formula(c(y,x,grouping_var), interactions = TRUE)
    stats<-pub_ready_stats(anova(lm(formula = formula, df)))
    stats<- paste(stats[3,],collapse = ' ')  #assuming we want interaction stat
    print(stats)
    #browser()
    final_plot <- final_plot + geom_smooth(method = "lm", se = TRUE, alpha = .1) +
      xlab(xlab) + ylab(ylab)+
      theme_bw(base_size = font_size) + 
      annotate('text',  x=Inf, y = Inf, label = stats, vjust=1.6, hjust=1.05) #,size = 4 )
  }
  if (!is.null(save_name)){
    ggsave(save_name,width = save_width,height = save_height,dpi = 500)
  }        
  return(final_plot)
  
}

plot_conditions <- function(df,cond_var, dv , xlab = NULL, ylab = NULL,
                            save_name = NULL, jitter_height = 0, jitter_width = 0,
                            save_width = 5, save_height = 4.5){
  #browser()
  summary_stats <- df %>% 
    group_by(!!as.name(cond_var)) %>%
    summarise(means = mean(!!as.name(dv),na.rm = TRUE),
              sd = sd(!!as.name(dv),na.rm = TRUE),
              n = n(),
              se = sd / sqrt(n),
              conf_int = se*2 )
  
  final_plot<-ggplot(data = df, 
                     aes(x = !!as.name(cond_var),y = !!as.name(dv))) +
    geom_point(position = position_jitter(
      width = jitter_width, height = jitter_height),
      alpha = .2) + 
    #geom_line(position = position_dodge(0.1)) + 
    geom_pointrange(data = summary_stats, aes(y = means, ymin = means-conf_int, 
                                              ymax = means+conf_int,), color = 'black',
                    position = position_nudge(x = .2))+
    theme_bw() +
    ylab(ylab) + xlab(xlab)  
  #annotate('text',  x=Inf, y = Inf, label = stats, vjust=1.6, hjust=1.05 
  #size = 4
  
}


plot_groups <- function(
    df =NULL, group_var=NULL, dv=NULL,
    condition = NULL, task_var = NULL, facet_var = NULL,
    xlab = NULL, ylab = NULL,
    dodge_width = 0.2,              # spacing of group summaries
    jitter_height = 0,
    save_name = NULL, save_width = 5, save_height = 4.5,
    
    # raw points
    show_points = TRUE,
    points_alpha = 0.35,
    points_size  = 1.2,
    points_jitter_width  = NULL,
    points_jitter_height = NULL,
    
    # summary appearance
    summary_shape_size = 2.5         # NEW: size of summary points
) {
  if (is.null(xlab)) xlab <- if (!is.null(condition)) condition else group_var
  if (is.null(ylab)) ylab <- dv
  font_size <- 15
  
  if (is.null(points_jitter_height)) points_jitter_height <- jitter_height
  
  # ---- basic cleaning ----
  df <- df %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[dv]])) %>%
    dplyr::mutate("{group_var}" := droplevels(as.factor(.data[[group_var]])))
  
  # ---- CONDITION branch ----
  if (!is.null(condition)) {
    
    df <- df %>%
      dplyr::filter(!is.na(.data[[condition]])) %>%
      dplyr::mutate("{condition}" := droplevels(as.factor(.data[[condition]])))
    
    facet_vars <- character(0)
    
    if (!is.null(task_var)) {
      df <- df %>%
        dplyr::filter(!is.na(.data[[task_var]])) %>%
        dplyr::mutate("{task_var}" := droplevels(as.factor(.data[[task_var]])))
      facet_vars <- c(facet_vars, task_var)
    }
    
    if (!is.null(facet_var)) {
      df <- df %>%
        dplyr::filter(!is.na(.data[[facet_var]])) %>%
        dplyr::mutate("{facet_var}" := droplevels(as.factor(.data[[facet_var]])))
      facet_vars <- c(facet_vars, facet_var)
    }
    
    if (is.null(points_jitter_width)) points_jitter_width <- 0.3 * dodge_width
    
    group_vars <- c(group_var, condition, facet_vars)
    
    summary_stats <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        mean = mean(.data[[dv]], na.rm = TRUE),
        n = dplyr::n(),
        sd = sd(.data[[dv]], na.rm = TRUE),
        se = sd / sqrt(n),
        ci = 2 * se,
        .groups = "drop"
      ) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.factor), droplevels))
    
    p <- ggplot2::ggplot(
      summary_stats,
      ggplot2::aes(
        x = .data[[condition]],
        y = mean,
        group = .data[[group_var]],
        color = .data[[group_var]],
        shape = .data[[group_var]]
      )
    )
    
    # raw points (background)
    if (show_points) {
      p <- p +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(
            x = .data[[condition]],
            y = .data[[dv]],
            color = .data[[group_var]],
            shape = .data[[group_var]]
          ),
          inherit.aes = FALSE,
          show.legend = FALSE,
          alpha = points_alpha,
          size = points_size,
          position = ggplot2::position_jitterdodge(
            jitter.width  = points_jitter_width,
            jitter.height = points_jitter_height,
            dodge.width   = dodge_width
          )
        )
    }
    
    # summaries
    p <- p +
      ggplot2::geom_line(
        position = ggplot2::position_dodge(dodge_width),
        show.legend = FALSE
      ) +
      ggplot2::geom_pointrange(
        ggplot2::aes(ymin = mean - ci, ymax = mean + ci),
        position = ggplot2::position_dodge(dodge_width),
        fatten = summary_shape_size / 2,   # controls point size
        linewidth = 0.8
      ) +
      ggplot2::scale_x_discrete(drop = TRUE) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(ylab) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      cb_scale_color(levels = levels(df[[group_var]])) +
      cb_scale_shape(levels = levels(df[[group_var]]))
    
    if (length(facet_vars) > 0) {
      p <- p + ggplot2::facet_wrap(
        stats::as.formula(paste("~", paste(facet_vars, collapse = "+"))),
        drop = TRUE, scales = "free_x"
      )
    }
    
  } else {
    
    # ---- NO-CONDITION branch ----
    if (is.null(points_jitter_width)) points_jitter_width <- 0.15
    
    summary_stats <- df %>%
      dplyr::group_by(.data[[group_var]]) %>%
      dplyr::summarise(
        mean = mean(.data[[dv]], na.rm = TRUE),
        n = dplyr::n(),
        sd = sd(.data[[dv]], na.rm = TRUE),
        se = sd / sqrt(n),
        ci = 2 * se,
        .groups = "drop"
      )
    
    stats_txt <- if (nlevels(df[[group_var]]) > 2) {
      pub_ready_stats(stats::aov(df[[dv]] ~ df[[group_var]]))
    } else {
      pub_ready_stats(stats::t.test(df[[dv]] ~ df[[group_var]], var.equal = TRUE))
    }
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[[group_var]],
        y = .data[[dv]],
        color = .data[[group_var]],
        shape = .data[[group_var]]
      )
    )
    
    if (show_points) {
      p <- p +
        ggplot2::geom_point(
          show.legend = FALSE,
          alpha = points_alpha,
          size = points_size,
          position = ggplot2::position_jitter(
            width  = points_jitter_width,
            height = points_jitter_height
          )
        )
    }
    
    p <- p +
      ggplot2::geom_pointrange(
        data = summary_stats,
        ggplot2::aes(
          x = .data[[group_var]],
          y = mean,
          ymin = mean - ci,
          ymax = mean + ci
        ),
        inherit.aes = FALSE,
        color = "black",
        fatten = summary_shape_size / 2,
        position = ggplot2::position_nudge(x = 0.2),
        show.legend = FALSE
      ) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(ylab) +
      ggplot2::annotate(
        "text", x = Inf, y = Inf, label = stats_txt,
        vjust = 1.6, hjust = 1.05
      ) +
      cb_scale_color(levels = levels(df[[group_var]])) +
      cb_scale_shape(levels = levels(df[[group_var]]))
  }
  
  if (!is.null(save_name)) {
    ggplot2::ggsave(save_name, p, width = save_width, height = save_height, dpi = 500)
  }
  
  return(p)
}





mediating_func <- function(x,m,y,control.value, treat.value, sims, boot = FALSE, data){
  library(mediation)
  data<-base::subset(data, select = c(x, m, y))
  #browser()
  data<-data[complete.cases(data),]
  mediator_model = lm(data[[m]] ~ data[[x]])
  summary(mediator_model)
  outcome_model = lm(data[[y]] ~ data[[m]] + data[[x]])
  summary(outcome_model)
  if (missing(sims)){
    sims = 5000
  }
  print(paste0('starting mediation with ', sims, ' simulations. If taking too long',
               ' set sims to something lower'))
  if (missing(control.value)) {
    result<-mediate(mediator_model, outcome_model, 
                    treat = 'data[[x]]',
                    mediator = 'data[[m]]',
                    sims = sims,
                    boot = boot)
  } else {
    result<-mediate(mediator_model, outcome_model, 
                    treat = 'data[[x]]',
                    mediator = 'data[[m]]',
                    control.value = control.value,
                    treat.value = treat.value,
                    sims = sims,
                    boot = boot)
  }
  return(result) 
}


pub_ready_stats<-function(x, one_tail = FALSE) {
  #this is for afex aov_ez
  #correction should either be none, hf or gg
  output = NULL
  if (class(x)[1]=="afex_aov") {
    anov_table<-x[["anova_table"]]
    pest <- afex::nice(x, 'pes')
    rownames(pest)<- pest$Effect
    for (j in rownames(anov_table)) {
      fstat = round(anov_table[j,"F"], digits = 2)
      pval = round(anov_table[j,"Pr(>F)"], digits = 3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      df1 = round(anov_table[j,"num Df"],digits =2)
      df2 = round(anov_table[j,"den Df"],digits =2)
      p_eta <- pest[j,'pes']
      pub_ready = paste0('F(',df1,',',df2,')=',fstat,', p',pval,', ',
                         greekLetters::greeks("eta^2"), '=', p_eta)
      pub_ready = unname(cbind(j,pub_ready))
      output = unname(rbind(output,pub_ready))
    }
  }
  if (all(class(x)==c("psych","fa"))){
    #x$STATISTIC is what we want to report unless something funky's goin on
    chi<-round(x$STATISTIC,digits = 2)
    ifelse(round(x$PVAL,digits = 3)==0, pval<-"<.001",pval<-paste0('=',round(x$PVAL,digits = 3)))
    chi_rep<- paste0(greekLetters::greeks("chi^2"),'(',x$dof,')=',chi,', p',pval)
    x$CFI<-((x$null.chisq-x$null.dof)-
              (x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
    other_fits<-paste0(chi_rep,', TLI=',round(x$TLI,digits=2),
                       ', CFI=', round(x$CFI,digits=2),
                       ', RMSEA=', round(x$RMSEA,digits=2)[1])
    output = other_fits
  }
  if ("method" %in% names(x)){
    if (x$method == "Pearson's product-moment correlation") {
      #browser()
      library(MBESS)
      r = round(x$estimate,digits = 2)
      df = x$parameter
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      cis <- paste0(round(x$conf.int[1],2),',',round(x$conf.int[2],2))
      output<-paste0("r(",df,")=",r,", p",pval,', 95% CI [',cis,']')}
    if (x$method == "Spearman's rank correlation rho") {
      #browser()
      library(MBESS)
      r = round(x$estimate,digits = 2)
      pval = round(x$p.value,3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      output<-paste0("r=",r,", p",pval)}
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
    if (grepl("Pearson's Chi-squared test",x$method)){
      #browser()
      pval = round(x$p.value,3)
      df = round(x$parameter,2)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      ifelse(is.na(df), df<-"",df<-paste0('(',df,')'))
      output<-paste0('X^2',df,'=',round(x$statistic,2),
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
    acme_b<-round(x$d.avg,digits = digits)
    ade_b = round(x$z.avg,digits = digits)
    while (acme_b == 0){
      digits = digits+1
      acme_b = round(x$d.avg,digits = digits)
    }
    while (ade_b == 0){
      digits = digits+1
      ade_b = round(x$z.avg,digits = digits)
    }
    
    acme_pval <- round(x$d.avg.p,digits = 3)
    ade_pval  <- round(x$z.avg.p,digits = 3)
    
    ifelse(acme_pval==0, 
           acme_pval<-"<.001",
           acme_pval<-paste0('=',acme_pval))
    ifelse(ade_pval==0, 
           ade_pval<-"<.001",
           ade_pval<-paste0('=',ade_pval))
    #browser()
    acme_ci<-round(x$d.avg.ci ,digits = digits)
    ade_ci<-round(x$z.avg.ci ,digits = digits)
    while (acme_ci[1] == 0 & acme_ci[2] == 0){
      digits = digits+1
      acme_ci<-round(x$d0.ci,digits = digits)
    }
    while (ade_ci[1] == 0 & ade_ci[2] == 0){
      digits = digits+1
      ade_ci<-round(x$z.avg.ci,digits = digits)
    }
    
    acme<-paste0("ACME ", greeks("beta"),"=",acme_b,", 95% CI [",acme_ci[1],", ",
                 acme_ci[2],"], p",acme_pval)
    ade <- paste0("ADE ", greeks("beta"),"=",ade_b,
                  ", 95% CI [",
                  ade_ci[1],", ",ade_ci[2],"], p",ade_pval)
    acme_output<-cbind('mediation',acme)
    ade_output <-cbind('direct effect',ade)
    model_m<-x$model.m
    x_to_m <- pub_ready_stats(model_m)
    x_to_m_output<- cbind('x to m',x_to_m[2,2])
    model_y<-x$model.y['model']
    model_data<-model_y$model
    m_to_y <- pub_ready_stats(lm(`data[[y]]` ~ `data[[m]]`, model_data))
    m_to_y_output <- cbind('m to y',m_to_y[2,2])
    
    output<-rbind(acme_output,ade_output, x_to_m_output,m_to_y_output)
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
    #browser()
    if (all(c("npar","AIC") %in% names(x))){
      #browser()
      reduced_model_aic<-round(x$AIC[1],1)
      full_model_aic<-round(x$AIC[2],1)
      chisq<-round(x$Chisq[2],2)
      df <- x$Df[2]
      pval <- round(x$`Pr(>Chisq)`[2],3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      output <- paste0('Full AIC=', full_model_aic, 
                       ', Reduced AIC=',reduced_model_aic,
                       ', X^2(',df,')=', chisq,
                       ', p',pval)
    } else if (all(c("Sum Sq","Mean Sq","NumDF","DenDF","F value","Pr(>F)") %in% names(x))){
      res_sum<-x
      idx = 1
      rownames = rownames(res_sum) 
      for (j in rownames) {
        f = round(res_sum[j,'F value'],2)
        df1 = res_sum[idx,'NumDF']
        df2 = round(res_sum[idx,'DenDF'],2)
        pval = round(res_sum[j,"Pr(>F)"], digits = 3)
        ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
        pub_ready = paste0('F(',df1,',',df2,')=',f,', p',pval)
        pub_ready = unname(cbind(j,pub_ready))
        output = unname(rbind(output,pub_ready))
        idx = idx+1
      }
      
    } else {
      res_sum<-x
      idx = 1
      rownames = rownames(res_sum) 
      rownames = rownames[-length(rownames)] #drop residual term for convenience
      for (j in rownames) {
        f = round(res_sum[j,'F value'],3)
        df1 = res_sum$Df[idx]
        df2 = res_sum$Df[nrow(res_sum)]
        pval = round(res_sum[j,"Pr(>F)"], digits = 3)
        ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
        pub_ready = paste0('F(',df1,',',df2,')=',f,', p',pval)
        pub_ready = unname(cbind(j,pub_ready))
        output = unname(rbind(output,pub_ready))
        idx = idx+1
      }
    }
  }
  if (length(class(x))==1){
    if (class(x)=='lm'){
      #browser()
      res_sum<-summary(x)$coefficients
      df<-summary(x)$df[2]
      if (nrow(res_sum)==2){
        #browser()
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
        pval = round(res_sum[j,"Pr(>|t|)"], digits = 3)
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
  #browser()
  if (class(x)[1] == "lmerModLmerTest" ){
    res_sum<-summary(x)$coefficients
    df<- res_sum
    idx = 1
    for (j in rownames(res_sum)) {
      tstat = round(res_sum[j,"t value"], digits = 2)
      pval = round(res_sum[j,"Pr(>|t|)"], digits = 3)
      df = round(res_sum[j,"df"], digits = 3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      b = round(res_sum[j,"Estimate"],3)
      #r2_coef = round(r2_coefs[idx],3)
      pub_ready = paste0('b=',b,', t(',df,')=',tstat,', p',pval)
      #r2_text,r2_coef)
      pub_ready = unname(cbind(j,pub_ready))
      output = unname(rbind(output,pub_ready))
      idx = idx+1
    }
  }
  if (class(x)[1] == "QuadTypeIndependenceTest" ){
    chi_stat = round(x@statistic@teststatistic,2)
    df = x@statistic@df 
    p_val = round(pchisq(chi_stat, df,lower.tail = FALSE),3)
    p_val = ifelse(p_val==0, '<.001',paste0('=',p_val))
    output <- paste0('X^2(',df,')=',chi_stat,', p', p_val)
  }
  #browser()
  if (class(x)[1] == "ScalarIndependenceTest" ){
    #browser()
    z_stat<-round(x@statistic@teststatistic,3)
    #assuming we want two-tailed
    p_val <- round(2 * (1 - pnorm(abs(z_stat))),3)
    p_val = ifelse(p_val==0, '<.001',paste0('=',p_val))
    output <- paste0('Z=',z_stat,', p=',p_val)
    #browser()
  }
  #browser()
  if (!isS4(x)){
    if (!is.null(x$family)){
      if (x$family["family"]=='binomial' & x$family["link"]=='logit'){
        
        res_sum <- summary(x)$coefficients
        
        #compute confidence intervals
        #ci_sum<-exp(cbind(OR = coef(x), confint(x)))
        ci_wald <- exp(confint.default(x))
        #browser()
        #loop through effects
        idx = 1
        for (j in rownames(res_sum)) {
          #log_odds coef
          log_odds_coef = res_sum[j,"Estimate"]
          log_odds_coef_rd<- round(log_odds_coef,digits = 2)
          z_stat <- round(res_sum[j,"z value"], digits = 2)
          odds_ratio_coef = round(exp(log_odds_coef),2)
          lower_ci<-round(ci_wald[j,'2.5 %'],2)
          upper_ci<-round(ci_wald[j,'97.5 %'],2)
          pval = round(res_sum[j,"Pr(>|z|)"], digits = 3)
          ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
          #r2_coef = round(r2_coefs[idx],3)
          pub_ready = paste0('b=',log_odds_coef_rd, ', p',pval,', OR=',odds_ratio_coef, 
                             ', 95% CI [',lower_ci,',',upper_ci,']')
          pub_ready = unname(cbind(j,pub_ready))
          output = unname(rbind(output,pub_ready))
          idx = idx+1 }
      }
    }
  }
  if(any(class(x)=='clm')){
    res_sum <- summary(x)$coefficients
    #can only get conf interval for predictor no threshold coefficients
    ci<-confint(x)
    #loop through effects
    idx = 1
    for (j in rownames(ci)) {
      #log_odds coef
      log_odds_coef = round(res_sum[j,"Estimate"], digits = 2)
      odds_ratio_coef = round(exp(log_odds_coef),2)
      lower_ci<-round(exp(ci[j,'2.5 %']),2)
      upper_ci<-round(exp(ci[j,'97.5 %']),2)
      pval = round(res_sum[j,"Pr(>|z|)"], digits = 3)
      ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      #r2_coef = round(r2_coefs[idx],3)
      pub_ready = paste0('b=',log_odds_coef, ', p',pval,', OR=',odds_ratio_coef, 
                         ', 95% CI [',lower_ci,',',upper_ci,']')
      pub_ready = unname(cbind(j,pub_ready))
      output = unname(rbind(output,pub_ready))
      idx = idx+1 }
    
  } 
  if(any(class(x)=='multinom')){
    summary_model <- summary(x)
    z_vals <- round(summary_model$coefficients / summary_model$standard.errors,3)
    p_vals <- round(2 * (1 - pnorm(abs(z_vals))),3)
    z_vals[] <- paste0('z = ',z_vals)
    p_vals[] <- paste0(', p = ',p_vals)
    z_vals[] <- paste0(z_vals,p_vals)
    output <- z_vals
  }
  #browser()
  if(any(class(x)=='emmGrid')){
    summary_model <- summary(x)
    idx = 1
    for (j in seq(nrow(summary_model))){
      contrast_name = as.character(summary_model[[j, 1]])
      contrast_name_factor = as.character(summary_model[[j, 2]])
      final_name = paste(contrast_name, contrast_name_factor, collapse = '')
      contrast_estimate = round(summary_model[j,'estimate'],3)
      t = round(summary_model[j,'t.ratio'],3)
      pval = round(summary_model[j,'p.value'],3)
      df = round(summary_model[j,'df'],2)
      pval <- ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      pub_ready = paste0('contrast estimate =',contrast_estimate,
                         ', t(',df,')=',t, ', p',pval)
      pub_ready = unname(cbind(unname(final_name),pub_ready))
      output = unname(rbind(output,pub_ready))
      #browser()
    }
    
  }
  if(any(class(x)=='summary_emm')){
    #browser()
    summary_model <- summary(x)
    idx = 1
    for (j in seq(nrow(summary_model))){
      final_name = as.character(summary_model[[j, 1]])
      contrast_estimate= round(summary_model[[j, 2]],3)
      t = round(summary_model[j,'t.ratio'],3)
      pval = round(summary_model[j,'p.value'],3)
      df = round(summary_model[j,'df'],2)
      pval <- ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      pub_ready = paste0('b =',contrast_estimate,
                         ', t(',df,')=',t, ', p',pval)
      pub_ready = unname(cbind(unname(final_name),pub_ready))
      output = unname(rbind(output,pub_ready))
      #browser()
    }
    
  }
  #browser()
  if(any(class(x)=='sim_slopes')){
    res<-x$slopes
    for (j in rownames(res)){
      estimate <- round(res[j,'Est.'],3)
      pval <- round(res[j,'p'],3)
      pval <- ifelse(pval==0, pval<-"<.001",pval<-paste0('=',pval))
      upper_ci <- round(res[j,'97.5%'],3)
      lower_ci <- round(res[j,'2.5%'],3)
      pub_ready = paste0('simple slope=',estimate, ', p',pval, 
                         ', 95% CI[',lower_ci,',',upper_ci,']' )
      final_name<- res[j,1]
      pub_ready = unname(cbind(unname(final_name),pub_ready))
      output = unname(rbind(output,pub_ready))
      #browser()
    }
  }
  if (all(class(x) == "BFBayesFactor")){
    #browser()
    bf = x@bayesFactor['bf']
    output = paste0('BF10 = ', round(exp(bf),3))
  }
  if (any(class(x) == "mipo")) {
    res_sum <- summary(x)
    res_sum <- res_sum[res_sum$term != "(Intercept)", ]
    fdr_vals <- p.adjust(res_sum$p.value, method = "fdr")
    rownames(res_sum) <- res_sum$term
    
    for (j in rownames(res_sum)) {
      b     <- round(res_sum[j, "estimate"], 3)
      tstat <- round(res_sum[j, "statistic"], 2)
      df    <- round(res_sum[j, "df"], 1)
      pval  <- round(res_sum[j, "p.value"], 3)
      fdr   <- round(fdr_vals[which(rownames(res_sum) == j)], 3)
      
      ifelse(pval == 0, pval <- "<.001", pval <- paste0("=", pval))
      ifelse(fdr  == 0, fdr  <- "<.001", fdr  <- paste0("=", fdr))
      
      pub_ready <- paste0("b=", b, ", t(", df, ")=", tstat, ", p", pval, ", p_FDR", fdr)
      pub_ready <- unname(cbind(j, pub_ready))
      output    <- unname(rbind(output, pub_ready))
    }
  }
  return(output)
}

vjp_demographics_row <- function(df, main_var, grouping_var, pretty_name, 
                                 nominal, central_tendency='mean'){
  
  if (is.character(df[[main_var]])){
    df[[main_var]] = as.factor(df[[main_var]])
  }
  
  if (is.na(central_tendency)){
    central_tendency = 'mean'
  }
  
  n_groups <- length(levels(as.factor(df[[grouping_var]])))  # <-- capture once
  
  if (nominal == FALSE){
    if (central_tendency == 'mean'){
      summary_df <- df %>% 
        group_by(!!as.name(grouping_var)) %>%
        get_summary_stats(!!as.name(main_var), type = "mean_sd")
      summary_df = paste0(round(summary_df$mean, 2), ' (', round(summary_df$sd, 2), ')')
      
      if (n_groups == 1){                        # <-- NEW: single group, no test
        paste_stats  <- ''
        sig_contrasts <- ''
      } else if (n_groups == 2){
        res <- t.test(base::get(main_var) ~ base::get(grouping_var), df, var.equal = TRUE)
        paste_stats  <- pub_ready_stats(res)
        sig_contrasts <- ''
      } else {
        res <- aov(base::get(main_var) ~ base::get(grouping_var), df)
        res_sum <- summary(res)[[1]]
        paste_stats  <- pub_ready_stats(res)
        if (res_sum$`Pr(>F)`[1] < .05){
          post_hoc_res <- TukeyHSD(res)
          post_hoc_res <- data.frame(post_hoc_res[[grouping_var]])
          groups_compared <- rownames(post_hoc_res)
          post_hoc_p_vals <- post_hoc_res$p.adj
          sig_contrasts <- paste0(groups_compared[post_hoc_p_vals < .05], collapse = ', ')
          if (length(sig_contrasts) == 0) sig_contrasts <- ''
        } else {
          sig_contrasts <- ''
        }
      }
    }
    
    if (central_tendency == 'median'){
      summary_df <- df %>% 
        group_by(!!as.name(grouping_var)) %>%
        get_summary_stats(!!as.name(main_var), type = "median")
      summary_df = round(summary_df$median, 2)
      
      if (n_groups == 1){                        # <-- NEW: single group, no test
        paste_stats  <- ''
        sig_contrasts <- ''
      } else if (n_groups == 2){
        res <- coin::median_test(formula(paste0(main_var, ' ~ ', grouping_var)), data = df)
        z_stat <- round(res@statistic@teststatistic, 3)
        p_val  <- round(2 * (1 - pnorm(abs(z_stat))), 3)
        paste_stats  <- paste0('Z=', z_stat, ', p=', p_val)
        sig_contrasts <- ''
      } else {
        res <- coin::median_test(formula(paste0(main_var, ' ~ ', grouping_var)), data = df)
        chi_stat <- res@statistic@teststatistic
        deg_free <- res@statistic@df
        p_val    <- pchisq(chi_stat, deg_free, lower.tail = FALSE)
        paste_stats <- pub_ready_stats(res)
        if (p_val < .05){
          groups_of_interest <- unique(df[[grouping_var]])
          group_pairs <- t(combn(groups_of_interest, 2))
          sig_contrasts <- NULL
          for (jj in seq(nrow(group_pairs))){
            group_one  <- as.character(group_pairs[jj, 1])
            group_two  <- as.character(group_pairs[jj, 2])
            df_subset  <- df[df[[grouping_var]] %in% c(group_one, group_two), ]
            ph_res     <- coin::median_test(formula(paste0(main_var, ' ~ ', grouping_var)),
                                            data = df_subset)
            z_stat     <- ph_res@statistic@teststatistic
            p_val_lower <- pnorm(z_stat, lower.tail = TRUE)
            p_val_upper <- pnorm(z_stat, lower.tail = FALSE)
            p_val       <- round(min(p_val_lower, p_val_upper), 3)
            if (p_val < .05)
              sig_contrasts <- paste0(sig_contrasts, group_one, '-', group_two, '; ')
          }
          if (is.null(sig_contrasts)) sig_contrasts <- ''
        } else {
          sig_contrasts <- ''
        }
      }
    }
    
    demo_row <- cbind(pretty_name, format(t(summary_df), scientific = FALSE),
                      unname(paste_stats), unname(sig_contrasts))
    
  } else if (nominal == TRUE){
    main_factor    <- droplevels(as.factor(df[[main_var]]))
    grouping_factor <- as.factor(df[[grouping_var]])
    counts          <- table(main_factor, grouping_factor)
    group_marginals <- colSums(counts)
    
    summary_df      <- sweep(counts, 2, group_marginals, `/`)
    summary_df[]    <- paste0(round(summary_df * 100, 1), '%')
    row_names       <- rownames(summary_df)
    summary_df      <- cbind(unname(row_names), summary_df)
    padding         <- rbind(rep('', nrow(summary_df)))  # kept original shape
    
    if (n_groups == 1){                          # <-- NEW: single group, no test
      paste_stats    <- ''
      sig_vals_vector <- ''
    } else {
      res         <- chisq.test(counts)
      paste_stats <- pub_ready_stats(res)
      if (!is.na(res$p.value) && res$p.value < .05){
        post_hoc_res      <- chisq.posthoc.test::chisq.posthoc.test(counts)
        post_hoc_res      <- post_hoc_res[post_hoc_res$Value == 'p values', ]
        post_hoc_res_long <- pivot_longer(post_hoc_res,
                                          cols = levels(df[[grouping_var]]),
                                          names_to = 'Group', values_to = 'pvals')
        post_hoc_res_long$pvals <- round(post_hoc_res_long$pvals, 3)
        if (any(post_hoc_res_long$pvals < .05)){
          sig_vals <- data.frame(post_hoc_res_long[post_hoc_res_long$pvals < .05, ])
          sig_vals_vector <- NULL
          for (j in seq(nrow(sig_vals))){
            sig_row <- sig_vals[j, ]
            sig_vals_vector <- paste0(sig_vals_vector,
                                      sig_row$Dimension, ' ', sig_row$Group, '; ')
          }
        } else {
          sig_vals_vector <- 'No Sig. Post Hoc Differences'
        }
      } else {
        sig_vals_vector <- ''
      }
    }
    
    if (length(levels(main_factor)) == 2){
      summary_df <- summary_df[1, , drop = FALSE]
      padding    <- padding[1]
    }
    
    demo_table <- cbind(unname(summary_df), t(padding), t(padding))
    empty_row  <- c(pretty_name, rep("", ncol(summary_df) - 1),
                    unname(paste_stats), sig_vals_vector)
    demo_row   <- rbind(unname(empty_row), demo_table)
  }
  
  return(demo_row)
}

vjp_build_demographics_table <- function(df, demo_vars, grouping_var){
  demo_table = NULL
  #browser()
  if (is.character(df[[grouping_var]])){
    df[[grouping_var]] = as.factor(df[[grouping_var]])
  }
  
  #browser()
  df[[grouping_var]]<-droplevels(df[[grouping_var]])
  for (j in demo_vars){
    #browser()
    if (length(j)==3){
      j = c(j,NULL)
    }
    #browser()
    row<-vjp_demographics_row(df,j[1],grouping_var,j[2],j[3],j[4])
    #browser()
    demo_table <- rbind(demo_table,unname(row))
  }
  
  #create colnames
  demo_df <- as.data.frame(demo_table)
  ns <- df %>% dplyr::filter(!is.na(!!as.name(grouping_var))) %>% 
    dplyr::count(!!as.name(grouping_var))
  group_names<-levels(df[[grouping_var]])
  #browser()
  colnames(demo_df)[-1] <- cbind(t(paste0(group_names,' (n=',ns$n,')')),
                                 "Statistics","Post Hocs")
  return(demo_df)
}
conf_int_proportion <- function(p, n){
  margin <- qnorm(0.975)*sqrt(p*(1-p)/n)
  upper_ci <- round(p + margin,4)*100
  lower_ci <- round(p - margin,4)*100
  output<-paste0(round(p*100,2),'%, 95% CI [',lower_ci,', ',upper_ci,']')
  return(output)
}

get_auc <- function(model){
  library(pROC)
  #browser()
  predicted_probs <- predict(model, type = "response")
  df = model$data
  dv = names(model$model[1])
  #
  roc_obj <- roc(model$model[[dv]], predicted_probs)
  # Print AUC
  auc_ci <- ci.auc(roc_obj)
  ci_lo<-round(auc_ci,3)[1]
  estimate<-round(auc_ci,3)[2]
  ci_hi<- round(auc_ci,3)[3]
  plot(roc_obj, main = paste("AUC =", estimate, ', [',ci_lo,',',ci_hi,']'))
}

plot_three_way <- function(model,
                           facet_var,
                           legend_var,
                           x_axis_var, 
                           x_lab       = NULL,
                           facet_lab   = NULL,
                           y_lab       = NULL,
                           legend_lab  = NULL,
                           drop_legend = FALSE,
                           save_name   = NULL,
                           save_height = NULL,
                           save_width  = NULL) {
  
  if (is.null(x_lab))      x_lab      <- x_axis_var
  if (is.null(y_lab))      y_lab      <- "Predicted outcome"
  if (is.null(facet_lab))  facet_lab  <- facet_var
  if (is.null(legend_lab)) legend_lab <- legend_var
  
  preds <- ggeffects::ggpredict(
    model,
    terms = c(
      x_axis_var,
      paste0(legend_var, " [meansd]"),
      facet_var
    )
  )
  preds <- as.data.frame(preds)
  
  # CI naming differences
  if (!("conf.low" %in% names(preds)) && "conf.lo" %in% names(preds)) {
    preds$conf.low  <- preds$conf.lo
  }
  if (!("conf.high" %in% names(preds)) && "conf.hi" %in% names(preds)) {
    preds$conf.high <- preds$conf.hi
  }
  
  # relabel groups: -1 SD, mean, +1 SD
  lvls <- sort(unique(preds$group))
  preds$group <- factor(
    preds$group,
    levels = rev(lvls),
    labels = rev(c("-1 SD", "mean", "+1 SD"))
  )
  
  # Named Okabe–Ito palette
  cb_vals <- c(
    "-1 SD" = "#0072B2",
    "mean"  = "#E69F00",
    "+1 SD" = "#009E73"
  )
  cb_vals <- cb_vals[levels(preds$group)]
  
  main_plot <- ggplot(preds,
                      aes(x = x,
                          y = predicted,
                          group = group)) +
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high,
          fill = group),
      alpha = 0.20,
      color = NA,
      show.legend = FALSE
    ) +
    geom_line(
      aes(color = group,
          linetype = group),
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = group,
          shape = group),
      size = 2
    ) +
    facet_wrap(~ facet) +
    scale_color_manual(values   = cb_vals) +
    scale_fill_manual(values    = cb_vals, guide = "none") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_shape_manual(values    = c(16, 17, 15)) +
    labs(
      x        = x_lab,
      y        = y_lab,
      color    = legend_lab,
      linetype = legend_lab,
      shape    = legend_lab,
      #title    = paste0(x_lab, " x ", facet_lab, " x ", legend_lab)
    ) +
    theme_minimal()
  
  # ---- DROP LEGEND IF REQUESTED ----
  if (drop_legend) {
    main_plot <- main_plot + theme(legend.position = "none")
  }
  
  if (!is.null(save_name)) {
    ggsave(
      filename = save_name,
      plot     = main_plot,
      bg       = "white",
      height   = save_height,
      width    = save_width,
      dpi      = 500
    )
  }
  
  return(main_plot)
}

plot_two_way <- function(model,
                         legend_var,
                         x_axis_var,
                         x_lab       = NULL,
                         y_lab       = NULL,
                         legend_lab  = NULL,
                         continuous_legend_var = TRUE,
                         drop_legend = FALSE,
                         drop_title  = FALSE,
                         save_name   = NULL,
                         save_height = 4,
                         save_width  = 4) {
  
  if (is.null(x_lab))      x_lab      <- x_axis_var
  if (is.null(y_lab))      y_lab      <- "Predicted outcome"
  if (is.null(legend_lab)) legend_lab <- legend_var
  
  # predicted values with CIs
  preds <- ggeffects::ggpredict(
    model,
    terms = c(
      x_axis_var,
      paste0(legend_var, " [meansd]")
    )
  )
  preds <- as.data.frame(preds)
  
  # handle possible CI naming differences
  if (!("conf.low" %in% names(preds)) && "conf.lo" %in% names(preds)) {
    preds$conf.low  <- preds$conf.lo
  }
  if (!("conf.high" %in% names(preds)) && "conf.hi" %in% names(preds)) {
    preds$conf.high <- preds$conf.hi
  }
  
  # if the legend var is continuous, relabel to -1 SD / mean / +1 SD
  if (continuous_legend_var) {
    lvls <- sort(unique(preds$group))
    preds$group <- factor(
      preds$group,
      levels = rev(lvls),
      labels = rev(c("-1 SD", "mean", "+1 SD"))
    )
  } else {
    preds$group <- factor(preds$group)
  }
  
  # base plot
  main_plot <- ggplot(preds,
                      aes(x = x,
                          y = predicted,
                          group = group)) +
    # CI ribbon (no legend)
    geom_ribbon(
      aes(ymin = conf.low,
          ymax = conf.high,
          fill = group),
      alpha = 0.20,
      color = NA,
      show.legend = FALSE
    ) +
    # lines & points (carry legend)
    geom_line(
      aes(color = group,
          linetype = group),
      linewidth = 0.9
    ) +
    geom_point(
      aes(color = group,
          shape = group),
      size = 2
    ) +
    labs(
      x        = x_lab,
      y        = y_lab,
      color    = legend_lab,
      linetype = legend_lab,
      shape    = legend_lab,
      title    = paste0(x_lab, " x ", legend_lab)
    ) +
    theme_minimal()
  
  # color scales
  if (continuous_legend_var && length(levels(preds$group)) == 3) {
    # Okabe–Ito palette for -1 SD / mean / +1 SD
    cb_vals <- c(
      "-1 SD" = "#0072B2",
      "mean"  = "#E69F00",
      "+1 SD" = "#009E73"
    )
    cb_vals <- cb_vals[levels(preds$group)]
    
    main_plot <- main_plot +
      scale_color_manual(values   = cb_vals) +
      scale_fill_manual(values    = cb_vals, guide = "none") +
      scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
      scale_shape_manual(values    = c(16, 17, 15))
  } else {
    # let ggplot choose colors; still hide separate fill legend
    main_plot <- main_plot +
      scale_fill_discrete(guide = "none")
  }
  
  # drop title / legend if requested
  if (drop_title) {
    main_plot <- main_plot + labs(title = NULL)
  }
  if (drop_legend) {
    main_plot <- main_plot + theme(legend.position = "none")
  }
  
  # save if requested
  if (!is.null(save_name)) {
    ggsave(
      filename = save_name,
      plot     = main_plot,
      bg       = "white",
      height   = save_height,
      width    = save_width,
      dpi      = 500
    )
  }
  
  return(main_plot)
}

cb_scale_color <- function(levels = NULL) {
  
  cb_palette <- c(
    "#000000", # black
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#009E73", # bluish green
    "#F0E442", # yellow
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#CC79A7"  # reddish purple
  )
  
  if (!is.null(levels)) {
    cb_palette <- cb_palette[seq_along(levels)]
    names(cb_palette) <- levels
  }
  
  ggplot2::scale_color_manual(values = cb_palette)
}

cb_scale_shape <- function(levels) {
  # 12 shapes; repeats if you somehow exceed 12
  base_shapes <- c(15, 16, 17, 18, 3, 7, 8, 4, 9, 10, 12, 13)
  
  lvls <- as.character(levels)
  n <- length(lvls)
  
  shps <- rep(base_shapes, length.out = n)
  names(shps) <- lvls
  
  ggplot2::scale_shape_manual(values = shps)
}

vjp_plot_ttest <- function(data, dv, group, paired = FALSE,
                           point_alpha = 0.5, jitter_width = 0.10,
                           x_label = NULL, y_label = NULL,
                           title = NULL, plot_save = NULL,
                           stats_type = "ttest",
                           bf_direction = "10",
                           ylims = NULL,
                           draw_line = TRUE) {
  
  library(ggplot2)
  library(ggtext)
  library(dplyr)
  library(rlang)
  
  dv    <- ensym(dv)
  group <- ensym(group)
  dv_name <- as_string(dv)
  g_name  <- as_string(group)
  
  # Validate arguments
  stats_type   <- match.arg(stats_type,   choices = c("ttest", "bayes", "both", "none"))
  bf_direction <- match.arg(bf_direction, choices = c("10", "01"))
  
  # Drop NAs
  df <- data %>%
    filter(!is.na(.data[[dv_name]]),
           !is.na(.data[[g_name]]))
  
  # Group levels (used by paired branches)
  groups <- levels(factor(df[[g_name]]))
  
  # --- t-test ---
  ttest_label <- NULL
  if (stats_type %in% c("ttest", "both")) {
    t_out <- if (paired) {
      t.test(
        x = df[[dv_name]][df[[g_name]] == groups[1]],
        y = df[[dv_name]][df[[g_name]] == groups[2]],
        paired = TRUE
      )
    } else {
      t.test(df[[dv_name]] ~ df[[g_name]])
    }
    
    ttest_label <- sprintf(
      "t(%s) = %.2f, p = %.3g",
      format(t_out$parameter, digits = 4),
      unname(t_out$statistic),
      t_out$p.value
    )
  }
  
  # --- Bayes Factor ---
  bf_label <- NULL
  if (stats_type %in% c("bayes", "both")) {
    library(BayesFactor)
    
    bf_out <- if (paired) {
      ttestBF(
        x = df[[dv_name]][df[[g_name]] == groups[1]],
        y = df[[dv_name]][df[[g_name]] == groups[2]],
        paired = TRUE
      )
    } else {
      ttestBF(formula = as.formula(paste(dv_name, "~", g_name)), data = df)
    }
    
    bf_value <- extractBF(bf_out)$bf[1]
    
    # Flip to BF01 if requested
    bf_value <- if (bf_direction == "01") 1 / bf_value else bf_value
    bf_sub   <- if (bf_direction == "01") "01" else "10"
    
    bf_label <- if (bf_value >= 1000 | bf_value < 0.001) {
      sprintf("BF<sub>%s</sub> = %.2e", bf_sub, bf_value)
    } else {
      sprintf("BF<sub>%s</sub> = %.2f", bf_sub, bf_value)
    }
  }
  
  # --- Build annotation label ---
  stats_label <- switch(stats_type,
                        "ttest" = ttest_label,
                        "bayes" = bf_label,
                        "both"  = paste0(ttest_label, "<br>", bf_label),
                        "none"  = NULL
  )
  
  # --- Annotation position: horizontal center, 99% up ---
  x_pos  <- (1 + length(groups)) / 2
  y_range <- range(df[[dv_name]], na.rm = TRUE)
  y_pos  <- y_range[1] + 0.99 * diff(y_range)
  
  # Fall back to variable names if labels not provided
  x_label <- if (!is.null(x_label)) x_label else g_name
  y_label <- if (!is.null(y_label)) y_label else dv_name
  
  p <- ggplot(df, aes(x = .data[[g_name]], y = .data[[dv_name]])) +
    geom_violin(
      trim = FALSE, scale = "width",
      alpha = 0.3, fill = "grey70", color = NA
    ) +
    geom_jitter(width = jitter_width, alpha = point_alpha, size = 1.6) +
    stat_summary(
      fun.data = mean_cl_normal, geom = "errorbar",
      width = 0.15, color = "red", linewidth = 1
    ) +
    stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
    theme_classic(base_size = 16) +
    labs(x = x_label, y = y_label, title = title)
  
  if (draw_line) {
    p <- p + stat_summary(fun = mean, geom = "line", aes(group = 1), color = "red", linewidth = 0.8)
  }
  
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  if (!is.null(stats_label)) {
    p <- p + annotate(
      "richtext",
      x = x_pos, y = y_pos,
      label = stats_label,
      hjust = 0.5, vjust = 0.5,
      lineheight = 1.4,
      color = "black",
      fill = NA, label.color = NA
    )
  }
  
  if (!is.null(plot_save)) {
    ggsave(plot_save, plot = p, width = 5, height = 4)
  }
  
  p
}

library(lme4)
library(performance)

compute_reliability <- function(df, id_col = "src_subject_id", t1_col, t2_col,
                                covariate = NULL,
                                x_label = t1_col, y_label = t2_col,
                                title = "Test-Retest Reliability",
                                plot_path = NULL) {
  
  v1 <- df[[t1_col]]
  v2 <- df[[t2_col]]
  
  valid <- complete.cases(v1, v2)
  v1 <- v1[valid]; v2 <- v2[valid]
  
  # ── Residualize against covariate if provided ─────────────────────────────
  adj_note <- ""
  if (!is.null(covariate)) {
    cov_vec  <- df[[covariate]][valid]
    v1       <- residuals(lm(v1 ~ cov_vec))
    v2       <- residuals(lm(v2 ~ cov_vec))
    adj_note <- paste0(" (adj. for ", covariate, ")")
  }
  # ─────────────────────────────────────────────────────────────────────────
  
  if (length(v1) < 3) stop("Not enough complete cases")
  if (sd(v1) == 0 | sd(v2) == 0) stop("No variance in one of the columns")
  
  r       <- cor(v1, v2)
  r_ci    <- cor.test(v1, v2)$conf.int
  
  icc_out <- irr::icc(cbind(v1, v2), model = "twoway", type = "agreement", unit = "single")
  icc_val <- icc_out$value
  icc_lb  <- icc_out$lbound
  icc_ub  <- icc_out$ubound
  
  message(
    "Test-retest reliability for ", t1_col, adj_note, ": ",
    "r = ", round(r, 2), " [", round(r_ci[1], 2), ", ", round(r_ci[2], 2), "], ",
    "ICC(2,1) = ", round(icc_val, 2), " [", round(icc_lb, 2), ", ", round(icc_ub, 2), "], ",
    "n = ", length(v1)
  )
  
  results <- tibble::tibble(
    variable      = t1_col,
    covariate     = if (!is.null(covariate)) covariate else NA_character_,
    r             = round(r, 2),
    r_95ci        = paste0("[", round(r_ci[1], 2), ", ", round(r_ci[2], 2), "]"),
    ICC_agreement = round(icc_val, 2),
    ICC_95ci      = paste0("[", round(icc_lb, 2), ", ", round(icc_ub, 2), "]"),
    n             = length(v1)
  )
  
  plot_data <- data.frame(t1 = v1, t2 = v2)
  p <- ggplot(plot_data, aes(x = t1, y = t2)) +
    geom_point(shape = 21, fill = NA, color = "black", size = 2, stroke = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#2196F3") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste0(title, adj_note),
         subtitle = paste0(
           "r = ", round(r, 2), " [", round(r_ci[1], 2), ", ", round(r_ci[2], 2), "]",
           "  |  ICC(2,1) = ", round(icc_val, 2), " [", round(icc_lb, 2), ", ", round(icc_ub, 2), "]",
           "  |  n = ", length(v1)),
         x = paste0(x_label, adj_note), 
         y = paste0(y_label, adj_note)) +
    theme_classic(base_size = 14)
  
  print(p)
  
  if (!is.null(plot_path)) {
    ggsave(filename = plot_path, plot = p, width = 6.5, height = 5)
    message("Saved plot to: ", plot_path)
  }
  
  return(invisible(results))
}


