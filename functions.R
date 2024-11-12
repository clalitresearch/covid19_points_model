#############################################################################################
print_total = function(.data, msg) {
#############################################################################################
  total = .data %>%
    dplyr::summarise(total = n()) %>% 
    pull(total)
  message(msg,": total:",total);
  .data
}


# this function finds the minimal significant and positive coefficient which will become the reference coef. for transforming coefficients into risk points
# e.g. if the coefficient for heart disease is 0.25 than all coefficients will be divided by 0.25

get_normalization_factor = function(glm_fit) {
  broom::tidy(glm_fit, conf.int=F) %>% 
    filter(!grepl('age', term, ignore.case = TRUE)) %>% 
    mutate(abs_estimate = abs(estimate)) %>% 
    select(term, abs_estimate) %>% 
    arrange(abs_estimate) %>% 
    head(1) %>% 
    as.list()
}

get_scores = function(glm_fit, ref_var_value, simplification_factor = 1.0, normalization_factor = 1.0, calc_initial_risk_points = F) {
  stopifnot (simplification_factor >= 1.0)
  multiplication_factor = 1.0/normalization_factor
  results = broom::tidy(glm_fit, conf.int=F) %>% 
    filter(!grepl('intercept', term, ignore.case = TRUE)) %>% 
    mutate(scores = multiplication_factor *(estimate)) %>% 
    mutate(divided = round(abs(scores) / simplification_factor)) %>% 
    mutate(verified = if_else(divided ==0, 1, divided)) %>% 
    mutate(simplified_scores= sign(scores)*verified ) %>% 
    mutate(original_scores = round(scores) ) %>% 
    select(term, original_scores,simplified_scores) 
  
  if (simplification_factor > 1.0){
    final =  (results %>%  select(term, simplified_scores) %>%  dplyr::rename(scores = simplified_scores))
  } else {
    final = (results %>%  select(term, original_scores) %>%  dplyr::rename(scores = original_scores))
  }
  if (calc_initial_risk_points) {
    lowest_possibe = get_lowest_possible_score(final)
    final = final %>% add_row(term = 'Inital Risk Points', scores = -lowest_possibe, .before = 1)
  }
  setNames(final$scores, final$term)
}

get_lowest_possible_score = function(df) {
  
  part1 = df %>% filter(grepl('months_from_last', term)) %>% dplyr::summarise(v=min(scores)) %>% pull(v)
  part2 = df %>%  filter(term == 'sw_previous_infection') %>% pull(scores)
  return (part1 + part2)
}
  


binom_ci <- function(row) {
  numer = row['numer']
  denom = row['denom']
  if (denom == 0) return(NA)
  this_test <- binom.test(x = numer, n = denom)
  estimate <- this_test$estimate[[1]]
  lower <- this_test$conf.int[1]
  upper <- this_test$conf.int[2]
  list(estimate=100*estimate, lower=100*lower, upper=100*upper)
}


# ----------------------------------------------------
add_classifier_metrics = function(df, y, y_hat_metrics) {
# ----------------------------------------------------
  stopifnot(length(y)==length(y_hat_metrics))
  
  df %>% mutate (
    TP    = lapply(score_threshold, function(x) sum(y_hat_metrics[y == 1] > x))  %>%  unlist(),
    TN    = lapply(score_threshold, function(x) sum(y_hat_metrics[y == 0] <= x)) %>%  unlist(),
    numer = lapply(score_threshold, function(x) sum(y[y_hat_metrics == x]))      %>%  unlist(),
    denom = lapply(score_threshold, function(x) length(y[y_hat_metrics == x]))   %>%  unlist()) %>%
    mutate(`%pop_>=` = 100*(1-cumsum(denom)/N + denom/N)) %>% 
    dplyr::mutate(
      FN = sum(y) - TP,
      FP = N - sum(y) - TN,
      Sensitivity = TP / (TP + FN),
      FPR = FP / (FP + TN),
      specificity = TN / (TN + FP),
      PPV = TP / (TP + FP),
      NPV = TN / (TN + FN),
      Lift = (TP / (TP + FN)) / ((TP + FP ) / N),
      predicted_positives = TP + FP,
      ppcr  = predicted_positives / N
      
    ) #%>%     select (-numer, -denom, -predicted_positives)
  
}


# -------------------------------------------------------------------------------
get_perfromace_table = function(all_df, n_digits = 5, add_observed_risk = F, add_extra_info = F, score_capping_quantiile = NA) {
# -------------------------------------------------------------------------------
  # auc_df = a sperate DataFrame holding the higer risk patients to show AUC for the more relevant population
  auc_df = all_df %>% 
    filter(y_hat >=0.1)
  
  y_hat = all_df$y_hat
  y     = all_df$y
  score = all_df$score
  if (score_capping_quantiile > 0 ) {
    capping_score = quantile(score, score_capping_quantiile, na.rm = T)
    print (glue("capping score ={capping_score} "))
    score = case_when(score > capping_score ~ capping_score, T ~ score)
  }
  
  performance_score =  tibble::tibble(score_threshold = sort(unique(score)),
                                      N = length(y))
  performance_score = add_classifier_metrics(performance_score, y, score) 
  
  performance_score = performance_score %>% 
    mutate( min_pred_risk = lapply(score_threshold, function(x) round(100*min(y_hat[score == x]),3)) %>%  unlist(),
            max_pred_risk = lapply(score_threshold, function(x) round(100*max(y_hat[score == x]),3)) %>%  unlist())
  if (add_observed_risk) {
    print ("calculating observed risk")
    results = apply(performance_score,1,binom_ci) #results is a list of lists`
    performance_score$obs.risk       = as.numeric(unlist(results)[attr(unlist(results),"names")=="estimate"])
    performance_score$obs.risk.lower = as.numeric(unlist(results)[attr(unlist(results),"names")=="lower"])
    performance_score$obs.risk.upper = as.numeric(unlist(results)[attr(unlist(results),"names")=="upper"])
  }
  performance_score = performance_score %>% 
    rename(points = score_threshold)  %>% 
    mutate(across(everything(), ~round(.x, n_digits))) 
  auc = pROC::auc(y, score)
  partial_auc = pROC::auc(auc_df$y, auc_df$score)
  min_negative_points=min(performance_score$points)
  

  performance_score = performance_score %>% 
    mutate(points = points + abs(min_negative_points))
  
  if (add_extra_info) {
    return (list(metrics=performance_score, auc=auc, partial_auc = partial_auc,min_negative_points=min_negative_points ))
  } else {
    return (performance_score)
  }
}




# -------------------------------------------------------------------------------
get_rtichoke_performance_table = function(df, add_observed_risk, glm_fit, scores, n_digits = 5, score_capping_quantiile = NA) {
# -------------------------------------------------------------------------------
  #print ("get_rtichoke_performance_table")
  #write.csv(df, glue("./debug/df.csv"))
  if (is.null(scores)) {
    probs  = predict(glm_fit, newdata = df, type = c("response"))
    y = df$sw_outcome
  } else {
    stats_df = calc_scores(df, glm_fit, scores)
    # we just convert the scores to numbers between [0.1] preserving the ranking. these are NOT probabilities and ofcourse are NOT calibrated!
    probs  =stats_df$score/1000.0
    y     = stats_df$y
    
  }
  rtichoke::create_performance_table(probs = probs, real = y, by=0.0002, enforce_riskpercentiles_symmetry = TRUE) %>% 
    rename(ppcr = percentpositives,Sensitivity = TPR, PPV=ppv, Lift=lift, predicted_positives = positives)
  
}

# -------------------------------------------------------------------------------
calc_metrics = function(df, add_observed_risk, glm_fit, raw_scores, n_digits=3, score_capping_quantiile = NA) {
# -------------------------------------------------------------------------------
  stats_df = calc_scores(df, glm_fit, raw_scores)
  metrics = get_perfromace_table(stats_df,n_digits, add_observed_risk, score_capping_quantiile= score_capping_quantiile)
  metrics
}

# -------------------------------------------------------------------------------
calc_scores = function(df, glm_fit, raw_scores) {
# -------------------------------------------------------------------------------
  all_predictors = tidy(glm_fit) %>% pull (term)
  predictors = all_predictors[-1]
  scores = rowSums(data.frame(mapply(`*`, df %>% select(all_of(predictors)), raw_scores,SIMPLIFY=FALSE)))
  y_hat  = predict(glm_fit, newdata = df, type = c("response"))
  all_df = tibble::tibble(score = scores,
                          y = df$sw_outcome ,
                          y_hat = y_hat)

  all_df = all_df %>% 
  group_by(score) %>% 
  dplyr::mutate(n=n()) %>% 
  ungroup()
  
  
  return (all_df)
}



# Corresponding Points model:
# -------------------------------------------------------------------------------
get_points_model = function(raw_scores, min_possible_score,prevent_negative_scores=F) {
# -------------------------------------------------------------------------------
  
  scores_df = data.frame(lapply(raw_scores, type.convert), stringsAsFactors=FALSE) %>% 
    pivot_longer(everything(), names_to = "var", values_to = "points")
  if (prevent_negative_scores) {
    scores_df = scores_df %>% 
      add_row(var = 'age_group10_18_30', points = 0, .before = 1) %>% 
      mutate(points = case_when(grepl('age_group10', var) ~ points + abs(min_possible_score), T ~ points))
  }

  scores_df %>% 
  mutate(var = stringr::str_replace(var, "age_group10", "Age")) %>% 
    mutate(var = stringr::str_replace(var, "plus", "+")) %>% 
    mutate(var = stringr::str_replace(var, "_cat", "")) %>% 
    mutate(var = stringr::str_replace(var, "sw_", "")) 
  
}
