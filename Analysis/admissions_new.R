####################################################
### TIME TRENDS IN ADMISSIONS TO INPATIENT CARE ###
###################################################

# loading packages
library(tidyverse) # data manipulation
library(ggplot2) # visualizations
library(segmented) # segmented (piece-wise) regression
library(purrr) # functional programming
library(MASS) # negative binomial regression
library(ciTools) # confidence intervals
library(AER) #dispersion test

# loading the data for yearly admissions to inpatient care stratified by diagnostic entity, sex and age-group
agg_data_adm <- read_csv("C:/Users/hp/Desktop/NUDZ/Time_Trends_Hospitalizations/Data/agg_data_f_groups_number_of_hosp.csv")

agg_data_adm <- agg_data_adm[order(agg_data_adm$year_adm),] %>% 
  mutate(across(.cols = c(sex, f_diag, age_group), ~ as.factor(.x)), # changing the variable type to factor for sex, age group and diagnostic group
         f_diag_stand = f_diag_abs/pop_size*100000, # calculating the admission rates
         year = year_adm,
         date_seq = as.numeric(year - 1994)) %>% # adding time variable (date_seq, coded as 0 to 21)
  filter(f_diag %in% c("F[0-6]|G30", "F00|G30",  "F10", "F1[1-9]", "F20", "F3[2-3]", "F4", "F50", "F60"), # selecting diagnoses
         !(age_group %in% c("1–4", "5–9", "10–14")),  # selecting adult population only (i.e. 15+ y.)
         sex != "total", 
         year != 2016) %>% # filtering out observations for the last year due to distortions
  mutate(age_cat = as.factor(case_when(age_group %in% c("15–19", "20–24", "25–29") ~ "adolescents & emerging adults", # adding new age categories
                                       age_group %in% c("65–69", "70–74", "75–79", "80+") ~ "seniors",
                                       TRUE ~ "adults"))) 



#################
### MODELING ###
################

# creating functions that prevent the code from stopping the execution upon error call:
# if no breakpoint/segm. model is identified, NA is returned
poss_segm <- possibly(.f = selgmented,
                      otherwise = NA)

poss_pred_segm <- possibly(.f = predict.segmented,
                           otherwise = NA)

poss_broken_fit <- possibly(.f = broken.line,
                            otherwise = NA)

poss_apc <- possibly(.f = slope,
                     otherwise = NA)

poss_aapc <- possibly(.f = aapc,
                      otherwise = c(NA, NA, NA))

poss_predict.segm <- possibly(.f = predict.segmented,
                              otherwise = NA)


# writing a function that prevent ifelse statement from turning factor into integer and preserves factor levels
safe.ifelse <- function(cond, yes, no) {
  class.y <- class(yes)
  if (class.y == "factor") {
    levels.y = levels(yes)
  }
  X <- ifelse(cond,yes,no)
  if (class.y == "factor") {
    X = as.factor(X)
    levels(X) = levels.y
  } else {
    class(X) <- class.y
  }
  return(X)
}

# defining a wrapper for glm.nb() that composes a custom call directly on the dataset,
# which is then captured by glm.nb() as an un-evaluated expression
glm.nb2 <- function(data, ...)
  eval(rlang::expr(glm.nb(data=!!rlang::enexpr(data), !!!list(...))) )

# running a series of models for each combination of diagnostic entity (diag. group/specific diagnosis) and sex (total/female/male)
models_adm <- agg_data_adm %>%
  filter(sex != "total") %>% 
  group_by(f_diag, sex, age_cat) %>%
  nest() %>%
  mutate(map_pois = map(.x = data,
                        ~ glm(f_diag_abs ~ date_seq + offset(log(pop_size/100000)), # Poisson models (pois)
                              family = "poisson",
                              data = .x)),
         dispersion = map_lgl(.x = map_pois,
                              ~ dispersiontest(.x)$p.value < 0.05), # dispersion test for Poisson models: all suggest over-dispersion issues
         map_nb  = map(data,
                       glm.nb2,
                       f_diag_abs ~ date_seq + offset(log(pop_size/100000))), # negative binomial models (nb)
         AIC = map2_lgl(.x = map_pois, # Aikaike Information criterion: Poisson vs. negative binomial
                        .y = map_nb,
                        ~ AIC(.x) > AIC(.y)),
         LR_test = map2_lgl(.x = map_pois, # Likelihood ratio test: poisson vs. negative binomial
                            .y = map_nb,
                            ~ pchisq(2 * (logLik(.y) - logLik(.x)), df = 1, lower.tail = FALSE) < 0.05),
         segm_lm = map(.x = map_nb, ~ poss_segm(olm = .x, # SEGMENTED regression
                                                type = "bic", # using Bayes Information Criterion to select the best fitting model
                                                Kmax = 2, # maximal number of breakpoints = 2, i.e. 0, 1 or 2 breakpoints can be extracted
                                                control = seg.control(fix.npsi = TRUE, 
                                                                      n.boot= 100, # number of bootstrap samples used in the bootstrap restarting algorithm (Wood, S. N. (2001) Minimizing model fitting objectives that contain spurious local minima by bootstrap restarting. Biometrics 57, 240–244.)
                                                                      seed = 2023))), 
         segm_lm = ifelse(is.na(segm_lm), map_nb, segm_lm)) # replacing missing values for segmented models (i.e. when 0 breakpoints identified) with values from negative binomial model


# as assessed by AIC, in 92% of cases provides a negative binomial model a better fit than a Poisson model
sum(models_adm$AIC)/length(models_adm$AIC) * 100

# similarly, likelihood ratio test points towards a better model fit for negative binomial models in 90%
sum(models_adm$LR_test)/length(models_adm$LR_test) * 100

# dispersion issues arise in 87% of the models
sum(models_adm$dispersion)/length(models_adm$dispersion) * 100

# -> altogether, this speaks for the use of negative binomial models, which allow accounting for over-dispersion


###################################
### EXTRACTING MODEL PARAMETERS ###
###################################

models_adm <- models_adm %>% 
  mutate(apc = map(.x = segm_lm,
                   ~ poss_apc(.x, APC = TRUE)), # annual percent change (APC) based on segmented regression
         aapc = map(.x = segm_lm,
                    ~ poss_aapc(.x, exp.it = TRUE) * 100), # average annual percent change (AAPC) based on segmented regression
         aapc_est_nb = map(.x = map_nb,
                           ~ round(unlist((exp(.x$coef[2]) - 1)*100),2)), # average annual percent change (AAPC) based on linear regression
         aapc_CI_nb = map(.x = map_nb, # confidence intervals for AAPC based on linear regression (based on robust "sandwich" standard errors)
                          ~ paste0("(", round((as.numeric(confint(.x, parm = "date_seq")[[1]]) * 100),2), ",", round((as.numeric(confint(.x, parm = "date_seq")[[2]]) * 100),2), ")")),
         aapc_est = map(.x = aapc,
                        ~ round(.x[[1]], 2)), # parameters estimate for AAPC based on segmented regression
         aapc_est = map2(.x = aapc_est,
                         .y = aapc_est_nb,
                         ~ ifelse(is.na(.x), .y, .x)), # if no breakpoint identified, AAPC takes value from linear regression, otherwise from segmented regression
         
         aapc_CI = map(.x = aapc,
                       ~ paste0("(", round(unlist(.x[[2]]), 2), ",", round(unlist(.x[[3]]), 2),")")), # confidence intervals for AAPC based on segmented regression
         aapc_CI = map2(.x = aapc_CI,
                        .y = aapc_CI_nb,
                        ~ ifelse(.x == "(NA,NA)", .y, .x)), # if no breakpoint identified, CIs take value from linear regression
         apcs = map(.x = apc,
                    ~ as.data.frame(.x[[1]])[[1]]), # exctracting individual APCs (APC 1-3)
         apc1 = map(.x = apcs,
                    ~ round(ifelse(is.na(.x[[1]]), NA, .x[[1]]), 2)),
         apc2 = map(.x = apcs,
                    ~ round(ifelse(is.na(.x[[1]]), NA, .x[[2]]), 2)),
         apc3 = map(.x = apcs,
                    ~ round(ifelse(length(.x) <= 2, NA, .x[[3]]), 2)),
         apc1_CI = map(.x = apc, # extracting CIs for individual APCs
                       ~ ifelse(is.na(.x), NA, paste0("(",round(as.numeric(as.data.frame(.x[[1]])["slope1","CI(95%).l"]), 2),",",round(as.numeric(as.data.frame(.x[[1]])["slope1","CI(95%).u"]), 2),")"))),
         apc2_CI = map(.x = apc,
                       ~ ifelse(is.na(.x), NA, paste0("(",round(as.numeric(as.data.frame(.x[[1]])["slope2","CI(95%).l"]), 2),",",round(as.numeric(as.data.frame(.x[[1]])["slope2","CI(95%).u"]), 2),")"))),
         apc3_CI = map(.x = apc,
                       ~ ifelse(nrow(as.data.frame(.x)) <= 2, NA, paste0("(",round(as.numeric(as.data.frame(.x[[1]])["slope3","CI(95%).l"]), 2),",",round(as.numeric(as.data.frame(.x[[1]])["slope3","CI(95%).u"]), 2),")"))),
         breakpoints = map(.x = segm_lm,
                           ~ if(nrow(as_tibble(.x$psi)) == 1) {add_row(as_tibble(.x$psi))} else {as_tibble(.x$psi)}), # extracting breakpoints from segmented regression
         breakpoints = map(.x = breakpoints,
                           ~ if(nrow(.x) == 0) {tibble(Initial = c(NA, NA), Est. = c(NA, NA), St.Err = c(NA, NA))} else {.x}),  # if no or just one breakpoint identified, add missing values (NAs)
         breakpoint1_raw_seq = map_dbl(.x = breakpoints, 
                                       ~ unlist(.x[1,"Est."])), # adding identified breakpoints
         breakpoint1 = map(.x = breakpoint1_raw_seq, 
                           ~ floor(.x + 1994)), # round down to the closest integer (i.e. to whole year)
         breakpoint2_raw_seq = map_dbl(.x = breakpoints,
                                       ~ unlist(.x[2,"Est."])),
         breakpoint2 = map(.x = breakpoint2_raw_seq,
                           ~ floor(.x + 1994)),
         time_period1 = map(.x = breakpoint1, 
                            ~ if_else(.x == "NA", as.character(.x), paste("1994 -", .x))), # time periods for APCs based on breakpoints
         time_period2 = map2(.x = breakpoint1,
                             .y = breakpoint2,
                             ~ gsub("- NA", "- 2015", if_else(.x == "NA", as.character(.x), paste(.x, .y, sep = " - ")))),
         time_period2 = map(.x = time_period2, 
                            ~ ifelse(str_detect(.x, "NA"), NA, .x)),
         time_period3 = map(.x = breakpoint2,
                            ~ if_else(.x == "NA", as.character(.x), paste(.x, "- 2015"))))



###################################
### PREPARING DATA FOR PLOTTING ###
###################################

# creating new dataframes
models_adm <- models_adm %>%
  mutate(new_data = map(.x = data,
                        ~ .x), 
         across(.cols = c(breakpoint1_raw_seq, breakpoint2_raw_seq),
                ~ ifelse(is.na(.x), -1, .x)),
         new_data = map2(.x = new_data,
                         .y = breakpoint1_raw_seq,
                         ~ add_row(.x, pop_size = .x[.x$date_seq == floor(.y),]$pop_size)), # adding the breakpoint-year population estimates (1. breakpoint)
         new_data = map2(.x = new_data,
                         .y = breakpoint1_raw_seq,
                         ~ .x %>% mutate(date_seq = ifelse(is.na(date_seq), .y, date_seq), # adding raw (1 or 2) time breakpoints derived from segmented regression to the time variable
                                         age_group = safe.ifelse(is.na(.x$age_group), as.character(.x[.x$date_seq == floor(.y),]$age_group), as.character(.x$age_group)))), # adding age variable
         new_data = map2(.x = new_data,
                         .y = breakpoint2_raw_seq,
                         ~ add_row(.x, pop_size = .x[.x$date_seq == floor(.y),]$pop_size)), # adding the breakpoint-year population estimates (1. breakpoint)
         new_data = map2(.x = new_data,
                         .y = breakpoint2_raw_seq,
                         ~ .x %>% mutate(date_seq = ifelse(is.na(date_seq), .y, date_seq), # adding raw (1 or 2) time breakpoints derived from segmented regression to the time variable
                                         age_group = ifelse(is.na(age_group), na.omit(unique(.x$age_group)), .x$age_group), # adding age variable
                                         year = date_seq + 1994)), # adding calendar year time variable (e.g. 2003.13)
         plot_data = map2(.x = new_data,
                          .y = segm_lm,
                          ~ .x %>% mutate(pred = exp(poss_predict.segm(.y, .x)))), # predicted values based on segmented regression
         pred_nb = map(.x = map_nb,
                       ~ fitted(.x)), # predicted values based on linear regression
         plot_data = map2(.x = plot_data,
                          .y = pred_nb,
                          ~ .x %>% mutate(pred = ifelse(is.na(pred), .y, pred))), # fiinal predicted values: takes values from linear model if no breakpoint detected
         plot_data = map2(.x = plot_data,
                          .y = segm_lm,
                          ~ add_ci(.x, .y, names = c("lcb", "ucb"))), # adding confidence intervals
         plot_data = map(.x = plot_data,
                         ~ .x %>% mutate(pred = pred/pop_size*100000, # computing fitted values for the whole strata
                                         lcb = lcb/pop_size*100000,
                                         ucb = ucb/pop_size*100000)))


####################
### DESCRIPTIVES ###
####################

# unnesting the dataframes
data_plotting_adm <- models_adm %>%
  dplyr::select(plot_data) %>%
  unnest() %>%
  
  # adding labels to diagnoses
  mutate(F_group = as.factor(case_when(f_diag == "F3[2-3]" ~ "depression",
                                       f_diag == "F00|G30" ~ "dementia in Alzheimer disease",
                                       f_diag == "F10" ~ "AUD",
                                       f_diag == "F1[1-9]" ~ "DUD",
                                       f_diag == "F20" ~ "schizophrenia",
                                       f_diag == "F4" ~ "anxiety disorders",
                                       f_diag == "F50" ~ "eating disorders",
                                       f_diag == "F60" ~ "specific personality disorders",
                                       f_diag == "F[0-6]|G30" ~ "any mental disorder"))) %>%
  mutate(across(.cols = c("lcb", "ucb"), 
                ~ ifelse(!year %% 1 == 0, NA, .x))) # remove confidence intervals around breakpoints

# age- and sex-specific admission rates for years 1994 and 2015 to illustrate the time trend
hosp_rates_94_15 <- data_plotting_adm %>%
  filter(year %in% c(1994, 2015),
         f_diag %in% c("F[0-6]|G30", "F00|G30",  "F10", "F1[1-9]", "F20", "F3[2-3]", "F4", "F50", "F60")) %>%
  dplyr::select(year,f_diag,sex,age_cat,pred,lcb,ucb) %>%
  mutate(across(c(pred, lcb, ucb),
                ~ round(.x, digits = 2))) %>%
  mutate(pred = paste0(pred, " (", lcb, ",", ucb, ")")) %>% 
  distinct() %>% 
  dplyr::select(year,f_diag,sex,age_cat,pred) %>% 
  pivot_wider(names_from = c(sex, year),
              values_from = c(pred)) %>%
  arrange(age_cat)




#####################
### MAIN RESULTS ###
####################

results_table_adm <- models_adm %>%
  mutate(across(.cols = c(apc1_CI, apc2_CI, apc3_CI),
                ~ na_if(.x, "(NULL, NULL)"))) %>%
  mutate(AAPC = paste(aapc_est, aapc_CI, sep = " "),
         APC1 = paste(apc1, apc1_CI, sep = " "),
         APC2 = paste(apc2, apc2_CI, sep = " "),
         APC3 = paste(apc3, apc3_CI, sep = " ")) %>% 
  dplyr::select(f_diag, sex, AAPC, APC1, time_period1, APC2, time_period2, APC3, time_period3) %>%
  mutate(across(.cols = c(APC1, APC2, APC3),
                ~ na_if(.x, "NA NA")),
         across(starts_with("time_"),
                ~ as.character(.x))) %>%
  pivot_wider(names_from = sex, values_from = c(AAPC, APC1, time_period1, APC2, time_period2, APC3, time_period3)) %>% 
  arrange(age_cat) %>% 
  dplyr::select(f_diag, age_cat, ends_with("_female"), ends_with("_male")) 

main_results_table_adm <- results_table_adm %>% 
  dplyr::select(f_diag, age_cat, starts_with("AAPC")) %>% 
  pivot_wider(names_from = age_cat, values_from = c(AAPC_female, AAPC_male)) %>% 
  dplyr::select(f_diag, ends_with("emerging adults"), ends_with("_adults"), ends_with("seniors")) 

append_results_table_adm <- results_table_adm %>% 
  dplyr::select(f_diag, age_cat, starts_with(c("APC1", "time_period1", "APC2", "time_period2", "APC3", "time_period3"))) %>%
  dplyr::select(f_diag, age_cat, ends_with(c("female", "male")))

write.csv(main_results_table_adm, "C:/Users/hp/Desktop/NUDZ/Time_Trends_Hospitalizations/Results/Tables/Results_admissions.csv")
write.csv(append_results_table_adm, "C:/Users/hp/Desktop/NUDZ/Time_Trends_Hospitalizations/Results/Tables/Results_admissions_Append.csv")





#####################
### VISUALIZATION ###
#####################

# selecting colorblind friendly colors
my_colors3   <- c("seniors" = "#56B4E9", 
                  "adults" = "#009E73", 
                  "adolescents & emerging adults" = "#CC79A7")

# selecting line-types
my_linetype   <- c("seniors" = "longdash",
                   "adults" = "dotdash",
                   "adolescents & emerging adults" = "dashed")


# segmented TIME TRENDS in admission rates by sex, age category and psychiatric condition
main_plot_adm <- data_plotting_adm %>%
  arrange(year) %>%
  ggplot(aes(x = year,
             color = age_cat)) +
  geom_line(aes(y = pred),
            size = 1.5) +
  geom_line(aes(y = lcb, color = age_cat, linetype = age_cat), size = 0.5, show.legend = F) +
  geom_line(aes(y = ucb, color = age_cat, linetype = age_cat), size = 0.5, show.legend = F) +
  scale_color_manual(values = my_colors3) +
  scale_linetype_manual(values = my_linetype) +
  labs(x = "Years", y = "Admissions per 100 000 person-years",
       title ="Age- and sex-specific temporal trends in admissions to inpatient care",
       subtitle = "from 1994 to 2015",
       color = "Age category") +
  theme(plot.title = element_text(size = 12,
                                  face = "bold")) +
  facet_wrap(~ F_group + sex,
             scales = "free",
             ncol = 2) +
  theme_bw()

ggsave(filename = "main_plot_adm.eps",
       path = "C:/Users/hp/Desktop/NUDZ/Time_Trends_Hospitalizations/Results/Graphs", 
       width = 21, 
       height = 29.7, 
       device='eps', 
       dpi=700)


# separate plots for each stratum (psychiatric condition x age category x sex)
main_diagnoses <- c("depression", "dementia", "DUD", "AUD", "schizophrenia", "anxiety", "eating disorders", "personality disorders")
plots_by_diagnosis_py <- list()
for (i in main_diagnoses) {
  plot <- data_plotting_adm %>%
    filter(F_group == i) %>%
    filter(year %% 1 == 0) %>% 
    arrange(year) %>%
    ggplot(aes(x = year,
               color = age_group)) +
    geom_line(aes(y = f_diag_stand),
              size = 1.5) +
    geom_line(aes(y = pred),
              size = 3,
              color = "black") +
    labs(x = "time", 
         y = "per 100 000 person-years",
         title = i,
         subtitle = "Admission rates",
         color = "age category") +
    theme(plot.title = element_text(size = 8,
                                    face = "bold")) +
    facet_wrap(~ sex + age_cat,
               scales = "free",
               ncol = 2) +
    theme_bw()  
  plots_by_diagnosis_py[[i]] <- plot
}

file_names_py <- stringr::str_c(names(plots_by_diagnosis_py), "_adm", ".png")
pwalk(list(file_names_py, plots_by_diagnosis_py),
      ggsave,
      width=8, 
      height=6,
      path = "C:/Users/hp/Desktop/NUDZ/Time_Trends_Hospitalizations/Results/Graphs/Supplement")
