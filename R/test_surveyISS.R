# load packages ----
# if you don't have the afscdata package installed, you will need to install this first:
# devtools::install_github("afsc-assessments/afscdata", force = TRUE)
# now install surveyISS:
# devtools::install_github("BenWilliams-NOAA/surveyISS", force = TRUE)
library(surveyISS)
library(tidyverse)
library(vroom)

# set iterations ----
# first, is this a full run?
full_run = FALSE
# set number of desired bootstrap iterations for full run
iters_full = 500
# set number of iterations for testing run time
iters_test = 5
# set number of iters for this run
if(isTRUE(full_run)){
  iters = iters_full
} else{
  iters = iters_test}

# load data ----
if (!dir.exists(here::here("data"))) {
  data_nebs <- surveyISS::query_data(survey = c(98, 143),
                                     region = 'nebs',
                                     species = c(21720, 21740),
                                     yrs = 1979)
  data_goa <- surveyISS::query_data(survey = 47,
                                    region = 'goa',
                                    species = c(21720, 21740),
                                    yrs = 1990)
} else{
  data_nebs <- list(lfreq = vroom::vroom(here::here('data', 'nebs', 'lfreq.csv')),
                    specimen = vroom::vroom(here::here('data', 'nebs', 'specimen.csv')),
                    cpue = vroom::vroom(here::here('data', 'nebs', 'cpue.csv')),
                    strata = vroom::vroom(here::here('data', 'nebs', 'strata.csv')))
  data_goa <- list(lfreq = vroom::vroom(here::here('data', 'goa', 'lfreq.csv')),
                   specimen = vroom::vroom(here::here('data', 'goa', 'specimen.csv')),
                   cpue = vroom::vroom(here::here('data', 'goa', 'cpue.csv')),
                   strata = vroom::vroom(here::here('data', 'goa', 'strata.csv')))
}

# start run time test ----
if(iters < iters_full){
  tictoc::tic()
}

# run ISS with full bootstrap ----
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_nebs$lfreq,
                    specimen_data = data_nebs$specimen,
                    cpue_data = data_nebs$cpue,
                    strata_data = data_nebs$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = TRUE,
                    boot_ages = TRUE,
                    al_var = TRUE,
                    al_var_ann = TRUE,
                    age_err = TRUE,
                    region = 'nebs',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'bootall')
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_goa$lfreq,
                    specimen_data = data_goa$specimen,
                    cpue_data = data_goa$cpue,
                    strata_data = data_goa$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = TRUE,
                    boot_ages = TRUE,
                    al_var = TRUE,
                    al_var_ann = TRUE,
                    age_err = TRUE,
                    region = 'goa',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'bootall')


# run ISS with tow-level bootstrap ----
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_nebs$lfreq,
                    specimen_data = data_nebs$specimen,
                    cpue_data = data_nebs$cpue,
                    strata_data = data_nebs$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = FALSE,
                    boot_ages = FALSE,
                    al_var = TRUE,
                    al_var_ann = TRUE,
                    age_err = TRUE,
                    region = 'nebs',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'boottow')
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_goa$lfreq,
                    specimen_data = data_goa$specimen,
                    cpue_data = data_goa$cpue,
                    strata_data = data_goa$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = FALSE,
                    boot_ages = FALSE,
                    al_var = TRUE,
                    al_var_ann = TRUE,
                    age_err = TRUE,
                    region = 'goa',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'boottow')

# run ISS with full bootstrap, no ageing error ----
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_nebs$lfreq,
                    specimen_data = data_nebs$specimen,
                    cpue_data = data_nebs$cpue,
                    strata_data = data_nebs$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = TRUE,
                    boot_ages = TRUE,
                    al_var = FALSE,
                    al_var_ann = FALSE,
                    age_err = FALSE,
                    region = 'nebs',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'bootall_noae')
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_goa$lfreq,
                    specimen_data = data_goa$specimen,
                    cpue_data = data_goa$cpue,
                    strata_data = data_goa$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = TRUE,
                    boot_ages = TRUE,
                    al_var = FALSE,
                    al_var_ann = FALSE,
                    age_err = FALSE,
                    region = 'goa',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'bootall_noae')

# run ISS with tow-level bootstrap, no ageing error ----
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_nebs$lfreq,
                    specimen_data = data_nebs$specimen,
                    cpue_data = data_nebs$cpue,
                    strata_data = data_nebs$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = FALSE,
                    boot_ages = FALSE,
                    al_var = FALSE,
                    al_var_ann = FALSE,
                    age_err = FALSE,
                    region = 'nebs',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'boottow_noae')
surveyISS::srvy_iss(iters = iters,
                    lfreq_data = data_goa$lfreq,
                    specimen_data = data_goa$specimen,
                    cpue_data = data_goa$cpue,
                    strata_data = data_goa$strata,
                    yrs = 1979,
                    boot_hauls = TRUE,
                    boot_lengths = FALSE,
                    boot_ages = FALSE,
                    al_var = FALSE,
                    al_var_ann = FALSE,
                    age_err = FALSE,
                    region = 'goa',
                    save_interm = FALSE,
                    save_stats = FALSE,
                    save = 'boottow_noae')

# stop run time test ----
if(iters < iters_full){
  end <- tictoc::toc(quiet = TRUE)
  runtime <- round((((as.numeric(strsplit(end$callback_msg, split = " ")[[1]][1]) / iters) * iters_full) / 60) / 60, digits = 1)
  cat("Full run of", crayon::green$bold(iters_full), "iterations will take", crayon::red$bold$underline$italic(runtime), "hours", "\u2693","\n")
} else{
  cat("All", crayon::green$bold$underline$italic('Done'), "\u2693","\n")
}

# plot results ----

# read results
res_age <- vroom::vroom(here::here('output', 'nebs', 'bootall_iss_ag.csv')) %>% 
  tidytable::select(-nss, -nhls) %>% 
  tidytable::mutate(region = 'N+EBS',
                    boot = 'All') %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'bootall_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'All')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'boottow_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'Tow only')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'boottow_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'Tow only')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'boottow_noae_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'Tow only, no ae')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'boottow_noae_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'Tow only, no ae')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'bootall_noae_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'All, no ae')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'bootall_noae_iss_ag.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'All, no ae'))

res_len <- vroom::vroom(here::here('output', 'nebs', 'bootall_iss_ln.csv')) %>% 
  tidytable::select(-nss, -nhls) %>% 
  tidytable::mutate(region = 'N+EBS',
                    boot = 'All') %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'bootall_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'All')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'boottow_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'Tow only')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'boottow_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'Tow only')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'boottow_noae_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'Tow only, no ae')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'boottow_noae_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'Tow only, no ae')) %>% 
  tidytable::bind_rows(vroom::vroom(here::here('output', 'nebs', 'bootall_noae_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'N+EBS',
                                           boot = 'All, no ae')) %>%
  tidytable::bind_rows(vroom::vroom(here::here('output', 'goa', 'bootall_noae_iss_ln.csv')) %>% 
                         tidytable::select(-nss, -nhls) %>% 
                         tidytable::mutate(region = 'GOA',
                                           boot = 'All, no ae'))

# plot age results
res_age %>% 
  tidytable::filter(sex == 4) %>% 
  ggplot(aes(x = year, y = iss, col = boot)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_grid(species_code~region, scales = 'free_y') +
  labs(x = 'Year', y = 'Age ISS') +
  theme_bw()

# plot length results
res_len %>% 
  tidytable::filter(sex == 4,
                    boot %in% c('All', 'Tow only')) %>% 
  ggplot(aes(x = year, y = iss, col = boot)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_grid(species_code~region, scales = 'free_y') +
  labs(x = 'Year', y = 'Length ISS') +
  theme_bw()

