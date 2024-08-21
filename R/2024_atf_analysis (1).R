# load ----
remotes::install_github("afsc-assessments/afscdata")
remotes::install_github("BenWilliams-NOAA/afscassess")
remotes::install_github("BenWilliams-NOAA/safe")
library(afscdata)
library(afscassess)
library(dplyr)

# globals ---- input the values needed for this year's projection model, also make sure appropriate files are in the "base" folder, needs proj.dat and std.dat from most recent full accepted model to run following code
year <- 2023
styr <- 1977
full_model_yr <- 2021
folder <- "update"
species <- "ARTH"
TAC <- c(96969, 97372, 96501) # last three years of tac year-3, year-2, year-1
rec_age = 1 # recruitmemt age

# setup ---- creates folders for output of running projection model
# setup_folders(year)
# accepted_model(2021, "model_name", 2023) # FLAG use when have previous full model structure setup

# query data ---- pulls data from akfin and cleans/prepares catch with expanded catch for final year(s)
goa_atf(year, off_yr = TRUE)
clean_catch(year, species, TAC = TAC, fixed_catch = "goa_atf_catch_1977_1990.csv")
bts_biomass(year)

# run proj ---- runs projection model and puts results in folder (defined above and below)
proj_ak(year=2023, last_full_assess = 2021, species="atf", region="goa", off_yr=TRUE, folder = "update", rec_age = rec_age)

# FLAG for now the exec summary info in "base/processed" is not transferring to the exec summary in "update/processed" so just hand type in the old data for last yr in the new exec_summ in "update/processed"

# setup for next year
setup_folders(year+1)
# accepted_model(2021, "model_name", 2024) # FLAG use to transfer this proj output to next year

# create figures ----
dir.create(here::here(year, "figs"))

# plot catch/biomass
std <- read.delim(here::here(year, "base", "goa_atf.std"), sep="", header = TRUE)
catch <- vroom::vroom(here::here(year, "data", "output", "fsh_catch.csv"))

#make sure catch length matches model output length
catch <- filter(catch, year>=styr)

#note that may need to change the first catch filter if on different update schedule
std %>%
  filter(name=="totalbiomass") %>%
  bind_cols(filter(catch, year <= full_model_yr)) %>%
  filter(year >= 1991) %>%
  dplyr::select(year, catch, value, std.dev) %>%
  bind_rows(filter(catch, year > full_model_yr) %>%
              left_join(vroom::vroom(here::here(year, folder, "proj", "author_f", "bigsum.csv")) %>%
                          filter(Year > full_model_yr, Alt == 2) %>%
                          mutate(value = Total_Biom * 1000) %>%
                          dplyr::select(year=Year, value))) %>%
  mutate(std.dev = ifelse(is.na(std.dev), std.dev[year==full_model_yr], std.dev)) %>%
  mutate(lci = value - std.dev * 1.96,
         uci = value + std.dev * 1.96) %>%
  mutate(perc = catch / value,
         lci = catch / lci,
         uci = catch / uci,
         mean = mean(perc)) %>%
  dplyr::select(year, value, mean, perc, lci, uci) -> df

png(filename=here::here(year, "figs", "catch_bio.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
df %>%
  ggplot2::ggplot(ggplot2::aes(year, perc)) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2) +
  ggplot2::geom_hline(yintercept = df$mean, lty = 3) +
  ggplot2::expand_limits(y = c(0, 0.05)) +
  afscassess::scale_x_tickr(data=df, var=year, start = 1990) +
  afscassess::theme_report() +
  ggplot2::xlab("\nYear") +
  ggplot2::ylab("Catch/Biomass\n")


dev.off()

# plot survey results ----

vast <- vroom::vroom(here::here(year, "data", "user_input", "vast.csv")) # make this file from GAP output files in species Results folder, usually named Index
sb <- vroom::vroom(here::here(year, "data", "output",  "goa_total_bts_biomass.csv"))
sb <- filter(sb, year>1990)

png(filename=here::here(year, "figs", "bts_biomass.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)

vast  %>%
  dplyr::mutate(t = biomass,
                Model = "VAST")  %>%
  dplyr::select(-biomass) %>%
  dplyr::bind_rows(sb %>%
                     dplyr::rename(t = biomass) %>%
                     dplyr::mutate(Model = "Design-based") %>%
                     dplyr::select(-se)) %>%
  dplyr::group_by(Model) %>%
  dplyr::mutate(mean = mean(t)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(ggplot2::aes(year, t, fill = Model, color = Model)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  afscassess::scale_x_tickr(data=vast, var=year) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::ylab("Survey biomass (t)\n") +
  ggplot2::xlab("\nYear") +
  ggplot2::expand_limits(y = 0) +
  scico::scale_fill_scico_d(palette = "roma", begin = 0.2) +
  scico::scale_color_scico_d(palette = "roma", begin = 0.2) +
  ggplot2::geom_line(ggplot2::aes(y = mean), lty = 3) +
  afscassess::theme_report() +
  ggplot2::theme(legend.position = c(0.2, 0.8))

dev.off()
