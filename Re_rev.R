## This code uses R0 to calculate R0 from the exponential growth phase of the epidemic

rm(list=ls())

setwd('~/USRe//')
require(EpiEstim)
require(tidyr)
require(readr)
require(dplyr)
require(ggplot2)
require(forcats)
require(lubridate)
require(googlesheets)
require(RCurl)
require(viridis)
require(flexdashboard)
#require(epuRate)
require(R0)
require(here)
require(rjson)
require(jsonlite)
require(RCurl)
require(highcharter)
require(here)
require(incidence)
require(purrr)
require(magrittr)

theme_set(theme_classic(base_size = 16))

require(RCurl)
data <- read.csv(text=getURL("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"),skip=1) #read from github

names(data) <- c('time','county','state','fips','cases','deaths')

data %>%
  group_by(state,county) %>%
  summarize() -> test


data %>%
  #subset(state == 'Missouri') %>%
  mutate(county = as.character(county))  -> data

data$time <- as.Date(as.character(data$time),format='%Y-%m-%d')


states <- unique(data$state)
ct <- 0

data$area <- paste0(data$state,' ',data$county, ' ','County')

sti = 0 
ct = 0
Re.dat <- NULL
for(st in states){
  
  sti <- sti + 1 
  
  print(sti / (length(unique(data$state))))
  
  data %>%
    mutate(Date = time) %>%
    mutate(Confirmed = deaths) %>%
    dplyr::select(time,county,state,Confirmed,deaths,Date) %>%
    subset(st == state) -> st.dat
  
  counties <- unique(st.dat$county)
  
  for(cty in counties){
  
    
    st.dat %>%
      mutate(Date = time) %>%
      mutate(Confirmed = deaths) %>%
      dplyr::select(time,county,state,Confirmed,deaths,Date) %>%
      subset(county == cty) -> dat
    
    dat %>%
      dplyr::select(Date,Confirmed) -> covid_pt
    
    covid_pt %>%
      mutate(
        Confirmed_lag = lag(x = Confirmed,
                            n = 1,
                            order_by = Date),
        Confirmed_var=Confirmed-Confirmed_lag,
        Confirmed_sign=if_else(Confirmed_var>=0,"+","-")
      ) %>%
      subset(Date >  first.date) -> covid_pt
    
    covid_pt  %>%
      dplyr::select(
        Date,Confirmed_var
      )  %>%
      dplyr::mutate(
        t_start = dplyr::row_number() %>% as.numeric(),
        t_end = t_start + 6
      ) -> covid_r
    
    covid_r$Confirmed_var[ covid_r$Confirmed_var < 0 ] <- 0
    
    
    data.cutoff <- 5
    cases.cutoff <- 10
    
    covid_r <- covid_r %>%
      subset(Date <= as.Date('2020-04-15'))
    
    if(nrow(covid_r) >= data.cutoff && max(covid_r$Confirmed_var, na.rm = T) >= cases.cutoff
       && cty != 'Colonial Heights city' && cty != 'Unknown' && cty != 'Tarrant'
       && cty != "Bristol"){
      
      #plot(covid_pt)
      gen_time = generation.time(type='gamma',val=c(4.8,2.3))
      begin_date = min(covid_r$Date[covid_r$Confirmed_var > 0], na.rm = T)
      if(begin_date + 14 > max(covid_r$Date)){
        end_date = as.Date(max(covid_r$Date), format = '%Y-%m-%d')}else{
        end_date = begin_date + 14
        }
      if(cty == 'Baltimore'){begin_date = as.Date("2020-04-04", format = '%Y-%m-%d')
      end_date = as.Date("2020-04-14", format = '%Y-%m-%d')}
      
      est_R <- estimate.R(covid_r$Confirmed_var,
                          t=covid_r$Date, method = 'EG', GT = gen_time,
                          begin = begin_date,
                          end = end_date)
      
      cty.Re <- data.frame(
        Date=covid_pt$Date,
        Confirmed = covid_pt$Confirmed,
        R_e_median = est_R$estimates$EG$R,
        R_e_q0025 = est_R$estimates$EG$conf.int[1],
        R_e_q0975 = est_R$estimates$EG$conf.int[2],
        growth_rate = est_R$estimates$EG$r,
        r_q0025 = est_R$estimates$EG$conf.int.r[1],
        r_q0975 = est_R$estimates$EG$conf.int.r[2],
        r_sq = est_R$estimates$EG$Rsquared,
        county = cty,
        state = st
      )
      
      
    }else{
      
      cty.Re <- data.frame(
        Date = covid_pt$Date,
        Confirmed = covid_pt$Confirmed,
        R_e_median = rep(NA,nrow(covid_pt)),
        R_e_q0025 = rep(NA,nrow(covid_pt)),
        R_e_q0975 = rep(NA,nrow(covid_pt)),
        growth_rate = rep(NA,nrow(covid_pt)),
        r_q0025 = rep(NA,nrow(covid_pt)),
        r_q0975 = rep(NA,nrow(covid_pt)),
        r_sq = rep(NA,nrow(covid_pt)),
        county = cty,
        state = st
      )
      
    }
    if(ct == 1){
      Re.dat <- cty.Re
    }else{
      Re.dat <- rbind(Re.dat,cty.Re)
    }
    
    save(Re.dat,file='county_rates_cases.RData')
    
  }
}

Re.dat$county_new <- paste0(Re.dat$county,' County')
save(Re.dat,file='county_rates_cases.RData')
Re.dat.nona <- Re.dat %>% filter(!is.na(R_e_median))
tmp.Re <- Re.dat
