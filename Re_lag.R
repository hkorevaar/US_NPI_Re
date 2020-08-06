## This code uses stochasti back-forecasting to reconstruct incidence from mortality data
## We calculate Re from this synthetic incidence data.

rm(list=ls())

## adapted from https://github.com/aperaltasantos/covid_pt

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
  filter(state == 'Michigan' & county == "Wayne") %>%
  mutate(ind = 1:n()) %>%
  ggplot() + geom_line(aes(x=ind,y=deaths)) 

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


## need to add reporting delay and backwards distribution to reconstruct incidence
## from symptom to death: mean 15, sd 6.9
## from infection to symptom mean 5.3 sd 3.2

delay_dist <- function(death_data){
  ## first back distribute by mean 15, sd 6.9
  date_vec <- rep(death_data$Date, death_data$Confirmed_var)
  len <- length(date_vec)
  
  ## back dist with delays for symptoms to death and infection to symptoms 
  dist <- rgamma(len, shape = 4.726, scale = 3.174) + rgamma(len, shape = 2.743, scale = 1.932)
  delta_date <- date_vec - dist
  
  t_start = min(delta_date)
  t_end = max(delta_date)
  
  date_table <- (table(delta_date))
  df <- data.frame(Date = names(date_table), Confirmed = as.vector(date_table))
  df %>%
    mutate(Date = as.Date(Date)) %>%
    tidyr::complete(Date = seq.Date(t_start, t_end, by='day')) -> df
  
  df$Confirmed[is.na(df$Confirmed)] <- 0
  
  
  ## inflate assuming CFR of .7%
  ## results are pretty consistent across a range of values (.7-1.3%)
  df$Inferred <- df$Confirmed*130
  
  ## smooth for incidence
  sp <- smooth.spline(df$Inferred, lambda = .0001)
  df$Inferred_Smooth <- round(sp$y)
  df$Inferred_Smooth[df$Inferred_Smooth < 0] <- 0
  
  out_df <- data.frame(Date = df$Date, Confirmed = df$Inferred_Smooth)
  
  return(out_df)
  }
ct=0
for(st in states){
  
  data %>%
    mutate(Date = time) %>%
    mutate(Confirmed = deaths) %>%
    dplyr::select(time,county,state,Confirmed,deaths,Date) %>%
    subset(st == state) -> st.dat
  
  counties <- unique(st.dat$county)
  
  for(cty in counties){
    
    ct <- ct + 1 
    
    print(ct / (length(unique(data$county))))
    
    st.dat %>%
      mutate(Date = time) %>%
      mutate(Confirmed = deaths) %>%
      dplyr::select(time,county,state,Confirmed,deaths,Date) %>%
      subset(county == cty) -> dat
    
    dat %>%
      dplyr::select(Date,Confirmed) -> covid_pt
    
    data.cutoff <- 5
    cases.cutoff <- 10
    
    if(nrow(covid_pt) >= data.cutoff && max(dat$Confirmed) >= cases.cutoff && cty != 'Colonial Heights city' && cty != 'Unknown'){
      
      #plot(covid_pt)
      
      covid_pt<-covid_pt  %>%
        #subset(Date >= '2020-03-05') %>%
        mutate(epiweek = epiweek(Date))
      
      first.date <- head(covid_pt$Date,1) 
      
      covid_pt %>%
        mutate(
          Confirmed_lag = lag(x = Confirmed,
                              n = 1,
                              order_by = Date),
          Confirmed_var=Confirmed-Confirmed_lag,
          Confirmed_sign=if_else(Confirmed_var>=0,"+","-")
        ) %>%
        subset(Date >  first.date) -> covid_pt
      
      
      # covid_pt  %>%
      #   group_by(epiweek) %>%
      #   dplyr::summarise(
      #     incidence=sum(Confirmed_var)
      #   ) -> covid_r
      
      covid_pt  %>%
        dplyr::select(
          Date,Confirmed_var
        )  %>%
        dplyr::mutate(
          t_start = dplyr::row_number() %>% as.numeric(),
          t_end = t_start + 6
        ) -> covid_r
      
      covid_r$Confirmed_var[ covid_r$Confirmed_var < 0 ] <- 0
      
      covid_r <- delay_dist(covid_r)
      
      covid_r  %>%
        dplyr::select(
          Date,Confirmed
        )  %>%
        dplyr::mutate(
          t_start = dplyr::row_number() %>% as.numeric(),
          t_end = t_start + 6
        ) -> covid_r
      
      
      
      ## Calculate Effective R (R_e or R_t)
      ## Authors A. Peralta-santos
      ## Based on https://cmmid.github.io/topics/covid19/current-patterns-transmission/global-time-varying-transmission.html
      
      ###Methods
      
      #Time-varying effective reproduction estimates were made with a 7-day sliding window using EpiEstim [4,5] adjusted for imported cases and assuming an uncertain serial interval with a mean of 4.7 days (95% CrI: 3.7, 6.0) and a standard deviation of 2.9 days (95% CrI: 1.9, 4.9) [6].
      #Time-varying estimates of the doubling time were made with a 7-day sliding window by iteratively fitting an exponential regression model.
      
      ### R_e calculation - Parametric SI method for
      ### Serial Interval - mean = 4.7 // sd = 2.9
      res_parametric_si <-
        estimate_R(
          covid_r$Confirmed,
          method ="parametric_si",
          config = make_config(
            list(
              mean_si = 4.7,
              std_si = 2.9
            )
          )
        )
      
      # plot(res_parametric_si, legend = FALSE)
      
      r_prt <- as.data.frame(res_parametric_si$R)
      
      # r_prt <- left_join(covid_r, r_prt, by="t_start")
      
      ### join by t-end
      r_prt <-
        left_join(
          x = covid_r,
          y = dplyr::select(
            r_prt,
            c("t_end", "Mean(R)", "Quantile.0.025(R)", "Quantile.0.975(R)")
          ),
          by = c("t_start" = "t_end")
        )
      
      r_prt <-
        r_prt %>%
        dplyr::rename(
          r_efect = "Mean(R)",
          r_low = "Quantile.0.025(R)",
          r_high = "Quantile.0.975(R)"
        )
      
      p1 <-
        ggplot()+
        geom_line(
          data = r_prt,
          aes(
            x = Date,
            y = r_efect
          ),
          alpha = 0.7,
          size = 1
        ) +
        geom_hline(
          yintercept=1,
          linetype="dashed",
          color = "black"
        ) +
        geom_hline(
          yintercept=0,
          color = "black"
        ) +
        geom_ribbon(
          data = r_prt,
          aes(
            ymin = r_low,
            ymax = r_high,
            x = Date
          ),
          alpha=0.5,
          fill = "grey70"
        ) +
        scale_x_date(
          breaks = "2 day",
          date_labels = "%b %d"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
        )+
        ggtitle("COVID-19 Effective reproduction")+
        ylab(bquote(R[e]))+xlab(NULL)
      
      # p1
      
      ### R_e calculation - Uncertainty method
      ### Serial Interval
      ### -- mean 4.7 (95% CrI: 3.7, 6.0)
      ### -- sd 2.9 (95% CrI: 1.9, 4.9)
      sens_configs <-
        make_config(
          list(
            mean_si = 4.7, std_mean_si = 0.7,
            min_mean_si = 3.7, max_mean_si = 6.0,
            std_si = 2.9, std_std_si = 0.5,
            min_std_si = 1.9, max_std_si = 4.9,
            n1 = 1000,
            n2 = 100,
            seed = 123456789
          )
        )
      
      Rt_nonparam_si <-
        estimate_R(
          covid_r$Confirmed,
          method = "uncertain_si",
          config = sens_configs
        )
      
      ### inspect R_e estimate
      #plot(Rt_nonparam_si, legend = FALSE)
      
      ## Posterio sample R_e estimate
      sample_windows <- seq(length(Rt_nonparam_si$R$t_start))
      #sample_windows <- Rt_nonparam_si$dates
      
      posterior_R_t <-
        map(
          .x = sample_windows,
          .f = function(x) {
            
            posterior_sample_obj <-
              sample_posterior_R(
                R = Rt_nonparam_si,
                n = 1000,
                window = x
              )
            
            posterior_sample_estim <-
              data.frame(
                window_index = x,
                window_t_start = Rt_nonparam_si$R$t_start[x],
                window_t_end = Rt_nonparam_si$R$t_end[x],
                date_point = covid_r[covid_r$t_start == Rt_nonparam_si$R$t_end[x], "Date"],
                Confirmed = covid_pt[covid_r$t_start == Rt_nonparam_si$R$t_end[x], "Confirmed"],
                R_e_median = median(posterior_sample_obj),
                R_e_q0025 = quantile(posterior_sample_obj, probs = 0.025,na.rm = T),
                R_e_q0975 = quantile(posterior_sample_obj, probs = 0.975,na.rm = T)
              )
            
            return(posterior_sample_estim)
            
          }
        ) %>%
        reduce(bind_rows)
      
      
      cty.Re <- posterior_R_t
      
      cty.Re  %>%
        mutate(county = cty) %>%
        mutate(state = st) %>%
        mutate(Date = date_point) %>%
       dplyr::select(Date,Confirmed,R_e_median,R_e_q0025,R_e_q0975,county,state) -> cty.Re
      
      covid_pt %>%
        subset(Date < min(cty.Re$Date)) %>%
        mutate(
          R_e_median = NA,
          R_e_q0025 = NA,
          R_e_q0975 = NA,
          county = cty,
          state = st
        ) %>%
        dplyr::select(Date,Confirmed,R_e_median,R_e_q0025,R_e_q0975,county,state) -> missing.dat
      
      
      cty.Re <- rbind(missing.dat,cty.Re)
      
      plot_posterior_R_t <-
        ggplot(data = posterior_R_t, mapping = aes(x = date_point, y = R_e_median)) +
        geom_line(alpha = 0.3, size = 1.2) +
        geom_ribbon(mapping = aes(ymin = R_e_q0025, ymax = R_e_q0975), alpha = 0.1) +
        geom_smooth(se = FALSE) +
        scale_x_date(
          date_breaks = "1 day",
          limits = c(min(covid_r$Date), max(posterior_R_t$date_point))
        ) +
        scale_y_continuous(
          breaks = 0:ceiling(max(posterior_R_t$R_e_q0975)),
          limits = c(0, NA)
        ) +
        geom_hline(yintercept = 1) +
        theme_classic()+ggtitle(cty)
      
      # plot_posterior_R_t %>% print()
      
      
    }else{
      
      cty.Re <- data.frame(
        Date = covid_pt$Date,
        Confirmed = covid_pt$Confirmed,
        R_e_median = rep(NA,nrow(covid_pt)),
        R_e_q0025 = rep(NA,nrow(covid_pt)),
        R_e_q0975 = rep(NA,nrow(covid_pt)),
        county = cty,
        state = st
      )
      
    }
    if(ct == 1){
      Re.dat <- cty.Re
    }else{
      Re.dat <- rbind(Re.dat,cty.Re)
    }
    
    save(Re.dat,file='county_Re.RData')
    
  }
}

Re.dat$county_new <- paste0(Re.dat$county,' County')
save(Re.dat,file='~/Desktop/county_Re.RData')

Re.dat %>%
  ggplot(aes(R_e_median)) + geom_histogram(bins=100)



