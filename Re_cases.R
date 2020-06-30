## This code uses EpiEstim to calculate Re from case or mortality data
## you will need to select the data you prefer in the code below

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
  group_by(state,county) %>%
  summarize() -> test


data %>%
  #subset(state == 'Missouri') %>%
  mutate(county = as.character(county))  -> data

data$time <- as.Date(as.character(data$time),format='%Y-%m-%d')


states <- unique(data$state)
ct <- 0

data$area <- paste0(data$state,' ',data$county, ' ','County')

## if you want to calculate from death data, you need to change "cases" to "deaths"
## in the sections below 

for(st in states){
  
  data %>%
    mutate(Date = time) %>%
    mutate(Confirmed = cases) %>%
    select(time,county,state,Confirmed,cases,Date) %>%
    subset(st == state) -> st.dat
  
  counties <- unique(st.dat$county)
  
  for(cty in counties){
    
    if(st == "Iowa" & which(counties == cty ) < 61){next}
    ct <- ct + 1 
    
    print(ct / (length(unique(data$county))))
    
    st.dat %>%
      mutate(Date = time) %>%
      mutate(Confirmed = cases) %>%
      select(time,county,state,Confirmed,cases,Date) %>%
      subset(county == cty) -> dat
    
    dat %>%
      select(Date,Confirmed) -> covid_pt
    
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
        select(
          Date,Confirmed_var
        )  %>%
        dplyr::mutate(
          t_start = dplyr::row_number() %>% as.numeric(),
          t_end = t_start + 6
        ) -> covid_r
      
      covid_r$Confirmed_var[ covid_r$Confirmed_var < 0 ] <- 0
      
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
          covid_r$Confirmed_var,
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
          covid_r$Confirmed_var,
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
        select(Date,Confirmed,R_e_median,R_e_q0025,R_e_q0975,county,state) -> cty.Re
      
      covid_pt %>%
        subset(Date < min(cty.Re$Date)) %>%
        mutate(
          R_e_median = NA,
          R_e_q0025 = NA,
          R_e_q0975 = NA,
          county = cty,
          state = st
        ) %>%
        select(Date,Confirmed,R_e_median,R_e_q0025,R_e_q0975,county,state) -> missing.dat
      
      
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
save(Re.dat,file='county_Re.RData')

Re.dat %>%
  ggplot(aes(R_e_median)) + geom_histogram(bins=100)


test <- read.csv('~/MOmodeling/mo_county_contacts.csv')

Re.dat$county %>% table()

Re.dat %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(county = factor(county)) %>%
  subset(county == 'St. Louis') %>%
  subset(state == 'Missouri') %>%
  group_by(county) %>%
  ggplot() +
  geom_line(aes(Date,Confirmed,color=county),size=2)+
  facet_wrap(~county,scales='free')+
  ylab('confirmed cases')+
  xlab(NULL)+
  xlim(min(Re.dat$Date),max(Re.dat$Date)+1)+
  theme(legend.position = 'none')


test %>%
  mutate(Date = as.Date(date,format='%Y-%m-%d')) %>%
  mutate(county = factor(county)) %>%
  subset(county %in% nms)  %>%
  group_by(county) %>%
  ggplot() +
  geom_line(aes(Date,google_contact,color=county),size=2)+
  facet_wrap(~county)+
  ylab('change in contact')+
  geom_hline(yintercept = 1,linetype=2)+
  xlab(NULL)+
  scale_y_continuous(labels=scales::percent)+
  xlim(min(Re.dat$Date),max(Re.dat$Date)+1)+
  theme(legend.position = 'none') -> google_p



Re.dat %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(county = factor(county)) %>%
  subset(county %in% nms) %>%
  group_by(county) %>%
  ggplot() +
  geom_line(aes(Date,R_e_median,color=county),size=2)+
  geom_ribbon(aes(Date,ymin=R_e_q0025,ymax=R_e_q0975,fill=county),alpha=0.8)+
  facet_wrap(~county,scales='free')+
  ylab('effective transmission rate')+
  xlab(NULL)+
  geom_hline(yintercept = 1,linetype=2)+
  xlim(min(Re.dat$Date),max(Re.dat$Date)+1)+
  theme(legend.position = 'none') -> Re_p


require(cowplot)
jpeg('Re_ex.jpeg',height=800,width=1200)
plot_grid(confirmed_p,google_p,Re_p,nrow=3,align='v',axis='l')
dev.off()

Re.dat %>%
  mutate(county = factor(county)) %>%
  group_by(county) %>%
  top_n(1, Date) -> test

setdiff(counties,test$county)

### the model can account for imported cases at the beginning of the outbreak
### perhaps consider the first 2 ??? or 4 on the first two days
### if so, it will most likely decrease the initial R_e but increase it afterwards due to undiagnosed community transmission

covid_r_inc <-
  rep(
    x = unlist(covid_r$Date),
    times = unlist(covid_r$Confirmed_var)
  ) %>%
  incidence(
    dates = .,
    interval = "1 day",
    standard = TRUE,
    first_date = min(covid_r$Date),
    last_date = max(covid_r$Date)
  )

### find peak for adjustment of trend on model by the split argument
covid_r_inc_peak <- find_peak(covid_r_inc)

### fit log-linear model
### fits two exponential models to incidence data,
### of the form: log(y) = r * t + b , where
### 'y' is the incidence,
### 't' is time (in days)
### 'r' is the growth rate
### 'b' is the origin
### function fit will fit one model by default,
### but will fit two models on either side of a splitting date
### (typically the peak of the epidemic) if the argument split is provided
covid_r_inc_model <-
  fit(
    x = covid_r_inc,
    # split = covid_r_inc_peak,
    NULL
  )

# check object entirely
covid_r_inc_model

# (daily growth rate)
covid_r_inc_model$info$r
covid_r_inc_model$info$r.conf

# (doubling time in days)
covid_r_inc_model$info$doubling
covid_r_inc_model$info$doubling.conf

# incidence predictions (fitted vs observed data)
plot(covid_r_inc, fit = covid_r_inc_model,ylab='daily incidence')

### predict number cases next 3 days maintaing current exponential growth
### model elements for forecast are in covid_r_inc_model$model
### structure of dataset for prediction can be checked with
# head(covid_r_inc_model$info$pred)
### must provide x-axis data as a mid-point from t_0
### create x vector for forecasting on the next 3 days (reasonable amount time)
case_pred_3_day <-
  data.frame(
    dates = covid_r_inc_model$info$pred$dates[nrow(covid_r_inc_model$info$pred)] + 1:7,
    dates.x = covid_r_inc_model$info$pred$dates.x[nrow(covid_r_inc_model$info$pred)] + 1:7
  )

n_case_pred_3_day <-
  predict(
    object = covid_r_inc_model$model,
    newdata = case_pred_3_day,
    se.fit = TRUE,
    # type = "response",
    interval = "prediction"
  )

### log-linear model
### predictions are in log scale
### anti-log to get final count predictions
n_case_pred_3_day <-
  exp(x = n_case_pred_3_day[["fit"]])

case_pred_3_day <-
  dplyr::bind_cols(
    case_pred_3_day,
    as.data.frame(n_case_pred_3_day)
  ) %>%
  mutate(
    type = "predict"
  )

case_obs_fit <-
  covid_r_inc_model$info$pred %>%
  mutate(
    type = "fit"
  )

### final prediction
covid_pred_3_day <-
  bind_rows(
    case_obs_fit,
    case_pred_3_day
  )

### plot not perfect
### points and lines not coinciding on the x axis with the geom_col
plot_pred_3_day <-
  ggplot() +
  geom_point(
    data = covid_r,
    mapping = aes(x = Date, y = Confirmed_var),
    fill = "grey90"
  ) +
  geom_point(
    data = covid_pred_3_day,
    mapping = aes(x = dates, y = fit, colour = type),
    alpha = 0.7,
    size = 1.5
  ) +
  geom_line(
    data = covid_pred_3_day,
    mapping = aes(x = dates, y = fit, colour = type),
    alpha = 0.7,
    size = 1.5
  ) +
  geom_ribbon(
    data = covid_pred_3_day,
    mapping = aes(x = dates, ymin = lwr, ymax = upr, fill = type),
    alpha = 0.25
  ) +
  scale_x_date(breaks = "1 day") +
  scale_y_continuous(
    limits = c(0, max(covid_pred_3_day$upr)),
    breaks = pretty(covid_pred_3_day$upr),
    labels = pretty(covid_pred_3_day$upr)
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('cases') + xlab(NULL)+
  ggtitle('Effective reproductive number projection')

### display plot
plot_pred_3_day

covid_pred_3_day %>%
  subset(type == 'predict') %>%
  select(c(fit, lwr,upr)) %>%
  colSums()


### computes the overall infectivity (lambda) due to previously infected individuals
### λ_t = ∑_{k=1}^{t-1}I_{t-k}w_k
lambda_covid_pt <-
  overall_infectivity(
    incid = data.frame(I = covid_r_inc$counts),
    si_distr = discr_si(k = c(100, 1:100), mu = 4.7, sigma = 2.9)
  )

plot(lambda_covid_pt)

