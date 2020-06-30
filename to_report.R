
## This is the code to analyze Re and R0 estiamtes, you will need esitmates 
## before running this code. 

## Please see Re_rev.R to calculate R0 from growth rates. 
## Please see Re_lag.R to calculate Re from stochastic back-forecasted mortality data.
## Please Re_cases.R to calculate Re from cases or mortality data without back-forecasting. We do not report 
## these estimates in the main text of the paper, but we used them as a robustness check and they 
## were qualitatively consistent. 

## NOTE if you load data which used EpiEstim to estimate Re directly from death data 
## (from Re_cases.R) we recommend using lag_days = 21 in the 'data_aggregator' function
## in order match dates with appropriate Re values.
## For cases (from Re_cases.R), the lag will be closer to 7-10 days.
## If you load data which used back-forecasted mortality to reconstruct incidence,
## (as in Re_lag.R), then the data has already been lagged appropriately so you should
## set lag_days = 0. 



## write function in order to facilitate faster data aggregation with multiple
## estimation data frames
require(tidyr)
require(dplyr)
require(ggplot2)

acts_iso <- read.csv('acts_iso.csv')
county_popden <- read.csv('county_popden.csv')

metro_areas <- c('New York New York City County','Massachusetts Suffolk County',
                 'Illinois Cook County', 'Pennsylvania Philadelphia County',
                 'Colorado Denver County','California Orange County','Wisconsin Milwaukee County',
                 'Texas Dallas County','Virginia Fairfax County','Texas Harris County',
                 'California Los Angeles County','Ohio Franklin County','Texas Tarrant County',
                 'Minnesota Hennepin County','North Carolina Mecklenburg County','California Sacramento County',
                 'Utah Salt Lake County','Rhode Island Providence County','Georgia Fulton County',
                 'Missouri St. Louis County','Florida Miami-Dade County')


data_aggregator <- function(covid_data, lag_days = 21){
  
  covid_data$county_new <- paste(covid_data$county, 'County')
  covid_data$county_state <- paste(covid_data$state, covid_data$county_new, sep = ' ')
  covid_data$state_county <- paste(covid_data$state, covid_data$county, sep = ' ')
  covid_data$county_state[covid_data$state == 'Louisiana'] <- paste(covid_data$state_county[covid_data$state == 'Louisiana'], 'Parish')
  nlag <- lag_days
  covid_data$date <- covid_data$Date - nlag
  
  Re_acts <- inner_join(covid_data, acts_iso, by = c('date', 'state'))
  Re_acts[is.na(Re_acts$level) & Re_acts$Date < '2020-02-26','level'] <- 0
  
  Re_acts <- left_join(Re_acts, county_popden, by='county_state')

  Re_acts$value[is.na(Re_acts$value)] <- 0
  Re_acts$state[Re_acts$county == 'New York City'] <- 'New York City'
  
  Re_acts <- Re_acts %>% filter(!is.na(R_e_median))
  
  begin <- as.Date("2020-02-15", format = '%Y-%m-%d')
  Re_acts$cal_time <- Re_acts$date - begin
  
  Re_acts <- Re_acts %>%
    group_by(state_county) %>%
    mutate(days = date - min(date),
           days1 = sum(level==1, na.rm = T),
           days2 = sum(level==2, na.rm = T),
           days3 = sum(level==3, na.rm = T),
           Re_percent = ((R_e_median - mean(R_e_median[level == 0]))/mean(R_e_median[level == 0]))*100)
  
  Re_acts$level1 <- ifelse(Re_acts$level == 1 | Re_acts$level == 2 | Re_acts$level == 3, 1, 0)
  Re_acts$level2 <- ifelse(Re_acts$level == 2 | Re_acts$level == 3, 1, 0)
  Re_acts$level3 <- ifelse(Re_acts$level == 3, 1, 0)
  
  
  return(Re_acts)
}


early_summary <- function(covid_acts, metros = TRUE, ndays = 5){
  metro_areas <- c('New York New York City County','Massachusetts Suffolk County',
                   'Illinois Cook County', 'Pennsylvania Philadelphia County',
                   'Colorado Denver County','California Orange County','Wisconsin Milwaukee County',
                   'Texas Dallas County','Virginia Fairfax County','Texas Harris County',
                   'California Los Angeles County','Ohio Franklin County','Texas Tarrant County',
                   'Minnesota Hennepin County','North Carolina Mecklenburg County','California Sacramento County',
                   'Utah Salt Lake County','Rhode Island Providence County','Georgia Fulton County',
                   'Missouri St. Louis County','Florida Miami-Dade County')
  if(metros == TRUE){sub = covid_acts %>% filter(county_state %in% metro_areas)
  }else{sub = covid_acts}
  
  sub_sum <- sub %>% group_by(state, county_state) %>%
    summarize(begRe  = mean(R_e_median[1:ndays]),
              rate = growth_rate[1],
              rmax = mean(r_q0975[1:ndays]),
              rmin = mean(r_q0025[1:ndays]),
              #bmax = mean(R_e_median[1:ndays]) + 1.96*sd(R_e_median[1:ndays])/ndays,
              bmax = mean(R_e_q0975[1:ndays]),
              #bmin = mean(R_e_median[1:ndays]) - 1.96*sd(R_e_median[1:ndays])/ndays,
              bmin = mean(R_e_q0025[1:ndays]),
              meanRe = mean(R_e_median),
              noRe = mean(R_e_median[level ==0], na.rm = T),
              time_high = length(level == 3),
              Re_high = mean(R_e_median[level==3], na.rm = T),
              hmax =  mean(R_e_median[level==3], na.rm = T) + 1.96*(sd(R_e_median[level==3], na.rm = T))/time_high,
              hmin =  mean(R_e_median[level==3], na.rm = T) - 1.96*(sd(R_e_median[level==3], na.rm = T))/time_high,
              pd = popden[1],
              pd_rob = pop[1]/(land[1]/1e6),
              rural = rural[1],
              urban = urban[1],
              pop=pop[1],
              start = min(cal_time),
              risk_prop = risk_prop[1])
  return(sub_sum)
  }

## load R0 estimates from Re_rev.R
load('county_rates.RData')
deaths_rates <- data_aggregator(Re.dat)

deaths_sub <- early_summary(deaths_rates, metros = FALSE)

deaths_sub$pd_rob[deaths_sub$state == 'District of Columbia'] <- 5047

## create summary table for supplement 
summary_table <- data.frame(county = deaths_sub$county_state, R0 = deaths_sub$begRe, density = deaths_sub$pd_rob)
summary_table <- summary_table[order(summary_table$density, decreasing = F),]

## log density
deaths_sub$log_den <- log10(deaths_sub$pd_rob)
deaths_sub$logpop <- log10(deaths_sub$pop)

## linear fit for density 
fit <- lm(begRe ~ log_den, data = deaths_sub)
## update with segmented fit
fit_seg <- segmented(fit, seg.Z = ~log_den, psi = 3)

## linear fit for density, metro counties only 
fit_sub <- lm(begRe ~ log10(pd_rob), data = deaths_sub %>% filter(county_state %in% metro_areas))

## linear fit for pop
deaths_sub$pop[deaths_sub$state == 'District of Columbia'] <- 705749
fit_pop <- lm(begRe ~ log10(pop), data = deaths_sub)


## run this to make Figure One 

## density plot with predicted values and ribbon CI
pdat <- with(deaths_sub, data.frame(log_den = log10(seq(min(pd_rob),
                                       max(pd_rob), length = 100))))
tmp2 <- predict(fit_seg, newdata = pdat, se.fit = TRUE)
tmp2$pop_den <- 10^(pdat$log_den)

denplot <- ggplot() + geom_point(aes(x=pd_rob, y=begRe), data = deaths_sub) +
  geom_linerange(aes(x=pd_rob, ymin = bmin, ymax = bmax), alpha = .25, data = deaths_sub)+
  geom_line(aes(x=tmp2$pop_den, y=tmp2$fit)) +
  geom_ribbon(aes(x=tmp2$pop_den, ymax=tmp2$fit + 1.96*tmp2$se.fit, ymin = tmp2$fit - 1.96*tmp2$se.fit), alpha = .5) +
  scale_x_log10() + theme_classic() + ylab(expression(Inferred ~R[0])) + xlab('log population density') +
  ggtitle('Early Transmission and Population Density') + ylim(0,9.5)


## do the same for metro counties 
pdat_sub <- with(deaths_sub %>% filter(county_state %in% metro_areas), data.frame(pd_rob = (seq(min(pd_rob),
                                                        max(pd_rob), length = 100))))

tmp2_sub <- predict(fit_sub, newdata = pdat_sub, se.fit = TRUE)
tmp2_sub$pd_rob <- pdat_sub$pd_rob

inset <- ggplot() +
  geom_point(aes(x=pd_rob, y=begRe), data = deaths_sub %>% filter(county_state %in% metro_areas), col = 'tomato2') +
  scale_x_log10() + theme_classic() +
  ylab(NULL) + xlab(NULL) +
  geom_line(aes(x=tmp2_sub$pd_rob, y=tmp2_sub$fit), col = 'tomato2') +
  geom_ribbon(aes(x=tmp2_sub$pd_rob, ymax=tmp2_sub$fit + 1.96*tmp2_sub$se.fit,
                  ymin = tmp2_sub$fit - 1.96*tmp2_sub$se.fit), alpha = .5, fill = 'tomato2') +
  geom_linerange(aes(x=pd_rob, ymin = bmin, ymax = bmax), alpha = .25,
                 data = deaths_sub %>% filter(county_state %in% metro_areas), col = 'tomato2') +
  ggtitle('Primary Metro Counties')


## population plot with fit and CI ribbon
pdat_pop <- with(deaths_sub, data.frame(pop = (seq(min(pop),
                                                        max(pop), length = 100))))
tmp_pop <- predict(fit_pop, newdata = pdat_pop, se.fit = TRUE)
tmp_pop$pop <- pdat_pop$pop

popplot <- ggplot() + geom_point(aes(x=pop, y=begRe), data = deaths_sub) +
  geom_linerange(aes(x=pop, ymin = bmin, ymax = bmax), alpha = .25, data = deaths_sub)+
  geom_line(aes(x=tmp_pop$pop, y=tmp_pop$fit)) +
  geom_ribbon(aes(x=tmp_pop$pop, ymax=tmp_pop$fit + 1.96*tmp_pop$se.fit, ymin = tmp_pop$fit - 1.96*tmp_pop$se.fit), alpha = .5) +
  scale_x_log10() + theme_classic() + ylab(expression(Inferred ~R[0])) + xlab('log population size') +
  ggtitle('Early Transmission and Population Size') + ylim(0,9.5)

## population inset for metro counties 
pdat_sub_pop <- with(deaths_sub %>% filter(county_state %in% metro_areas), data.frame(pop = (seq(min(pop),
                                                                                                 max(pop), length = 100))))
fit_sub_pop <- fit_sub <- lm(begRe ~ log10(pop), data = deaths_sub %>% filter(county_state %in% metro_areas))
tmp_pop_sub <- predict(fit_sub_pop, newdata = pdat_sub_pop, se.fit = TRUE)
tmp_pop_sub$pop <- pdat_sub_pop$pop

inset_pop <- ggplot() +
  geom_point(aes(x=pop, y=begRe), data = deaths_sub %>% filter(county_state %in% metro_areas), col = 'tomato2') +
  scale_x_log10() + theme_classic() +
  ylab(NULL) + xlab(NULL) +
  geom_line(aes(x=tmp_pop_sub$pop, y=tmp_pop_sub$fit), col = 'tomato2') +
  geom_ribbon(aes(x=tmp_pop_sub$pop, ymax=tmp_pop_sub$fit + 1.96*tmp_pop_sub$se.fit, ymin = tmp_pop_sub$fit - 1.96*tmp_pop_sub$se.fit),
              fill = 'tomato2', alpha = .5) +
  geom_linerange(aes(x=pop, ymin = bmin, ymax = bmax), alpha = .25,
                 data = deaths_sub %>% filter(county_state %in% metro_areas), col = 'tomato2') +
  ggtitle('Primary Metro Counties')

## make Figure One 
den_inset <- ggdraw(denplot ) +
  draw_plot(inset, .5, .55, .45, .4) 
pop_inset <- ggdraw(popplot ) +
  draw_plot(inset_pop, .5, .55, .45, .4) 
plot_grid(den_inset, pop_inset, labels = c('A','B'))

#### continuous estimates of Re from Re_lag.R  
load('county_Re.RData')

Re_acts <- data_aggregator(Re.dat, lag_days = 0, post = TRUE)
Re_acts$state[Re_acts$county == 'New York City'] <- 'New York City'

## drop outliers (generally single day spikes resulting from instability in estimates)
Re_acts <- Re_acts %>% filter(R_e_median < 20)

## we will drop observations from post opening 
Re_acts$post_open <- ifelse(Re_acts$date > as.Date('2020-05-17', format = '%Y-%m-%d') & Re_acts$level != 3, 1, 0)

## run this to make Figure 2
## want county level means, will plot these with lines connecting counties
Re_grouped <- Re_acts %>%
  filter(post_open == 0) %>%
  group_by(state_county,level, state) %>%
  summarize(meanRe = mean(R_e_median, na.rm = T),
            pd = popden[1],
            pop = pop[1])


ggplot(Re_grouped) + geom_point(aes(x=level,meanRe,group=level,col=level)) +
  geom_line(aes(x=level,meanRe,group=state_county), alpha = .2) +
  geom_hline(yintercept = 3, lty = 3)+
  geom_hline(yintercept = 1, lty = 3)+
  facet_wrap(vars(state), ncol = 7, nrow = 7) +
  scale_y_log10() +
  scale_color_manual(values = c('grey','chartreuse3','dodgerblue','tomato2'),
                                      labels = c('none','low','med','high'), name = 'response') +
  ggtitle(expression(paste('Lagged', ~R[e], ' estimates by state response level',sep = ' '))) +
  xlab('response level') + ylab(expression(Inferred ~R[e])) + theme_classic() +
  scale_x_discrete(breaks = c(0,1,2,3), labels = c('none','low','med','high')) +
  theme(axis.text.x = element_text(angle = 45))

### examine outliers, find they are mostly rural and start late 
Re_delta <- Re_acts %>% 
  filter(post_open == 0) %>%
  group_by(state_county, county_state, state) %>%
  summarize(pd = popden[1],
            Re0 = if(length(R_e_median[level==0])>5){mean(R_e_median[level == 0][(length(R_e_median[level==0])-5):length(R_e_median[level==0])])}
            else{NA},
            Re1 = if(length(R_e_median[level==1])>5){mean(R_e_median[level == 1][(length(R_e_median[level==1])-5):length(R_e_median[level==1])])}
            else{NA},
            Re2 =if(length(R_e_median[level==2])>5){mean(R_e_median[level == 2][(length(R_e_median[level==2])-5):length(R_e_median[level==2])])}
            else{NA},
            Re3 = if(length(R_e_median[level == 3])>5){mean(R_e_median[level == 3][(length(R_e_median[level==3])-5):length(R_e_median[level==3])])}
            else{NA})

Re_outlier <- Re_grouped %>%
  group_by(state_county) %>%
  dplyr::summarize(diff = if(length(meanRe[level == 2] - meanRe[level == 3])>0){meanRe[level == 2] - meanRe[level == 3]}else{NA})
Re_outlier$out <- ifelse(Re_outlier$diff < -1, 1, 0)
Re_outlier$out[is.na(Re_outlier$out)] <- 0
outliers <- Re_outlier$state_county[Re_outlier$out == 1]

## look at the distribution of density in places with increase Re between level 2 and 3
Re_acts %>% filter(state_county %in% outliers) %>%
  plyr::summarize(dist = quantile(popden, probs = seq(from = .1, to = 1, by = .1)))

#### Now we look at change in Re by state and map these 
Re_grouped_wide <- Re_grouped %>% pivot_wider(names_from = 'level',values_from = 'meanRe')

Re_grouped_wide %>% ggplot() + geom_point(aes(x=pd,y=`2`-`3`))

Re_state_delta <- Re_grouped_wide %>%
  group_by(state) %>% 
  summarize(delta_one = mean(`1`-`0`, na.rm = T),
            delta_two =  mean(`2`-`1`, na.rm = T),
            delta_three = mean(`3`-`2`, na.rm = T))

States <- map_data('state')
Re_state_delta$region <- tolower(Re_state_delta$state)
States_data <- left_join(States, Re_state_delta, by = 'region')

map1 <- ggplot() + 
  geom_polygon( data=States_data, aes(x=long, y=lat, group=group, fill = delta_one),
                color="white" ) + theme_void() + 
  scale_fill_gradient2(high= 'tomato2',low= 'chartreuse3', limits = c(-4.1,2), na.value = 'grey70',
                       name = expression(paste('change in',~R[e]))) +
  ggtitle('No NPI to Low') +
  theme(legend.position=c(0.89,0.25),
        legend.text = element_text(size = 20),title = element_text(size = 23))

map2 <- ggplot() + 
  geom_polygon( data=States_data, aes(x=long, y=lat, group=group, fill = delta_two),
                color="white" ) + 
  scale_fill_gradient2(high= 'tomato2',low= 'chartreuse3', limits = c(-4.1,2), na.value = 'grey70') + 
  theme_void() +guides(fill = FALSE) +
  theme(title = element_text(size = 23))+
  ggtitle('Low to Medium')


## we take places that have increased more that 50% from
## https://www.npr.org/sections/health-shots/2020/03/16/816707182/map-tracking-the-spread-of-the-coronavirus-in-the-u-s
## accessed june, 25

spikes <- c('oklahoma','florida','arizona','texas','idaho','kansas','oregon','georgia','tennessee',
            'washington','arkansas','california','ohio','alabama')

States_data$spike <- ifelse(States$region %in% tolower(spikes), 1, 0)

map3 <- ggplot() + 
  geom_polygon( data=States_data, aes(x=long, y=lat, group=group, fill = delta_three,
                color=as.factor(spike))) + theme_void() + 
  scale_color_manual(values = c('white','tomato2'),  na.value = 'grey70',name = 'June cases + > 50%', labels = c('no','yes')) +
  scale_fill_gradient2(high= 'tomato2',low= 'chartreuse3', limits = c(-4.1,2), na.value = 'grey70',
                       name = expression(paste('change in',~R[e]))) +
  geom_polygon( data=States_data %>% filter(spike == 1), aes(x=long, y=lat, group=group, fill = delta_three),
                                      color='tomato2') +
  theme(legend.position=c(0.15,0.2),
        legend.text = element_text(size = 20), title = element_text(size = 23)) + ggtitle('Medium to High') + guides(fill = FALSE)


plot_grid(map1,map2,map3, nrow = 1)


#### We also look at county by county differences 
#### First create a data frame to store results from t tests
#### Then calculate t tests for counties with 5 obs in the two levels being compared

stat_sig <- data.frame(state_county = unique(Re_acts$state_county))
stat_sig$t.val1 <- stat_sig$t.val2 <-stat_sig$t.val3 <- NA
stat_sig$mean1 <- stat_sig$mean2 <- stat_sig$mean3 <- NA
stat_sig$sig1 <- stat_sig$sig2 <- stat_sig$sig3 <- NA
stat_sig$pd <- NA
Re_pre <- Re_acts %>% filter(post_open == 0)

for(county in stat_sig$state_county){
  sub = Re_pre[Re_pre$state_county == county, ]
  if(length(sub$R_e_median[sub$level == 0]) > 5 & length(sub$R_e_median[sub$level == 1]) > 5){
  t = t.test(sub$R_e_median[sub$level == 0], sub$R_e_median[sub$level == 1])
  stat_sig[stat_sig$state_count == county, 't.val1'] <- t$statistic
  stat_sig[stat_sig$state_county == county, 'sig1'] <- ifelse(sign(t$conf.int[1]) == sign(t$conf.int[2]),
                                                        1, 0)
  stat_sig[stat_sig$state_county == county, 'mean1'] <- t$estimate[2] - t$estimate[1]
  }
  if(length(sub$R_e_median[sub$level == 1]) > 5 & length(sub$R_e_median[sub$level == 2]) > 5){
    t = t.test(sub$R_e_median[sub$level == 1], sub$R_e_median[sub$level == 2])
    stat_sig[stat_sig$state_count == county, 't.val2'] <- t$statistic
    stat_sig[stat_sig$state_county == county, 'sig2'] <- ifelse(sign(t$conf.int[1]) == sign(t$conf.int[2]),
                                                               1, 0)
    stat_sig[stat_sig$state_county == county, 'mean2'] <- t$estimate[2] - t$estimate[1]
  }
  if( length(sub$R_e_median[sub$level == 2]) > 5 & length(sub$R_e_median[sub$level == 3]) > 5){
    t = t.test(sub$R_e_median[sub$level == 2], sub$R_e_median[sub$level == 3])
    stat_sig[stat_sig$state_count == county, 't.val3'] <- t$statistic
    stat_sig[stat_sig$state_county == county, 'sig3'] <- ifelse(sign(t$conf.int[1]) == sign(t$conf.int[2]),
                                                               1, 0)
    stat_sig[stat_sig$state_county == county, 'mean3'] <- t$estimate[2] - t$estimate[1]
  }
  stat_sig[stat_sig$state_county == county, 'pd'] <- sub$popden[1]
}


## Now we look at how change in Re is associated with starting values

Re_state_delta <- Re_grouped_wide %>%
  group_by(state_county, state) %>% 
  summarize(delta_one = mean(`1`-`0`, na.rm = T),
            r0 = `0`,
            r1 = `1`,
            r2 = `2`,
            r3 = `3`,
            delta_two =  mean(`2`-`1`, na.rm = T),
            delta_three = mean(`3`-`2`, na.rm = T),
            pd = pd[1])

d1 <- (lm(delta_one ~ r0 , data = Re_state_delta))
d2 <- (lm(delta_two ~ r1, data = Re_state_delta))
d3 <- (lm(delta_three ~ r2, data = Re_state_delta))


## robustness checks for regression
metro_areas <- c('New York New York City','Massachusetts Suffolk',
                 'Illinois Cook', 'Pennsylvania Philadelphia',
                 'Colorado Denver','California Orange','Wisconsin Milwaukee',
                 'Texas Dallas','Virginia Fairfax','Texas Harris',
                 'California Los Angeles','Ohio Franklin','Texas Tarrant',
                 'Minnesota Hennepin','North Carolina Mecklenburg','California Sacramento',
                 'Utah Salt Lake','Rhode Island Providence','Georgia Fulton',
                 'Missouri St. Louis','Florida Miami-Dade')

d1_check <- (lm(delta_one ~ log10(pd), Re_state_delta %>% filter(state_county %in% metro_areas)))

Re_state_delta %>% filter(state_county %in% metro_areas) %>%
  ggplot() + geom_point(aes(x=log10(pd),y=delta_one))
Re_state_delta %>% filter(state_county %in% metro_areas) %>%
  ggplot() + geom_point(aes(x=log10(pd),y=r3))

d1_check <- lm(r1 ~ r0, Re_state_delta)
confint(d1_check)
d2_check <- lm(r2 ~ r1, Re_state_delta)
confint(d2_check)
d3_check <- lm(r3 ~ r2, Re_state_delta)
confint(d3_check)

## You can also examine these, but the r-squared is very low
f1 <- (lmer(R_e_median ~ level + (1|state), data = Re_acts %>% filter(level %in% c(0,1))))
f2 <- (lmer(R_e_median ~ level + (1|state), data = Re_acts %>% filter(level %in% c(1,2))))
f3 <- (lmer(R_e_median ~ level + (1|state), data = Re_acts %>% filter(level %in% c(2,3))))

## make up some data where Re2<Re1, and Re2a is smaller than Re1- but crucially also smaller for big ones
Re1 <- rnorm(1000,3,2)
Re2 <- Re1-rnorm(length(Re1),1,1)
Re2a <- Re1-Re1*0.5+rnorm(length(Re1),0,0.5)

plot(Re1,Re2)
points(Re1,Re2a,col=4)
abline(0,1)  ## ~ all points are below the 0,1

## for the first slope is not different from 1.
fit <- lm(Re2~Re1)
summary(fit)
abline(fit, col=2)

## for the second it is
fit1 <- lm(Re2a~Re1)
summary(fit1)
abline(fit1, col=2)

## and confidence interval is below 1
fit1$coeff[2]+c(-1,1)*1.96*summary(fit)$coeff[2,2]
