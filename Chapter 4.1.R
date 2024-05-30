################################################################################
#Section 4.1
################################################################################

library(haven)
library(xts)
library(zoo)
library(haven)
library(dplyr)
library(lubridate)
library(DFA)
library(ggplot2)
#install.packages("DFA")

#################################################################### Declustered

data <- read_dta("../dst_final.dta") %>%
  mutate(t = dst_time)


summary(data)
data <- data %>% filter(!is.na(t))

data <- data %>%
  mutate(event_time = if_else(Dst <= -130, t, as.POSIXct(NA, origin = "1970-01-01")))

data <- data %>%
  mutate(
    episode = cumsum(if_else(lag(Dst, default = first(Dst)) > -130 & Dst <= -130, 1, 0)),
    is_new_episode = if_else(lag(Dst, default = Inf) > -130 & Dst <= -130 & row_number() > 1, 1, 0)
  )

data <- data %>%
  group_by(episode) %>%
  mutate(
    episode_start = dplyr::first(event_time,na_rm=TRUE),
    episode_end = dplyr::last(event_time,na_rm=TRUE)
  ) %>%
  ungroup()

episode_details <- data %>%
  group_by(episode) %>%
  summarise(
    start = first(episode_start),
    end = last(episode_end),
    min_Dst = min(Dst, na.rm = TRUE),
    .groups = 'drop'  # Ensure grouping is dropped after summarisation
  ) %>%
  mutate(
    gap_to_next = as.numeric(difftime(lead(start), end, units = "hours")),
    valid_episode = ifelse(is.na(gap_to_next) | gap_to_next >= 48, TRUE, FALSE)
  ) %>%
  filter(episode > 0)

#This is the runs declustering approach, also explained here https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007SW000329 
for (i in 1:(nrow(episode_details) - 1)) {

  if (!episode_details$valid_episode[i] & i < nrow(episode_details)) {
   
    current_min_Dst <- episode_details$min_Dst[i]
    current_start <- episode_details$start[i] 
    
    next_index <- i + 1
    
    while (next_index <= nrow(episode_details) & !episode_details$valid_episode[next_index]) {
      #If the current min_Dst is lower, update the next episode's min_Dst and start
      if (current_min_Dst < episode_details$min_Dst[next_index]) {
        episode_details$min_Dst[next_index] <- current_min_Dst
        episode_details$start[next_index] <- current_start
      }
      
      next_index <- next_index + 1
    }
    
    if (next_index <= nrow(episode_details)) {
      if (current_min_Dst < episode_details$min_Dst[next_index]) {
        episode_details$min_Dst[next_index] <- current_min_Dst
        episode_details$start[next_index] <- current_start  
      }
    }
  }
}

episode_details <- episode_details %>%
  rename(t = start) %>%
  mutate(date = as.Date(t),
         year = year(t),
         week = week(t)) %>%
  dplyr::select(episode, t, date, year, week, Dst = min_Dst, valid_episode) %>%
  filter (valid_episode == TRUE)

episode_details <- episode_details %>%
  mutate(
    cycle = sin((as.numeric(date - as.Date("2008-01-01")) / (11 * 365.25)) * 2 * pi),
    negDst = -Dst
  ) %>%
  arrange(negDst) %>%
  mutate(
    n = row_number(),
    compcdf = 1 - (n - 1) / n(),
    extremep = if_else(negDst >= 130, max(compcdf), NA_real_)
  )

#The solar cycles' dates are historical and widely available on the web. only used for the GAM approach (Section 5.5)
solar_cycles <- data.frame(
  cycle = 1:7,
  start = as.Date(c("1954-04-01", "1964-10-01", "1976-03-01", "1986-09-01", "1996-08-01", "2008-12-01", "2019-12-01")),
  end = as.Date(c("1964-10-01", "1976-03-01", "1986-09-01", "1996-08-01", "2008-12-01", "2019-12-01", "2030-12-01"))
)

episode_details <- episode_details %>%
  rowwise() %>%
  mutate(solar_cycle = {
    cycle <- NA
    for (i in 1:nrow(solar_cycles)) {
      if (date >= solar_cycles$start[i] && date < solar_cycles$end[i]) {
        cycle <- solar_cycles$cycle[i]
        break
      }
    }
    cycle
  }) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(cyc = {
    start_date <- solar_cycles$start[solar_cycles$cycle == solar_cycle]
    as.numeric(difftime(t, start_date, units = "hours"))
  }) %>%
  ungroup()

#episode_details$negDst is the main dataset for all our analysis that did not include the lower and higher estimates

############################################### Unclustered for later EI and GAM

#Remove the noise, -20 was given in multiple articles, see https://www.stce.be/news/616/welcome.html
high_intensity_events <- data %>%
  filter(Dst < -20) %>%
  mutate(negDst=abs(Dst))


#We can decluster it here but, given our approach, we selected episode_details as main declustered dataset
#library(climex)
#options(xts_check_TZ = FALSE)
#high_intensity_events_xts <- xts(high_intensity_events$negDst, order.by = high_intensity_events$t)
#high_intensity_events_xts <- climex::decluster(high_intensity_events_xts, threshold=131, cluster.distance = NULL, silent = FALSE, mc.cores = NULL)
#high_intensity_events_df <- data.frame(
 # t = index(high_intensity_events_xts),
  #negDst = coredata(high_intensity_events_xts)
#)
#high_intensity_events <- high_intensity_events_df
#sum(!is.na(high_intensity_events$negDst))


high_intensity_events <- high_intensity_events %>%
  mutate(
    date = as.Date(t),  
    year = lubridate::year(t),     
    week = lubridate::week(t)      
  ) %>%
  dplyr::select(t, date, year, week, negDst)  


solar_cycles <- data.frame(
  cycle = 1:7,
  start = as.Date(c("1954-04-01", "1964-10-01", "1976-03-01", "1986-09-01", "1996-08-01", "2008-12-01", "2019-12-01")),
  end = as.Date(c("1964-10-01", "1976-03-01", "1986-09-01", "1996-08-01", "2008-12-01", "2019-12-01", "2030-12-01"))
)

high_intensity_events <- high_intensity_events %>%
  rowwise() %>%
  mutate(solar_cycle = {
    cycle <- NA
    for (i in 1:nrow(solar_cycles)) {
      if (date >= solar_cycles$start[i] && date < solar_cycles$end[i]) {
        cycle <- solar_cycles$cycle[i]
        break
      }
    }
    cycle
  }) %>%
  ungroup()

high_intensity_events <- high_intensity_events %>%
  rowwise() %>%
  mutate(cyc = {
    start_date <- solar_cycles$start[solar_cycles$cycle == solar_cycle]
    as.numeric(difftime(t, start_date, units = "hours"))
  }) %>%
  ungroup()

############################################################### Diagnostic Plots
library(evmix)
evmix::mrlplot(high_intensity_events$negDst, tlim=c(0, 500), try.thresh = c(130,150, 200), p.or.n = FALSE)
evmix::tshapeplot(high_intensity_events$negDst, legend.loc = "topleft", tlim=c(0, 500), try.thresh = c(130,150, 200))
evmix::tscaleplot(high_intensity_events$negDst, tlim=c(0, 500), try.thresh = c(130,150, 200), legend.loc = "bottomleft")

library(ReIns)
H <- Hill(high_intensity_events$negDst, plot=FALSE)

M <- Moment(high_intensity_events$negDst)

gH <- genHill(high_intensity_events$negDst, gamma=H$gamma)
plot(H$k[1:5000], M$gamma[1:5000], xlab="k", ylab=expression(gamma), type="l", ylim=c(-0.2,0.5), main="Generalised Hill Estimator for the EVI")
lines(H$k[1:5000], gH$gamma[1:5000], lty=2, col="blue")
legend("bottomright", c("Moment", "Generalised Hill"), lty=1:2)

ReIns::Hill.2oQV(episode_details$negDst, plot=TRUE, col="blue")
mtext("Shape parameter", side=2, line=2)  # y-axis label
legend("bottomright", c("Hill 2oQV"), col=c("blue"), lty=1:1)

threshold <- 130


#MAX-TO-SUM ratio (following Taleb (2016))
library(evir)
summary(high_intensity_events[c("t", "negDst")])
MSplot <- function(data, p = 4) {
  data <- abs(data)
  plot_data <- data.frame(x = integer(), R = numeric(), p = integer())
  
  for (i in 1:p) {
    y <- data^i
    S <- cumsum(y)
    M <- cummax(y)
    R <- M / S
    plot_data <- rbind(plot_data, data.frame(x = 1:length(data), R = R, p = i))
  }
  
  ggplot(plot_data, aes(x = x, y = R, color = factor(p))) +
    geom_line(size = 1) +
    scale_color_brewer(palette = "Set1", name = "Moments") +
    labs(title = "MSplot for various moments", x = "n", y = "Rn") +
    ylim(0, 0.25) +
    theme_minimal()
}

MSplot(high_intensity_events$negDst)

# QQ plot and Mean Excesses plot
evir::qplot(high_intensity_events$negDst, xi=0)
evir::emplot(high_intensity_events$negDst, 'xy')
evir::meplot(high_intensity_events$negDst)

################################################################################
#Section 4.2
################################################################################

################################################################# Extremal Index
library(evd)
#The only R package offering to plot the output. Results were compared to other packages (extRemes, etc)
#You can switch the variable with episode_details$negDst

theta <- evd::exi(high_intensity_events$negDst, u=130, r = 0)
tlim <- quantile(high_intensity_events$negDst, probs = c(0.01,0.99), na.rm=T)
evd::exiplot(high_intensity_events$negDst, r=0, tlim=c(131,400))

par(mfrow = c(2, 2))  
acf(na.omit(episode_details$negDst), main = "ACF for the (runs) declustered dataset")
pacf(episode_details$negDst, main = "PACF for the (runs) declustered dataset")
acf(na.omit(high_intensity_events$negDst), main = "ACF for the full clustered dataset")
pacf(high_intensity_events$negDst, main = "PACF for the full clustered dataset")
par(mfrow = c(1, 1))

################################################################################
#Additional preprocessing steps needed
################################################################################

########################################################################## TALEB
library(dplyr)
str(episode_details)
L <- 130
H <- 32000

episode_details <- episode_details %>%
  mutate(Dual = L - H * log((H - negDst) / (H - L)))

#################################################################### Missing Obs

new_row <- data.frame(
  episode = NA,
  t = NA,
  date = NA,
  year = NA,
  week = NA,
  Dst = -1600,
  valid_episode = NA,
  cycle = NA,
  negDst = 1600,
  n = NA,
  compcdf = NA,
  extremep = NA,
  solar_cycle = NA,
  cyc = NA,
  Dual=1600
)
episode_details_upper <- bind_rows(episode_details, new_row)

new_row <- data.frame(
  episode = NA,
  t = NA,
  date = NA,
  year = NA,
  week = NA,
  Dst = -1182,
  valid_episode = NA,
  cycle = NA,
  negDst = 1182,
  n = NA,
  compcdf = NA,
  extremep = NA,
  solar_cycle = NA,
  cyc = NA,
  Dual=1182
)
episode_details_upper <- bind_rows(episode_details_upper, new_row)
#episode_details <- episode_details[-nrow(episode_details), ]

#We apply the Taleb approach for Truncation, see Section 2.2.4
H <- 32000
L <- 130
episode_details_upper <- episode_details_upper %>%
  mutate(Dual = L - H * log((H - negDst) / (H - L)))

new_row <- data.frame(
  episode = NA,
  t = NA,
  date = NA,
  year = NA,
  week = NA,
  Dst = -850,
  valid_episode = NA,
  cycle = NA,
  negDst = 850,
  n = NA,
  compcdf = NA,
  extremep = NA,
  solar_cycle = NA,
  cyc = NA,
  Dual=850
)
episode_details_lower <- bind_rows(episode_details, new_row)

new_row <- data.frame(
  episode = NA,
  t = NA,
  date = NA,
  year = NA,
  week = NA,
  Dst = -500,
  valid_episode = NA,
  cycle = NA,
  negDst = 500,
  n = NA,
  compcdf = NA,
  extremep = NA,
  solar_cycle = NA,
  cyc = NA,
  Dual=500
)
episode_details_lower <- bind_rows(episode_details_lower, new_row)
#episode_details <- episode_details[-nrow(episode_details), ]
episode_details_lower <- episode_details_lower %>%
  mutate(Dual = L - H * log((H - negDst) / (H - L)))


summary(episode_details_lower$Dual)
summary(episode_details_upper$Dual)
summary(episode_details$Dual)
