em.rain_caps_phase <- function(){
  
  # phase dates
  impactful_rain <- ymd(c("2007-01-01", "2010-02-01", "2015-12-01"))
  
  increase_dates <- impactful_rain + months(2)
  
  decrease_dates <- ymd(c("2008-03-01", "2011-09-01", "2016-11-01")) 
  
  low_dates <- c(ymd("2000-12-01"), decrease_dates[-3] + c(months(13), months(9)))
  
  # session info
  session <- read.csv("../data/em_sessiontable.csv") %>% 
    group_by(trip) %>%
    summarise(rain = mean(rain),
              captures = mean(phAbundance, na.rm = T),
              capturesSy = mean(syAbundance, na.rm = T)) %>%
    ungroup() %>% 
    mutate(year = lubridate::year(lubridate::ymd(trip)),
           trip = ymd(trip)) %>% 
    bind_rows(data.frame(trip = max(.$trip)+months(1:12)))
  
  ## phase data
  phase <- data.frame(low =  low_dates, 
                      decrease = ymd(decrease_dates), 
                      increase = increase_dates) %>% 
    pivot_longer(cols = c(low, decrease, increase)) %>% 
    filter(complete.cases(.)) %>% arrange(value) %>% 
    mutate(no = rep(1:3, each = 3),
           period = paste('period', 1:9),
           phaseNo = paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3)))
  
  session$phase <- NA
  session$period <- NA
  session$phaseNo <- NA
  for(i in 1:nrow(phase)){
    index <- which(session$trip >= phase$value[i])
    session$phase[index] <- phase$name[i]
    session$period[index] <- phase$period[i]
    session$phaseNo[index] <- phase$phaseNo[i]
    
  }
  
  return(session)
}

em.visualisation_subset <- function(rain_caps, ph){
  
  max_y <- max(rain_caps$captures, na.rm = T)
  
  subset_dates <- rain_caps %>% 
    filter(trip >= min(ymd(ph@other$ind.metrics$trip))) %>%
    mutate(phaseNo = factor(phaseNo, levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                     rep(1:3, each = 3)))) %>% 
    group_by(phase, phaseNo) %>% 
    summarise(trip = min(trip)) %>% 
    filter(complete.cases(phase)) %>% 
    arrange(trip) %>% 
    rename(subset = phaseNo)
  
  
  vis_dates <- c(subset_dates$trip, 
                 max(ymd(ph@other$ind.metrics$trip))+months(2))
  
  vis_comb <- combn(1:length(vis_dates),2) 
  vis_comb_seq <-t(vis_comb[,vis_comb %>% diff == 1])
  
  subset_data_visualise <- data.frame(subset = subset_dates$subset,
                                      trip1 = ymd(vis_dates[vis_comb_seq[,1]]),
                                      trip2 = ymd(vis_dates[vis_comb_seq[,2]])) %>% 
    mutate(trip3 = ymd(trip1) - days(1),
           trip4 = ymd(trip2) - days(1)) %>% 
    pivot_longer(cols = trip1:trip4, values_to = "trip") %>% 
    arrange(subset, trip) %>% 
    mutate(rain = rep(c(0, rep(max_y, 2), 0), nrow(subset_dates))) %>% 
    dplyr::select(-name) %>% 
    mutate(phase = sub("[[:digit:]]+", "", subset),
           phase = factor(ifelse(phase == 'L', 'low', ifelse(phase == 'I',
                                                             'increase', 
                                                             'decrease')),
                          levels = c('low','increase','decrease'))) %>% 
    relocate(phase)
  
  return(list(subset = subset_data_visualise,
              dates = subset_dates))
  
}


em.conceptual_phase <- function(startdate = "2007-01-01",
                                rep123 = c(0.07, 0.402,0.615)){
  
  ph <- readRDS('../data/pherm_genotypes.rds')
  rain_caps <- em.rain_caps_phase()
  rain_caps$phase <- factor(rain_caps$phase, 
                            levels = c('low','increase', 'decrease'))
  subphase <- em.visualisation_subset(rain_caps, ph)
  
  subset_data_visualise <- subphase$subset 
  start_dates <- subphase$dates 
  
  myfillsubtle <- c(low = "#FFF3CD", increase = "#D2E7FA",
                    decrease = "#F7D1DF")
  
  myfill <- rev(c(increase = '#1E88E5',
                  decrease = '#D81B60', 
                  low = '#FFC107'))
  
  myfill <- c(decrease = terrain.colors(10)[4],
              increase = terrain.colors(10)[1],
              low = terrain.colors(10)[7])
  
  
  
  grob_repeats <- grobTree(textGrob('Periods:',
                                    x=0.002,  
                                    y= 0.93,
                                    hjust=0, rot = 0,
                                    gp=gpar(col="black", 
                                            fontsize=11)))
  
  # grob_repeats <- grobTree(textGrob(paste("period"),
  #                                   x=rep123[1] +0.005,  
  #                                   y= 0.98,
  #                                   hjust=0, rot = 0,
  #                                   gp=gpar(col="black", 
  #                                           fontsize=9, 
  #                                           fontface="italic")))
  
  dd <- filter(rain_caps, 
               trip > ymd(startdate) & complete.cases(captures)) %>% 
    pivot_longer(cols = captures:capturesSy, names_to = 'species', values_to = 'captures') %>% 
    mutate(trip = ymd(trip),
           phaseNo = factor(phaseNo, levels = unique(rain_caps$phaseNo)))
  
  sdv <- subset_data_visualise %>% 
    filter(rain > 0) %>% 
    mutate(phaseNo = factor(subset, levels = unique(rain_caps$phaseNo)),
           period = phaseNo,
           species = 'captures') %>%
    group_by(phaseNo, phase, period, species) %>% 
    summarise(rain = mean(rain),
              trip2 = as.Date(mean(trip), origin = "1970-01-01") - days(125),
              trip3 = as.Date(mean(trip), origin = "1970-01-01")-days(220)) %>% 
    arrange(trip2) %>% 
    filter(!duplicated(phaseNo)) %>% 
    ungroup()  %>% 
    mutate(trip = ifelse(phase == 'increase', 
                         as.character(trip2), 
                         as.character(trip3)),
           trip = ymd(trip)+months(2)) 
  
  subset_data_visualise$species <- 'captures'
  subset_data_visualise2 <- subset_data_visualise
  subset_data_visualise2$species <- 'capturesSy'
  subset_data_visualise2$rain <- ifelse(subset_data_visualise2$rain > 0,
                                        10, 0)
  
  phaseVis <- rbind(subset_data_visualise,
                    subset_data_visualise2) %>% 
    filter(subset != 'L1')
  #  sdv$period <- (sdv$phaseNo)
  #sdv2 <- rbind(sdv[1,], sdv)
  # sdv2$period[1] <- 'Period: 1'
  # sdv2$trip[2] <- ymd('2005-11-15')
  # sdv2$trip[1] <- ymd('2003-12-01')
  # sdv2$species <- dd$species[1]
  # sdv2 <- sdv2[-2,]
  conceptual <- ggplot(filter(dd, species == 'captures'),
                       aes(ymd(trip), captures))+
    geom_area(data = filter(phaseVis, species == 'captures'), 
              aes(x = trip, y = rain*1.07, group = subset, fill = phase), 
              alpha = 0.40) +
    geom_area(alpha = 1, fill = "grey90", col = "grey80")+
    geom_bar(stat = "identity", alpha = 0.8, col="grey20", fill = "black")+
    scale_fill_manual(values = myfill, 
                      breaks = c('low', 'increase', 'decrease'),
                      name = 'population phase') +
    #facet_wrap(~species, nrow = 2, )+
    #facet_grid(rows = vars(species),  scale = 'free_y', space = 'free_y')+
     geom_text(data = sdv[-1,], aes(x = trip, y = rain + 5, 
                               label = period), hjust = -0.5, size = 4.5)+
    theme_classic()+
    xlab("Year") +
    ylab("Captures (per 100 trap nights)")+
    geom_vline(xintercept = start_dates$trip[-1],#[seq(1,9,3)], 
               lty = 2, col = "black", alpha = 0.4, lwd = 0.5)+
    # annotation_custom(grob_repeats)+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          # legend.position = c(0.11, 0.7),
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.margin = margin(0.01, 0.5, 0.01, 0.01, "cm")) +
    scale_x_date(date_breaks = "years" , date_labels = "%Y")
  conceptual
  return(conceptual)
}



# load -----------

#setwd('../phd-simulations/')
#source('./code/libraries.R')
# source('./code/functions_conceptual.R')
# source('./code/functions_genetic-data-creation.R')

# conceptual -------
con <- em.conceptual_phase('2007-01-01', c(0.24, 0.507,0.68)) 
#  theme(legend.position = 'none')

phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])

conditions  <- read.csv('../phd-simulations/output/conditions_fourth_sim.csv')

initialise_conditions$gen <- 0
initialise_conditions$year <- 0
initialise_conditions$pops <- 24
initialise_conditions$loci <- '24388'
initialise_conditions$fst <- 0.01684324
initialise_conditions$nInd <- 2400
initialise_conditions$phaseNo <- 'I1'

conditionsx <- bind_rows(conditions[1,],initialise_conditions[c(1:5,9:ncol(initialise_conditions))],
                         conditions) %>% 
  rename(Nmax = N)
conditionsx <- conditionsx[-1,] 
conditionsx$loci[2] <- '-'
# conditionsx$gen[1]<- 0
# conditionsx$year[1]<- 0
conditionsx$phase[1] <- 'initialise'

simCaps <- ggplot(conditionsx, aes(year, nInd/100, fill = phase))+
  geom_bar(stat = 'identity', alpha = 0.40) +
  theme_classic() +
  scale_x_continuous(breaks = conditionsx$year, labels = conditionsx$gen)+
  scale_fill_manual(values = c(phaseCol, initialise = 'grey60'),
                    name = 'Population Phase',
                    labels = c('Low', 'Increase', 'Decrease', 'Sim initialised'))+
  theme(legend.position = 'bottom',
        plot.margin = margin(0.01, 0.5, 0.01, 0.01, "cm"),
        axis.title.y = element_text(size = 12))+
  xlab('Generation')+
  ylab(expression('N'[' / 100']))

fig1 <- grid.arrange(con + labs(title = 'Sandy inland mouse captures')+ theme(legend.position = 'none',
                                 axis.title.y = element_text(size = 12)),
                     simCaps+labs(title = 'Simulated N'), ncol = 1)

ggsave('./figures/fig1_captures_sim.png',fig1,  units = 'cm',
       height = 16, width = 16, dpi = 300)


fx <- conditionsx[,-c(4,5, 6,7,8)] %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  flextable() %>% 
  autofit() %>% 
  align(j = 5, align = 'center') %>% 
  bold(part = 'header') %>% 
  bold(j = 7, i = 2:20) %>% 
  italic(j = 8) %>% 
  bg(j = 1:5, i = which(conditionsx$phase == 'L'), bg = muted(phaseCol[1], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionsx$phase == 'I'), bg = muted(phaseCol[2], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionsx$phase == 'D'), bg = muted(phaseCol[3], c = 40, l = 90))

fx
fx2 <- data.frame(name= c('pops:', 'offspring:', 'loci:', ''),
                  value = c(24,4,24388, NA)) %>% 
  flextable() %>% 
  autofit() %>% 
  delete_part() %>% 
  theme_alafoli() %>%
  bold(j = 1) %>% 
  align(j = 1, align = 'right') %>% 
  align(j = 2, align = 'left') %>% 
  bg(i = 1:3, j = 2, bg = 'grey90')
  
fx2

flextable::save_as_docx(fx2, fx, path = './figures/generation_parameters.docx')
