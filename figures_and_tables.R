
source('./code/libraries.R') # change library locations
source('./functions_figures_tables.R')
# source('./code/functions_conceptual.R')
 source('./code/functions_genetic-data-creation.R')

# load ------------
# dataframes created for plotting and analysing 
phases <- readRDS('./output/phase_df.rds') # population phases

phfst <- readRDS('./output/ph_vis_fst.rds') # ph = P. hermmansburgensis
syfst <- readRDS('./output/sy_vis_fst.rds') # sy = S. youngsoni
                                            # vis = visualise

phalf <- readRDS('./output/ph_vis_alf.rds')
syalf <- readRDS('./output/sy_vis_alf.rds')

phmodel <- readRDS('./output/ph_vis_slopes.rds')
symodel <- readRDS('./output/sy_vis_slopes.rds')

ph <- readRDS('./output/pherm_filtered_genotypes_phases.rds')

phinit <- read.csv('./output/dataframes/initialise_conditions.csv') %>% 
  rename(Nmax = N, N = nInd) %>% 
  mutate(phase = phaseNo, year = gen, species = 'ph')
syinit <- read.csv('./output/dataframes/initialise_conditions_sy.csv')%>% 
  rename(Nmax = N, N = nInd)%>% 
  mutate(phase = phaseNo, year = gen,species = 'sy')

realFst <- data.frame(realfst = c(0.032,0.017, 0.034, 0.027, 0.005, 0.012, 0.02, 0.01,0.03,
                              0.007,0.003, 0.007, 0.003, 0.001, NA,   0.005, 0, 0),
                      lower = c(0.028, 0.016, 0.028,0.006, 0.004, 0.01, 0.017, 0.009, 0.028,
                                0.005, 0.002, NA, -0.001, 0, NA, 0.003, NA, 0),
                      upper = c(0.036, 0.018, 0.04, 0.047, 0.006, 0.015, 0.023, 0.011,0.032,
                                0.008, 0.005, NA, 0.007, 0.001, NA, 0.006, NA, 0),
                      species = rep(c('ph', 'sy'), each = 9),
                      phaseNo =  factor(rep(paste0(rep(c('L', 'I', 'D'), 3),
                                               rep(1:3, each = 3)),2),
                                        levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                        rep(1:3, each = 3))))
realFst

phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])

# analysis ----------
## fst ------------
## regression
phCor <- phfst$fstphase %>% 
  ungroup() %>% 
  group_by(name) %>% 
  mutate(id = 1:160) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  filter(complete.cases(fstch2))

syCor <- syfst$fstphase %>% 
  ungroup() %>% 
  group_by(name) %>% 
  mutate(id = 1:160) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  filter(complete.cases(fstch2))

# numbers
  cor(phCor$fstch2, phCor$fstphase)
  cor(syCor$fstch2, syCor$fstphase)
  lm(fstphase ~ fstch2, data = phCor) %>% summary
  lm(fstphase ~ fstch2, data = syCor) %>% summary
# plot  
  rbind(phCor, syCor) %>% 
    mutate(species = ifelse(species == 'ph', 'P. hermannsburgensis', 
                            'S. youngsoni')) %>% 
  ggplot(aes(fstch2, fstphase, colour = phase))+
    facet_wrap(~species, scale = 'free')+
  scale_colour_manual(values = phaseCol) +
  theme_bw()+
    geom_abline(slope = 1, intercept = 0,
                colour = 'pink', lwd = 2, lty =1)+
  geom_smooth(method = 'lm', colour = 'grey90')+
  geom_point(size = 3, alpha = 0.5)+
    theme(strip.background = element_blank(),
          strip.text = element_text(face = 'italic', size =12),
          panel.grid = element_line(colour = 'grey97'))
## alf -----------
  c(phalf$real$meandiff %>% mean %>% round(.,10),
  phalf$nsim$meandiff %>% mean %>% round(.,10),
  phalf$sim$meandiff %>% mean %>% round(.,10), 
 
  syalf$real$meandiff %>% mean %>% round(.,10),
  syalf$nsim$meandiff %>% mean %>% round(.,10),
  syalf$sim$meandiff %>% mean %>% round(.,10)) %>% 
    plot(., pch = 16, col = c('red','pink', 'grey70')) # index 1-3 is ph, 4-6 is sy 
  abline(h=0, lty = 3, col = 'black')                  # red = empirical data
                                                       # pink = sampled sim
                                                       # grey = full sim
  
  c(phalf$real$vardiff %>% mean %>% round(.,10),
    phalf$nsim$vardiff %>% mean %>% round(.,10),
    phalf$sim$vardiff %>% mean %>% round(.,10), 
    
    syalf$real$vardiff %>% mean %>% round(.,10),
    syalf$nsim$vardiff %>% mean %>% round(.,10),
    syalf$sim$vardiff %>% mean %>% round(.,10)) %>% 
    plot(., pch = 16, col = c('red','pink', 'grey70'))# index 1-3 is ph, 4-6 is sy 
  abline(h=0, lty = 3, col = 'black') 

  
# FIG conceptual -------

## captures --------
con <- em.conceptual_phase('2007-01-01', c(0.24, 0.507,0.68), syoung = F) 
conSy <- em.conceptual_phase('2007-01-01', c(0.24, 0.507,0.68), syoung = T) 
#  theme(legend.position = 'none')


## sim N ------------
### ph -------
conditions <- phfst$conditionstb 
conditions$year <- 0.5
for(i in 2:nrow(conditions)){
  x <- ifelse(conditions$phase[i] == 'L', 1, 0.5)
  conditions$year[i] <- conditions$year[i-1] + x
}


conditionsx <- bind_rows(conditions[1,],phinit,conditions) 
conditionsx <- conditionsx[-1,] %>% 
  left_join(realFst)

simCaps <- ggplot(conditionsx, aes(year, N/48, fill = phase))+
  geom_bar(stat = 'identity', alpha = 0.40) +
  theme_classic() +
  scale_x_continuous(breaks = conditionsx$year, labels = conditionsx$gen)+
  scale_fill_manual(values = c(phaseCol, initialise = 'grey60'),
                    name = 'Population Phase',
                    labels = c('Low', 'Increase', 'Decrease', 'Sim initialised'))+
  theme(legend.position = c(0.91,0.76),
        legend.background = element_rect(colour = 'grey'),
        legend.key.size = unit(0.90, "lines"),
        plot.margin = margin(0.01, 1, 0.01, 0.01, "cm"),
        axis.text.y = element_text(margin = margin(l = 0, r = 2)),
        axis.title.y = element_text(size = 12))+
  xlab('Generation')+
  ylab(expression(italic('N')[e]*' (per subpopulation)')); simCaps

### sy ---------
conditions <- syfst$conditionstb %>% 
  mutate(year = ifelse(year == 0, 0, year - 0.5))
  
conditionsxy <- bind_rows(conditions[1,],syinit,conditions) 
conditionsxy <- conditionsxy[-1,] %>% 
  left_join(realFst)

simCapsSy <- ggplot(conditionsxy, aes(year, N/48, fill = phase))+
  geom_bar(stat = 'identity', alpha = 0.40, width =0.4) +
  theme_classic() +
  scale_x_continuous(breaks = conditionsxy$year, labels = conditionsxy$gen)+
  scale_y_continuous(breaks = c(0,100, 200))+
  scale_fill_manual(values = c(phaseCol, initialise = 'grey60'),
                    name = 'Population Phase',
                    labels = c('Low', 'Increase', 'Decrease', 'Sim initialised'))+
  theme(legend.position = c(0.92,0.8),
        legend.background = element_rect(colour = 'grey'),
        plot.margin = margin(0.01, 1, 0.01, 0.01, "cm"),
        axis.text.y = element_text(margin = margin(l = 0, r = 2)),
        axis.title.y = element_text(size = 12))+
  xlab('Generation')+
 # ylim(0, 250)+
  ylab(expression(italic('N')[e])); simCapsSy

## save FIG -------

fig1 <- grid.arrange(con +theme(legend.position = 'none')+
                       ggtitle(expression('A)'*italic('  P. hermannsburgensis'))), 
             simCaps,
             conSy+ylab('Captures')+
               ggtitle(expression('B)'*italic('  S. youngsoni')))+
               theme(legend.position = 'none',
                     axis.text.y = element_text(margin = margin(l = 14, r = 2))), 
             simCapsSy+theme(legend.position = 'none'),
             ncol = 1, heights = c(2.2,1.75,1.25,0.6))

#ggsave('./figures/fig1_captures_sim.png',fig1,  units = 'cm',
#       height = 21, width = 16, dpi = 300)
# TABLE generations ---------
realFst$realfst %>% sprintf("%f0", .)

### ph ---------
conditionstbxx <- conditionsx[-1,] %>% 
  bind_rows(conditionsxy[-1,])
fx <-conditionstbxx %>% 
  mutate(realfstx = sprintf(realfst, fmt = '%#.3f'),
         N = round(N/48),
         realfstx = paste0(realfstx,
                           '  (', lower, '-',upper, ')'),
         realfstx = ifelse(is.na(lower), realfst, realfstx),
         realfstx = ifelse(is.na(realfst), '', realfstx),
         species = ifelse(species == 'ph', 'P. hermannsburgensis', 
                          'S. youngsoni')) %>% 
  dplyr::select(species, gen, phaseNo, phase, Nmax,  migration, realfstx, N) %>% 
  mutate(#fst = ifelse(duplicated(phaseNo), NA, fst),
         realfstx= ifelse(duplicated(paste(phaseNo, species)), NA, realfstx),
         
         #CIL = ifelse(duplicated(phaseNo), NA, CIL),
         #CIU = ifelse(duplicated(phaseNo), NA, CIU)
         phase = case_when(
           phase == 'Initial' ~  'initialise',
           phase == 'L' ~ 'low',
           phase == 'I' ~ 'increase',
           phase == 'D' ~ 'decrease'
         ),
         phase = ifelse(duplicated(paste(phaseNo,species)), NA, phase),
         Nmax = ifelse(duplicated(paste(phaseNo, species, Nmax)), NA, Nmax),
         migration = ifelse(duplicated(paste(phaseNo,species)), NA, migration),
         phaseNo= ifelse(duplicated(paste(phaseNo,species)), NA, phaseNo),
         species= ifelse(duplicated(species), NA, species)
  ) %>% 
  mutate_if(is.numeric, round, digits = 2) %>% 
  #mutate(N = round(N/48)) %>% 
  #dplyr::select(-phaseNo) %>% 
  #relocate(realfst, .before = N) %>% 
  rename(realfst = realfstx, t= gen,
         Phase = phase, m = migration) -> fx 
  
conditionstb <- conditionstbxx[1:19,]
rbind(names(fx)[-1], fx[1:19,-1]) %>% 
  flextable() %>% 
  autofit() %>% 
  flextable::border_remove() %>% 
  bold(part = 'body',i = 1, j = 1:7) %>%
  italic(i = 1,j = c(1,5,7)) %>%
 # bold(j = c(7,8), i = 1:nrow(conditionsx)) %>% 
  italic(2:20, j =6) %>% 
  # hline(i = which(fx$species == 'S. youngsoni')-1,
  #       border = fp_border(color = 'grey40', width = 3)) %>% 
  hline(i = c(1,20),
        border = fp_border(color = 'grey40', width = 2)) %>% 
  hline(i = (which(!duplicated(conditionstb$phaseNo)))[-1],
        j = 1:7,
        border = fp_border(color = "grey80",
                           width = 1, style = "dashed")) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'L')+1, bg = scales::muted(phaseCol[1], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'I')+1, bg = scales::muted(phaseCol[2], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'D')+1, bg = scales::muted(phaseCol[3], c = 40, l = 90)) %>% 
  compose(
    part = "header", 
    value = as_paragraph(('A) '), as_i('P. hermannsburgensis'),
                         ' simulation parameters')
  ) %>% 
  compose(
    part = "body", j = 6, i = 1,
    value = as_paragraph(('Real '),as_i('F'), as_sub('ST'))
  ) %>% 
  compose(
    part = "body", j = 4, i = 1,
    value = as_paragraph(as_i('N'), ' max')
  ) %>% compose(
    part = "body", j = 2, i = 1,
    value = as_paragraph('Phase', as_sub('id'))
  ) %>% 
  merge_at(part = 'header', j=1:7) %>% 
  align(j = 2, align = 'right') %>% 
  align(j = 4:5, align = 'center') %>% 
  align(j = c(3,6), align = 'left')  ->fxtb;fxtb 
### sy ------------
conditionstb <- conditionstbxx[20:30,]
rbind(names(fx)[-1], fx[20:30,-1]) %>% 
  flextable() %>% 
  autofit() %>% 
  flextable::border_remove() %>% 
  bold(part = 'body',i = 1, j = 1:7) %>%
  italic(i = 1,j = c(1,5,7)) %>%
  # bold(j = c(7,8), i = 1:nrow(conditionsx)) %>% 
  italic(2:12, j =6) %>% 
  # hline(i = which(fx$species == 'S. youngsoni')-1,
  #       border = fp_border(color = 'grey40', width = 3)) %>% 
  hline(i = c(1,12),
        border = fp_border(color = 'grey40', width = 2)) %>% 
  hline(i = (which(!duplicated(conditionstb$phaseNo)))[-1],
        j = 1:7,
        border = fp_border(color = "grey80",
                           width = 1, style = "dashed")) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'L')+1, bg = scales::muted(phaseCol[1], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'I')+1, bg = scales::muted(phaseCol[2], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'D')+1, bg = scales::muted(phaseCol[3], c = 40, l = 90)) %>% 
  
  compose(
    part = "header", 
    value = as_paragraph(('B) '), as_i('S. youngsoni'),
                         ' simulation parameters')
  ) %>% 
  compose(
    part = "body", j = 6, i = 1,
    value = as_paragraph(('Real '),as_i('F'), as_sub('ST'))
  ) %>% 
  compose(
    part = "body", j = 4, i = 1,
    value = as_paragraph(as_i('N'), ' max')
  ) %>% 
  compose(
    part = "body", j = 2, i = 1,
    value = as_paragraph('Phase', as_sub('id'))
  ) %>% 
  merge_at(part = 'header', j=1:7) %>% 
  align(j = 2, align = 'right') %>% 
  align(j = 4:5, align = 'center') %>% 
  align(j = c(3,6), align = 'left') ->fxtbSy;fxtbSy 

## save table -------
# flextable::save_as_docx(fxtb,fxtbSy,
#                         path = './figures/generation_conditions.docx')


# FIG fst ---------------
realFst$name <- 'fstch2'
dfs_phase <- bind_rows(phfst$fstphase,
                       syfst$fstphase) %>% 
  left_join(realFst) %>% 
  mutate(species = ifelse(species == 'ph', 'P. hermannsburgensis',
                          'S. youngsoni'))

dfs_phase %>%   
  ggplot(aes(phaseNo, value, colour = name))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  scale_color_manual(values = c('coral3', 'grey50'),
                     name = NULL,
                     labels = c('Observed', 
                                expression('Simulated')),
                     guide = guide_legend(override.aes = list(size = c(3,3),
                                                              alpha = 1)))+
  scale_size_manual(values = c(3,3), guide = 'none')+
 # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.75) +
  geom_point(aes(size = name), alpha = 0.4, position = position_dodge(0.5))+
  theme_bw()+
  facet_wrap(~species, ncol = 1)+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.85,0.35),
        legend.text.align = 0,
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )+
  ylim(0,0.04)+
  ylab(expression(italic('F')[ST]))+
  xlab('Sampling period') -> figb;figb

#ggsave('./figures/fig3_fst_simulations.png',
#      figb, units = 'cm', height = 14, width = 14, dpi = 600)

# FIG alf -----------

simAlf <- bind_rows(phalf$nsim,syalf$nsim) %>%
  bind_rows(phalf$sim, syalf$sim) %>% 
  mutate(species = ifelse(species== 'sy', 'S. youngsoni',
                          'P. hermannsburgensis'),
         type = ifelse(type == 'Sim', 'Simulation', 'n Simulation'),
         type = factor(type, levels = c('Simulation', 'n Simulation')),
         phasepair = gsub('[[:digit:]]+', '', name),
         phasepair = case_when(
           phasepair == 'D_L' ~ 'Decrease - Low',
           phasepair == 'L_I' ~ 'Low - Increase',
           phasepair == 'I_D' ~ 'Increase - Decrease'
           
         )) %>% 
  filter(type == 'Simulation')

## real var -----------------
  
  realAF <- phalf$real[-1,] %>% 
    bind_rows(syalf$real[-1,]) %>% 
    mutate(name = factor(name, levels =levels(simAlf$name)),
           species = ifelse(species == 'sy', 'S. youngsoni',
                            'P. hermannsburgensis'))
  
  phalf$nsim %>% 
    bind_rows(syalf$nsim) %>% 
  mutate(name = factor(name, levels = levels(simAlf$name)),
         species = ifelse(species == 'sy', 'S. youngsoni',
                          'P. hermannsburgensis')) %>% 
  ggplot(aes(name, vardiff, colour = type, group = i))+
    geom_line(size = 0.3, alpha = 0.2)+
    #geom_line(data = changeSim2, linewidth = 0.2)+
    geom_line(data = realAF,
              size = 1)+
    theme_bw()+
    scale_x_discrete(labels = sub('_', ' - ', levels(simAlf$name)))+
    facet_wrap(~species, ncol = 1,scale = 'free')+
    theme(legend.position = c(0.85,0.90),
          legend.background = element_rect(colour = 'grey'),
          panel.grid = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          legend.text.align = 0,
          strip.text.x = element_text(face = 'italic', size = 11),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
    scale_colour_manual(values = c('grey50', 'coral3'),  #viridis_pal()(10)[c(4,1,7)],
                        labels = c(expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                   'Empirical',
                                   'Simulated'),
                        name = 'Data')+
  #  scale_x_discrete(labels = sub('_', ' - ', changeReal$name[-1]))+
    labs(x ="Population Phases")+
    ylab('Variation of Change in Allele Frquency')-> fig_var;fig_var

  #ggsave('./figures/fig4_alf_phaseNo.png',
  #       fig_var, units = 'cm', height = 16, width = 14, dpi = 300)
  # 
# FIG model slopes ------------------
  
  
  ## slopes ------------------
  dfslopes.raw <- bind_rows(phmodel$real, phmodel$sim) %>% 
    bind_rows(symodel[[1]]) %>%
    bind_rows(symodel[[2]]) %>% 
    filter(pvalue < 0.05,
           # grepl('01', i) | i == '51'
           #R2 > 0.8
    ) %>% 
    mutate(direction = ifelse(slope < 0, 'negative slope', 'positive slope'),
           dslope = slope,
           slope = abs(slope),
           sim_i = round(i),
           species = ifelse(is.na(species), 'S. youngsoni',
                            'P. hermannsburgensis'))
  
  dfslopes <- dfslopes.raw %>% 
    group_by(species, type,  i) %>% 
    summarise(p = mean(pvalue, na.rm = T),
              plower = p - sd(pvalue)/sqrt(n())*2,
              pupper = p + sd(pvalue)/sqrt(n())*2,
              s = mean(slope, na.rm = T),
              slower = s - sd(slope)/sqrt(n())*2,
              supper = s + sd(slope)/sqrt(n())*2,
              r = mean(R2, na.rm = T),
              rlower = r - sd(R2)/sqrt(n())*2,
              rupper = r + sd(R2)/sqrt(n())*2) %>% 
    ungroup() %>% 
    mutate(sim_i = round(i)) %>% 
    group_by(species, type, sim_i) %>% 
    summarise(s = mean(s),
              slower = mean(slower),
              supper = mean(supper))
    
  # R2 is not a good way to find loci acting abnormally - ignore it. 
  # use slope and pvalue...?
  
  ## 3rd quartile ----
  sppdfslopes.raw <- filter(dfslopes.raw, species == 'P. hermannsburgensis')
  phq75 <- mean(tapply(sppdfslopes.raw$slope, sppdfslopes.raw$i,
                function(x) quantile(x, probs = 0.75))[501])
  
  sppdfslopes.raw <- filter(dfslopes.raw, species != 'P. hermannsburgensis')
  syq75 <- mean(tapply(sppdfslopes.raw$slope, sppdfslopes.raw$i,
                  function(x) quantile(x, probs = 0.75))[501])
  
  ## mean slope ---------
  sppdfslopes.raw <- filter(dfslopes.raw, species == 'P. hermannsburgensis')
  mslope <- tapply(sppdfslopes.raw$slope, sppdfslopes.raw$i,
                       mean)
  
  mslope[-length(mslope)] %>% summary
  mslope[length(mslope)]
  
  sppdfslopes.raw <- filter(dfslopes.raw, species != 'P. hermannsburgensis')
  mslope <- tapply(sppdfslopes.raw$slope, sppdfslopes.raw$i,
                       mean)
  mslope[-length(mslope)] %>% summary
  mslope[length(mslope)]
  
  
  ## save FIG ----------
  ggplot(dfslopes, aes(x = sim_i, y = s, colour = type))+
    geom_errorbar(aes(ymin = slower, ymax = supper), width = 0)+
    geom_point(size = 1)+
    facet_grid(rows = vars(species), scale = 'free')+
    theme_bw()+
    scale_color_manual(values = c('grey40', 'coral3'),
                       name = 'Data',
                       labels = c('Simulated', 'Real'))+
    theme(panel.grid = element_blank(),
          plot.margin = margin(l = 1, r = 1, b = 5),
          strip.background = element_blank(),
          strip.text.y = element_text(face = 'italic',size = 11),
          axis.title.x = element_blank(),
        #  axis.text.y = element_text(margin = margin(l = 5, r = 2)),
          axis.text.x = element_blank()
        )+
    ylab(expression('Mean model | '* italic(slope) *' |'))+
    xlab('Simulation number') -> slopesMean
  
  
  ggplot(dfslopes.raw, aes(x = sim_i, y = slope, colour = type,
                           group = sim_i))+
     geom_hline(data = data.frame(species = c('P. hermannsburgensis',
                                             'S. youngsoni'),
                                 slope = c(phq75, syq75),
                                 type = 'real'),
               aes(yintercept = slope), lty = 3, colour = 'coral2')+
    geom_boxplot(fill = 'grey90', width = 0.5, size = 0.3, outlier.size = 0.5)+
    facet_grid(rows = vars(species), scale = 'free')+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(face = 'italic', size = 11),
          #axis.title.y = element_blank(),
          axis.text.y = element_text(margin = margin(l = 0, r = 2)),
          plot.margin = margin(l = 1, r = 1, b = 5))+
    xlab('Simulation number')+
    scale_color_manual(values = c('grey40', 'coral3'),
                       name = 'Data',
                       labels = c('Simulated', 'Real'))+
    ylab(expression('Model | '* italic(slope) *' |')) -> slopesDist;slopesDist
  
  
  yleft <- textGrob(expression('Model | '* italic('β'[2]) *' |'), 
                    rot = 90, gp = gpar(fontsize = 14))
  
  slopeFig <- grid.arrange(slopesMean + theme(legend.position = 'none')+
                             labs(tag = 'A)'),
                           slopesDist + theme(legend.position = 'none')+
                             labs(tag = 'B)'),
                           heights = c(5,6))
  
  #ggsave('./figures/fig_model_slopes_npp.png',
   #      slopeFig, units = 'cm', height = 21, width = 14, dpi = 300)

  # FIG pherm outliers --------------
  sppdfslopes.raw <- filter(dfslopes.raw, species == 'P. hermannsburgensis')
  q75 <- tapply(sppdfslopes.raw$slope, sppdfslopes.raw$i,
                function(x) quantile(x, probs = 0.75))
  
  q75[-length(q75)] %>% summary
  q75[length(q75)]
  
  outlierSlopes <- mean(q75[-length(q75)]) #quantile(q75[-length(q75)], probs = 0.75)
  
  ## extra outliers -----
  sum(phmodel$realL1$pvalue < 0.05& abs(phmodel$realL1$slope) > outlierSlopes)
  sum(phmodel$real$pvalue < 0.05 & abs(phmodel$real$slope) > outlierSlopes)
  sum(phmodel$real$loci[phmodel$real$pvalue < 0.05] %in% phmodel$realL1$loci[phmodel$realL1$pvalue < 0.05])
  sum(phmodel$real$loci[phmodel$real$pvalue < 0.05 & abs(phmodel$real$slope) > outlierSlopes] %in% phmodel$realL1$loci[phmodel$realL1$pvalue < 0.05 & abs(phmodel$realL1$slope) > outlierSlopes])
  
  ## 95% slopes 
  simslopes <- filter(dfslopes.raw, species == 'P. hermannsburgensis',
                      type != 'Real')
 
 # outlierSlopes <- mean(simslopes$slope) + sd(simslopes$slope)*2
  sum(phmodel$real$loci[abs(phmodel$real$slope) > outlierSlopes &
                          phmodel$real$pvalue < 0.05] %in% 
      phmodel$realL1$loci[abs(phmodel$realL1$slope) > outlierSlopes &
                            phmodel$realL1$pvalue < 0.05])
  
    ## data  ------------
  outliers_withoutL1 <- phmodel$real %>% mutate(abslope = abs(slope)) %>% 
    filter(abslope > outlierSlopes, pvalue < 0.05) %>% 
    dplyr::select(loci) %>% unlist
  
  locmeta_ph <- ph@other$loc.metrics[,c(1,28,5,6,9,10, 29)]
  names(locmeta_ph) <- c('AlleleID', 'uid', 'mus_genome', 'mus_position',
                         'Pse_genome', 'Pse_position', 'rdepth')  
  locmeta_ph %>% head
  locmeta <- locmeta_ph %>% 
    mutate(scaffold = sub('HiC_scaffold_', '', Pse_genome),
           scaffold = ifelse(scaffold == '', 0, scaffold),
           scaffold = as.numeric(scaffold)) %>% 
    group_by(scaffold) %>% 
    mutate(position = Pse_position - min(Pse_position))
  
  negslopes <- phmodel$realL1 %>%
    mutate(abslope = abs(slope),
           direction = ifelse(slope < 0, 'negative slope', 'positive slope')) %>%
    filter(abslope > outlierSlopes, pvalue < 0.05, loci %in% outliers_withoutL1) %>% 
    left_join(locmeta) %>% 
    mutate(Pse_genome = factor(Pse_genome),
           mus_genome = factor(mus_genome)) %>% 
    rename(pse_genome = Pse_genome)
  
  dfplot <- phmodel$alfreal %>%
    left_join(negslopes) %>% 
    filter(uid %in% negslopes$uid) %>% 
    filter(phase != 'decrease2') %>%   
    mutate(alf = ifelse(direction == 'negative slope', alf1, alf2),
           phase = factor(phase, levels = c('low',
                                            'increase',
                                            'decrease'))) %>%
    filter(phase != 'decrease2',
           #Pse_position != 0
           ) %>% 
    # filter(R2 > 0.75) %>% 
    arrange(pse_genome) %>% 
    mutate(noloci = 1:nrow(.))
  ## n loci ------
  dfplot$loci %>% unique %>% length
  ## Af phase ------
  phaseColour <- phaseCol
  names(phaseColour) <- c('low','increase','decrease')

    ggplot(dfplot, aes(phaseNo, alf, colour = phase, group = loci)) +
    #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
    geom_line(size = 0.5, alpha = 0.5)+
    #geom_hline(yintercept = c(0,1), colour = 'grey', lty = 3)+
    # geom_boxplot(aes(group = phaseNo, fill = phase), 
    #              colour = 'grey50',
    #              alpha = 0.5, width = 0.1, size = 0.4)+#
    geom_point(size = 1, alpha = 0.5, aes(group = phaseNo))+
    # facet_wrap(~direction, ncol = 1, scale = 'free')+
    scale_colour_manual(values = phaseColour,#rep('grey',3),
                        # labels = c('low to increase', 
                        #            'increase to decrease',
                        #            'decrease to low'),
                        name = 'Phase')+
    guides(alpha='none',
           colour = guide_legend(override.aes = list(linewidth=1)))+
    theme_bw()+
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_text(colour = 'grey50'),
          axis.ticks.length.y = unit(0.3, units = 'cm'),
          axis.line.y = element_line(colour = 'black'),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_blank())+
    ylab('Allele Frequency')+
    xlab('Population phase #') -> fig_af_phaseNo;fig_af_phaseNo

  ## lm npp --------------
  
    ggplot(dfplot, aes(npp, alf, 
               group = loci))+
    # geom_point(size = 1, pch = 21, colour = 'grey50') +
    geom_jitter(aes(colour = phase),size = 0.5, 
                width = 0.01) +
    #geom_line()+
    # geom_vline(xintercept = 0.001)+
    scale_colour_manual(values = c(phaseColour, Models = 'grey'),
                        name = c('Population phase'),
                        guide = guide_legend(override.aes = list(
                          size = c(2,2,2,1),
                          linetype = c(rep("blank", 3), "solid"),
                          shape = c(rep(16, 3), NA))))+
    # scale_colour_manual(values = virid(5)[1:4],
    #                     labels = c('unassigned', 'scaffold 18', 'scaffold 22',
    #                                'scaffold 23'))+
    
    geom_smooth(method = 'lm', se = F, aes(colour = 'Models'), 
                alpha = 0.5, size = 0.5)+
    theme_bw()+
    #facet_wrap(~loci, scale = 'free')+
    scale_x_log10()+
    #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
    theme(panel.grid = element_blank(),
          #legend.position = 'none',
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.text.x = element_text(colour = 'grey50'),
          axis.ticks.length.y = unit(0.3, units = 'cm'),
          axis.line.y = element_line(colour = 'black'),
          axis.line.x = element_line(colour = 'black'),
          strip.background = element_blank(),
          # axis.title.x = element_blank(),
          strip.text = element_blank()
    )+
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
    xlab(expression('NPP '[log-scale]))+
    ylab('Allele Frequency') -> fig_af_npp;fig_af_npp
  
  ## save FIG ------------
  fig_realAF <- grid.arrange(fig_af_npp, 
                       fig_af_phaseNo +
                         theme(legend.position = 'none'), ncol =1)
  
  #ggsave('./figures/fig6_alf_change_real.png',fig_realAF,
  #       units = 'cm', height = 14, width = 12, dpi = 300)
  
    
  # FIG histogram of counts ------
  dfslopes.raw$i[dfslopes.raw$i == 21 & dfslopes.raw$species==
                   'P. hermannsburgensis'] %>% table
  df2 <- dfslopes.raw$i[dfslopes.raw$slope > outlierSlopes]
  
  df2 <-as.data.frame(table(dfslopes.raw$i[dfslopes.raw$slope > outlierSlopes],
                            dfslopes.raw$species[dfslopes.raw$slope > outlierSlopes])) 
    
  
  ## slopes ------------------
  histdata_REALnpp <- as.data.frame(table(dfslopes.raw$i, dfslopes.raw$species)) %>% 
  mutate(type = ifelse(Var1 == 21, 'Real', 'Sim')) %>% 
    rename(Freq1 = Freq) %>% 
    left_join(df2) %>% 
    pivot_longer(cols = c(Freq1, Freq)) %>% 
    mutate(name = ifelse(name == 'Freq1', 
                         'pvalue < 0.05',
                         'pvalue < 0.05 and slope > 0.021')) %>% 
    filter(name == 'pvalue < 0.05')
  
  # Real npp Hist -------------
  
  histdata_REALnpp %>%
    mutate(sim = as.numeric(as.character(Var1))*100,
           div5 = sim %% 3,
           subset = ifelse((div5 == 0),TRUE, FALSE)) %>% 
    filter(subset) %>%
    ggplot(aes(x = value, fill = type, colour = type))+ 
    geom_histogram() + 
    facet_grid(~Var2, scales = 'free_x'#, strip.position = 'bottom'
    )+
    scale_fill_manual(values = c('coral2', 'grey70'),
                      labels = c('Empir.', 'Sim.'),
                      name = 'Data')+
    scale_colour_manual(values = c('coral2', 'black'),
                        labels = c('Empir.', 'Sim.'),
                        name = 'Data')+
    theme_bw()+
    theme(strip.background = element_blank(),
          #  strip.placement = 'outside',
          strip.text.x = element_text(face = 'italic'),
          legend.position = c(0.92,0.83), 
          plot.margin = margin(10,10,10,10),
          legend.key.size = unit(0.3, units = 'cm'),
          legend.text = element_text(size = 8),
          strip.text.y = element_blank(),#element_text(face = 'bold', size = 10),
          legend.background = element_rect(colour = 'grey'),
          panel.grid = element_blank())+
    xlab('Number of loci with p-values < 0.05') +
    ylab('Count')->fig_hist2
  fig_hist2
  
  # ggsave('./figures/fig5_loci_count_hist_real.png',fig_hist2,
  #        units = 'cm', height = 8, width = 14, dpi = 300)
  # 
 
  ## 95% n pvalues --------------
 m <- mean(histdata_REALnpp$value[1:500])
 s <- sd(histdata_REALnpp$value[1:500])
 x <- histdata_REALnpp$value[501]
q <- qnorm(0.95, mean = m, sd = s)
x-q
pnorm(x, mean = m, sd = s)
(x-m)/s
max(histdata_REALnpp$value[1:500])

m <- mean(histdata_REALnpp$value[502:1001])
s <- sd(histdata_REALnpp$value[502:1001])
x <- histdata_REALnpp$value[1002]
qnorm(0.95, mean = m, sd = s)
pnorm(x, mean = m, sd = s)
(x-m)/s

# APPENDIX table 
## table loci position -----------
candidates <- dfplot %>% 
  dplyr::select(loci, slope, mus_genome, mus_position,
                pse_genome, Pse_position, scaffold, position) %>%
  unique(.) %>% 
  filter(#mus_position != 0 | 
           Pse_position != 0)




loctb <- table(candidates$scaffold)#[-1]
locmeta_scaff <- locmeta %>% 
  filter(scaffold %in% as.numeric(names(loctb)),
         uid %in% phmodel$real$uid,
         #mus_position != 0 | 
         Pse_position != 0)
nSNP <- 24388
allloctb <- table(locmeta_scaff$scaffold)

nSNPprob <- (allloctb/nSNP)

em.likelyhood <- function(p, n, obs){
  precision <- 100000
  xx <- sapply(1:precision, function(x) sum(sample(0:1,
                                               n, prob = c(1-p, p), 
                                               replace = T)))
  sum(xx > obs)/precision
}
   
nSNPexp <- round(nSNPprob * nrow(candidates),1)
expdiff <- loctb - nSNPexp
expdiff[expdiff>0]
mean(expdiff)

dfscaff <- as.data.frame(allloctb, responseName = 'Total') %>% 
  rename(Scaffold = Var1) %>% 
    mutate(Prob = nSNPprob, Expected = nSNPexp, Observed = unname(loctb),
           Difference = unname(expdiff))
for (i in 1:nrow(dfscaff)) {
  dfscaff$Proportion[i] <- round(em.likelyhood(dfscaff$Prob[i],
                                         nrow(candidates),
                                         dfscaff$Observed[i]), 3)
}
dfscaff$Corrected <- round(p.adjust(dfscaff$Proportion, method = 'fdr'),2)

  flextable(dfscaff[,-c(2,3,6,7)]) %>% 
  autofit() %>% 
  border_remove() %>% 
    color(i = which(dfscaff$Corrected < 0.05), j = 1,color = 'coral') %>% 
    bold(i =  which(dfscaff$Corrected < 0.05), j = 1) %>% 
    align(j = 1, align = 'center') %>% 
    hline(part = 'header',
    border = fp_border(color = "grey50",
                       width = 3, style = "solid")) %>% 
    hline(i = nrow(dfscaff),
          border = fp_border(color = "grey50",
                             width = 4, style = "solid")) %>% 
    vline(j = 1, part = 'body',
          border = fp_border(color = "grey70",
                             width = 2, style = "solid")) %>% 
  bg(i = which(dfscaff$Corrected < 0.05), j = 2:4, bg = 'mistyrose') %>%
    compose(part = 'header', j = 4,
            value = as_paragraph(('p-value'), as_sub(' FDR'))) ->fxtbScaff;fxtbScaff

### save fx -----
#  flextable::save_as_docx(fxtbScaff,
 #                         path = './figures/scaffold_candidate_loci.docx')

