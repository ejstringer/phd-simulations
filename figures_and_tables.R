
source('./code/libraries.R')
# source('./code/functions_conceptual.R')
 source('./code/functions_genetic-data-creation.R')

# load ------------
phases <- readRDS('./output/phase_df.rds')
phfst <- readRDS('./output/ph_vis_fst.rds')
syfst <- readRDS('./output/sy_vis_fst.rds')
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
# fig 1 conceptual -------

con <- em.conceptual_phase('2007-01-01', c(0.24, 0.507,0.68), syoung = F) 
conSy <- em.conceptual_phase('2007-01-01', c(0.24, 0.507,0.68), syoung = T) 
#  theme(legend.position = 'none')

phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])

## sim gener ---------
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
  theme(legend.position = c(0.91,0.75),
        legend.background = element_rect(colour = 'grey'),
        plot.margin = margin(0.01, 1, 0.01, 0.01, "cm"),
        axis.text.y = element_text(margin = margin(l = 0, r = 2)),
        axis.title.y = element_text(size = 12))+
  xlab('Generation')+
  ylab(expression(italic('N')[e]*' (per subpopulation)')); simCaps

## sy sim gener ---------
conditions <- syfst$conditionstb %>% 
  mutate(year = ifelse(year == 0, 0, year - 0.5))
  
conditionsxy <- bind_rows(conditions[1,],phinit,conditions) 
conditionsxy <- conditionsxy[-1,] %>% 
  left_join(realFst)

simCapsSy <- ggplot(conditionsxy, aes(year, N/48, fill = phase))+
  geom_bar(stat = 'identity', alpha = 0.40, width =0.4) +
  theme_classic() +
  scale_x_continuous(breaks = conditionsxy$year, labels = conditionsxy$gen)+
  scale_y_continuous(breaks = c(0,50,100))+
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

# concept fig -------

fig1 <- grid.arrange(con +theme(legend.position = 'none')+
                       ggtitle(expression('A)'*italic('  P. hermannsburgensis'))), 
             simCaps,
             conSy+ylab('Captures')+
               ggtitle(expression('B)'*italic('  S. youngsoni')))+
               theme(legend.position = 'none',
                     axis.text.y = element_text(margin = margin(l = 14, r = 2))), 
             simCapsSy+theme(legend.position = 'none'),
             ncol = 1, heights = c(2.2,1.75,1.25,0.6))

ggsave('./figures/fig1_captures_sim.png',fig1,  units = 'cm',
       height = 21, width = 16, dpi = 300)
# table 1 ---------
realFst$realfst %>% sprintf("%f0", .)


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
         Nmax = ifelse(duplicated(paste(phaseNo, species)), NA, Nmax),
         migration = ifelse(duplicated(paste(phaseNo,species)), NA, migration),
         species= ifelse(duplicated(species), NA, species)
  ) %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  #mutate(N = round(N/48)) %>% 
  dplyr::select(-phaseNo) %>% 
  #relocate(realfst, .before = N) %>% 
  rename(realfst = realfstx, t= gen,
         Phase = phase, m = migration) -> fx 
  
conditionstb <- conditionstbxx[1:19,]
rbind(names(fx)[-1], fx[1:19,-1]) %>% 
  flextable() %>% 
  autofit() %>% 
  flextable::border_remove() %>% 
  bold(part = 'body',i = 1, j = 1:6) %>%
  italic(i = 1,j = c(1,4,6)) %>%
 # bold(j = c(7,8), i = 1:nrow(conditionsx)) %>% 
  italic(2:20, j =5) %>% 
  # hline(i = which(fx$species == 'S. youngsoni')-1,
  #       border = fp_border(color = 'grey40', width = 3)) %>% 
  hline(i = c(1,20),
        border = fp_border(color = 'grey40', width = 2)) %>% 
  hline(i = (which(!duplicated(conditionstb$phaseNo)))[-1],
        j = 2:6,
        border = fp_border(color = "grey80",
                           width = 1, style = "dashed")) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'L')+1, bg = muted(phaseCol[1], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'I')+1, bg = muted(phaseCol[2], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'D')+1, bg = muted(phaseCol[3], c = 40, l = 90)) %>% 
  compose(
    part = "header", 
    value = as_paragraph(('A) '), as_i('P. hermannsburgensis'),
                         ' simulation parameters')
  ) %>% 
  compose(
    part = "body", j = 5, i = 1,
    value = as_paragraph(('Real '),as_i('F'), as_sub('ST'))
  ) %>% 
  compose(
    part = "body", j = 3, i = 1,
    value = as_paragraph(as_i('N'), ' max')
  ) %>% 
  merge_at(part = 'header', j=1:6) %>% 
  align(j = 2, align = 'right') %>% 
  align(j = 3:4, align = 'center') %>% 
  align(j = 5, align = 'left')  ->fxtb;fxtb 

conditionstb <- conditionstbxx[20:30,]
rbind(names(fx)[-1], fx[20:30,-1]) %>% 
  flextable() %>% 
  autofit() %>% 
  flextable::border_remove() %>% 
  bold(part = 'body',i = 1, j = 1:6) %>%
  italic(i = 1,j = c(1,4,6)) %>%
  # bold(j = c(7,8), i = 1:nrow(conditionsx)) %>% 
  italic(2:12, j =5) %>% 
  # hline(i = which(fx$species == 'S. youngsoni')-1,
  #       border = fp_border(color = 'grey40', width = 3)) %>% 
  hline(i = c(1,12),
        border = fp_border(color = 'grey40', width = 2)) %>% 
  hline(i = (which(!duplicated(conditionstb$phaseNo)))[-1],
        j = 2:6,
        border = fp_border(color = "grey80",
                           width = 1, style = "dashed")) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'L')+1, bg = muted(phaseCol[1], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'I')+1, bg = muted(phaseCol[2], c = 40, l = 90)) %>% 
  bg(j = 1:5, i = which(conditionstb$phase == 'D')+1, bg = muted(phaseCol[3], c = 40, l = 90)) %>% 
  
  compose(
    part = "header", 
    value = as_paragraph(('B) '), as_i('S. youngsoni'),
                         ' simulation parameters')
  ) %>% 
  compose(
    part = "body", j = 5, i = 1,
    value = as_paragraph(('Real '),as_i('F'), as_sub('ST'))
  ) %>% 
  compose(
    part = "body", j = 3, i = 1,
    value = as_paragraph(as_i('N'), ' max')
  ) %>% 
  merge_at(part = 'header', j=1:6) %>% 
  align(j = 2, align = 'right') %>% 
  align(j = 3:4, align = 'center') %>% 
  align(j = 5, align = 'left') ->fxtbSy;fxtbSy 

flextable::save_as_docx(fxtb,fxtbSy,
                        path = './figures/generation_conditions.docx')



# FST fig ---------------
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
                     labels = c('Observed Fst', 'Simulated Fst'),
                     guide = guide_legend(override.aes = list(size = c(3,3),
                                                              alpha = 1)))+
  scale_size_manual(values = c(3,3), guide = 'none')+
 # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 0.75) +
  geom_point(aes(size = name), alpha = 0.1)+
  theme_bw()+
  facet_wrap(~species, ncol = 1)+
  #guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(#legend.position = 'bottom', #c(0.8,0.9),
        # legend.background = element_rect(colour = 'grey50'),
        # legend.key.size = unit(0.25, units = 'cm'),
        # legend.box.spacing = unit(0.1, units = 'cm'),
        # legend.text = element_text(size = 6),
        # legend.title = element_text(size = 8)
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.85,0.35),
        legend.background = element_rect(colour = 'grey'),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )+
  ylim(0,0.04)+
  ylab(expression(italic('F')[ST]))+
  xlab('population phase #') -> figb;figb

#figab<- grid.arrange(figa, figb, ncol = 1)  
ggsave('./figures/fig3_fst_simulations.png',
       figb, units = 'cm', height = 14, width = 14, dpi = 600)

# ALF -----------
phalf <- readRDS('./output/ph_vis_alf.rds')
syalf <- readRDS('./output/sy_vis_alf.rds')

simAlf <- bind_rows(phalf$sim,syalf$sim) %>%
  mutate(species = ifelse(is.na(species), 'S. youngsoni',
                          'P. hermannsburgensis'),
         phasepair = gsub('[[:digit:]]+', '', name),
         phasepair = case_when(
           phasepair == 'D_L' ~ 'Decrease - Low',
           phasepair == 'L_I' ~ 'Low - Increase',
           phasepair == 'I_D' ~ 'Increase - Decrease'
           
         ))
## phase -----------
ggplot(simAlf, aes(name, meandiff,fill = species)) +
  theme_bw()+
  geom_hline(yintercept = 0,lty = 2, colour = 'grey50')+
  geom_boxplot()+
  guides(colour = 'none')+
  facet_wrap(~species, scale = 'free_x', ncol = 2)+
  scale_colour_manual(values = rep('black', 7))+
  scale_fill_manual(values = virid(10)[c(7, 5, 10)])+
  theme(legend.position = c(0.85,0.85),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
       # axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.00001,
                                   decimal.mark = '.'))+
  ylab(expression('Mean Change in AF'))+
  xlab('Sampling Phase') -> figMean


ggplot(simAlf, aes(name, vardiff, fill = species)) +
  theme_bw()+
  geom_hline(yintercept = 0,lty = 2, colour = 'grey50')+
  geom_boxplot()+
  guides(colour = 'none')+
facet_wrap(~species, scale = 'free_x', ncol = 2)+
  scale_colour_manual(values = rep('black', 7))+
  scale_fill_manual(values = virid(10)[c(7, 5, 10)])+
  theme(legend.position = 'none',#c(0.85,0.8),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
        #axis.title = element_text(size = 13),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.y = element_text(margin = margin(l = 3, r = 2))
        )+
  ylab(expression('Variance of Change in AF'))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.00001,
                                   decimal.mark = '.'))+
  xlab('Sampling Phase') -> figVar

grid.arrange(figMean+theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           legend.position = 'none'),
             figVar + theme(strip.text.x = element_blank()))  
## species ----
ggplot(simAlf, aes(species, meandiff, fill = species)) +
  theme_bw()+
  geom_hline(yintercept = 0,lty = 2, colour = 'grey50')+
  geom_boxplot()+
  guides(colour = 'none')+
  #facet_wrap(~species, scale = 'free_x', ncol = 2)+
  scale_colour_manual(values = rep('black', 7))+
  scale_fill_manual(values = virid(10)[c(7, 5, 10)])+
  theme(legend.position = 'none',#c(0.85,0.8),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
        #axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_text(face = 'italic'),
        axis.text.y = element_text(margin = margin(l = 3, r = 2))
  )+
  ylab(expression('Mean Change in AF'))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.00001,
                                   decimal.mark = '.'))+
  xlab('Sampling Phase') -> meanAF 
  ggplot(simAlf, aes(species, vardiff, fill = species)) +
  theme_bw()+
  geom_hline(yintercept = 0,lty = 2, colour = 'grey50')+
  geom_boxplot()+
  guides(colour = 'none')+
  #facet_wrap(~species, scale = 'free_x', ncol = 2)+
  scale_colour_manual(values = rep('black', 7))+
  scale_fill_manual(values = virid(10)[c(7, 5, 10)])+
  theme(legend.position = 'none',#c(0.85,0.8),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(colour = 'grey'),
        #axis.title = element_text(size = 13),
        axis.title.x = element_blank(),
        strip.text.x = element_text(face = 'italic', size = 11),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.text.x = element_text(face = 'italic'),
        axis.text.y = element_text(margin = margin(l = 7, r = 2))
  )+
  ylab(expression('Variance of Change in AF'))+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.00001,
                                   decimal.mark = '.'))+
  xlab('simulation')->varAF

  figAF <- grid.arrange(meanAF, varAF, ncol = 1, heights = c(2,3))
  
  ggsave('./figures/fig4_AF_simulations.png',
         figAF, units = 'cm', height = 14, width = 10, dpi = 600)
  
# AF real -----------------