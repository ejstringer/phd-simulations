
# fst confidence intervals
phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])
dfs <- lapply(list.files('/data/scratch/emily/simulations/meta/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs %>% head

phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
  
 vis<-  data.frame(gen = c(rep(which(!duplicated(dfs$phaseNo)), each = 4)-0.5, 19,19)[-c(1,2)],
             phaseNo = rep(phaseNo, each = 4),
             lower = 0,
             upper = 0.04,
             fst = c(0, rep(rep(c(0.04, 0), each = 2), 8))[-32]) %>% 
   mutate(phase = str_sub(phaseNo, 1,1),
          phase = factor(phase, levels = c('L','I', 'D')),
          phaseNo= factor(phaseNo, 
                 levels = paste0(rep(c('L','I', 'D'), 3), 
                                 rep(1:3, each = 3))[-1]))

vis
dfs %>% 
ggplot(aes(gen, fst, fill = phase))+
  #geom_area(data = vis, aes(x = gen, y = upper), alpha = 0.60)+
  geom_area(data = vis, 
            aes(x = gen, y = fst, group = phaseNo, fill = phase), 
            alpha = 0.60) +
  geom_point(size = 1, #position = position_nudge(x = 0.1, y = 0), 
             alpha = 0.5, colour = 'grey50')+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  geom_line(aes(group = simulation), linewidth = 0.2, 
            alpha = 0.5, colour = 'grey50')+
 # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  
 # geom_boxplot(aes(group = gen), colour = 'black', width = 0.5)+
  #  geom_point(aes(y = fstch2), colour = 'grey60', alpha = 0.5)+
  theme_classic()+
  theme(legend.position = 'bottom')+
  scale_x_continuous(breaks = seq(1,19,2))+
  scale_fill_manual(values = phaseCol,
                     name = 'Population Phase',
                     labels = c('low', 'increase', 'decrease'))+
  ylim(0,0.04)+
  ylab('Fst')+
  xlab('generation') -> figa;figa 
#ggsave('./test_fig.png')

dfs %>% group_by(phaseNo, phase, fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) -> dfs_phase
  #filter(!duplicated(name) | name == 'fstphase') %>% 
dfs_phase %>%   
ggplot(aes(phaseNo, value, colour = name))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  scale_color_manual(values = c('coral3', 'grey50'),
                     name = '',
                     labels = c('Observed Fst', 'Simulated Fst'),
                     guide = guide_legend(override.aes = list(size = c(3,3),
                                                              alpha = 1)))+
  scale_size_manual(values = c(3,3), guide = 'none')+
  geom_point(aes(size = name), alpha = 0.91)+
  theme_classic()+
  #guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = 'bottom' #c(0.8,0.9),
        # legend.background = element_rect(colour = 'grey50'),
        # legend.key.size = unit(0.25, units = 'cm'),
        # legend.box.spacing = unit(0.1, units = 'cm'),
        # legend.text = element_text(size = 6),
        # legend.title = element_text(size = 8)
  )+
  ylim(0,0.04)+
  ylab('Fst')+
  xlab('population phase #') -> figb;figb

figab<- grid.arrange(figa, figb, ncol = 1)  
ggsave('./figures/fig3_fst_simulations2.png',
       figab, units = 'cm', height = 20, width = 16, dpi = 300)


# alternate ---
dfs %>% group_by(phaseNo, phase, fstch2, gen) %>% 
  summarise(fstgen = mean(fst),
            n = n(),
            sd = sd(fst),
            se = sd/n,
            lower = fstgen - sd*1.96,
            upper = fstgen + sd*1.96,
            grp = 'cool') %>% 
  ggplot(aes(gen, fstgen, colour = phase, fill = phase))+
  geom_line(aes(group = grp), linewidth = 0.75)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  geom_point(size = 1.5, pch = 21, colour = 'black')+
  #  geom_point(aes(y = fstch2), colour = 'grey60', alpha = 0.5)+
  theme_classic()+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  theme(legend.position = 'none')+
  scale_x_continuous(breaks = seq(1,19,2))+
  scale_fill_manual(values = phaseCol)+
  scale_color_manual(values = phaseCol)+
  ylim(0,0.04)+
  ylab('Fst')+
  xlab('generation')

# table -----
conditions  <- lapply(list.files('/data/scratch/emily/simulations/meta/', full.names = T),
                      read.csv) %>% 
  do.call('rbind', .)

conditionsx <-conditions %>% rename(Nmax = N) %>% 
  group_by(gen,year, phaseNo, phase, Nmax, migration) %>% 
  summarise(N = mean(nInd),
            fstgen = mean(fst),
            cil = fstgen - (sd(fst)/n())*1.96,
            ciu = fstgen + (sd(fst)/n())*1.96)

conditionsxphase <- conditions %>% 
  filter(!duplicated(paste0(simulation, phaseNo))) %>% 
  group_by(phaseNo) %>% 
  summarise(fst = mean(fstphase),
            CIL = fst - (sd(fstphase)/n())*1.96,
            CIU = fst+ (sd(fstphase)/n())*1.96)

conditionsTB <- left_join(conditionsx, conditionsxphase) %>% as.data.frame

conditionsTB %>% 
  mutate(fst = ifelse(duplicated(phaseNo), NA, fst),
         CIL = ifelse(duplicated(phaseNo), NA, CIL),
         CIU = ifelse(duplicated(phaseNo), NA, CIU))
write.csv(conditionsTB, './output/dataframes/ph_sim_conditons.csv', row.names = F) 

# savedata -------
dfs$species <- 'ph'
vis$species <- 'ph'
dfs_phase$species <-  'ph'
conditionsTB$species <-  'ph'
conditionsTB$year <- conditionsTB$gen
saveRDS(list(fstgen = dfs, phasevis = vis, fstphase = dfs_phase,
             conditionstb = conditionsTB),
        './output/ph_vis_fst.rds')
