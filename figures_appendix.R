# Sup figure 


#setup ------------------
phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]



# phermann --------------
dfs <- lapply(list.files('./simulations/meta_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N50 <- lapply(list.files('./simulations/meta_mAdjust/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N30 <- lapply(list.files('./simulations/meta_decreaseN30/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))




dfs %>% group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) -> dfs_phase

dfs_N50 %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase2') -> dfs_N50_phase

dfs_N30 %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase3') -> dfs_N30_phase

dfs_phaseData <- bind_rows(dfs_phase, dfs_N50_phase) %>% 
  bind_rows(dfs_N30_phase)
dfs_phaseData_Phermann <- dfs_phaseData
#filter(!duplicated(name) | name == 'fstphase') %>% 


# Syoung--------------------

dfs <- lapply(list.files('./simulations/sy/meta_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N50 <- lapply(list.files('./simulations/sy/meta_N200_mEquation/',
                             full.names = T), read.csv) %>% 
  do.call('bind_rows', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N30 <- lapply(list.files('./simulations/sy/meta_N200_mAdjust/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))


dfs %>% group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) -> dfs_phase

dfs_N50 %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase2') -> dfs_N50_phase

dfs_N30 %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase3') -> dfs_N30_phase

dfs_phaseData <- bind_rows(dfs_phase, dfs_N50_phase) %>% 
  bind_rows(dfs_N30_phase)

dfs_phaseData_Syoung <- dfs_phaseData

# plot -------------------


dfs_phaseData_Phermann$species <- 'Pseudomys hermannsburgensis'
dfs_phaseData_Syoung$species <- 'Sminthopsis youngsoni'

dfs_phaseData_Phermann %>% 
  bind_rows(dfs_phaseData_Syoung) %>% 
  ggplot(aes(phaseNo, value, colour = name))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  
  scale_color_manual(values = c('black', 'yellow2','orange', 'red'),
                     name = 'Scenario',
                     labels = c('Observed', 'Simulated step 1', 'Simulated step 2','Simulated step 3'),
                     guide = guide_legend(override.aes = list(size = c(3,4,3,2),
                                                              alpha = 1)))+
  scale_size_manual(values = c(3,4,3,2), guide = 'none')+
  scale_alpha_manual(values = c(0.1, 0.2,0.2, 0.3), guide = 'none')+
  facet_wrap(~species, ncol = 1)+
  geom_point(aes(size = name, alpha = name), position = position_dodge(0.5))+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.81,0.3),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = 'italic', size =12),
        legend.background = element_rect(colour = 'grey50'),
      #  legend.key.size = unit(0.25, units = 'cm'),
       # legend.box.spacing = unit(0.1, units = 'cm'),
    #    legend.text = element_text(size = 6),
     #   legend.title = element_text(size = 8)
  )+
  ylim(0,0.04)+
  ylab(expression(italic('F')[ST]))+
  xlab('Sampling period') -> figb;figb

#ggsave('./figures/figS1_fst_simulation_scenarios.png',
#       figb, units = 'cm', height = 14, width = 12, dpi = 300)
