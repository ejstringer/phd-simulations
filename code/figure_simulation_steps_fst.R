# fst confidence intervals
phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])

phaseNo <- paste0(rep(c('L','I', 'D'), 3),
                  rep(1:3, each = 3))[-1]

# load ph files ------------
dfs <- lapply(list.files('/data/scratch/emily/simulations/meta_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N50 <- lapply(list.files('/data/scratch/emily/simulations/meta_mAdjust/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N30 <- lapply(list.files('/data/scratch/emily/simulations/meta_decreaseN30/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

## explore -----
dfs %>% head
dfs$N %>% table
dfs %>% group_by(phaseNo, migration) %>% 
  summarise(N = mean(N))

## summarise ----------

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

## join data ----------
dfs_phaseDataPh <- bind_rows(dfs_phase, dfs_N50_phase) %>% 
  bind_rows(dfs_N30_phase) %>% 
  mutate(species = 'ph')


dfs_phaseDataPh$name %>% table()


# load dunnart ------

dfsSy <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N200 <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_N200_mEquation/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N200m <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_N200_mAdjust/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

## summarise -----------
dfsSy %>% group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) -> dfsSy_phase

dfs_N200 %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase2') -> dfs_N200_phase

dfs_N200m %>%
  group_by(phaseNo, phase, migration,fstch2, fstphase) %>% 
  summarise(n = n()) %>% 
  pivot_longer(cols = c(fstch2, fstphase)) %>% 
  filter(name != 'fstch2') %>% 
  mutate(name = 'fstphase3') -> dfs_N200m_phase

## join ----------
dfs_phaseDataSy <- bind_rows(dfsSy_phase, dfs_N200_phase) %>% 
  bind_rows(dfs_N200m_phase) %>% 
  mutate(species = 'sy')
dfs_phaseDataSy$name %>% table()

# join SPECIES data------

dfs_phaseData <- bind_rows(dfs_phaseDataPh, dfs_phaseDataSy) %>% 
  mutate(species = ifelse(species == 'ph', 
                          'Pseudomys hermannsburgensis',
                          'Sminthopsis youngsoni'))



# PLOt -------

dfs_phaseData %>% 
  ggplot(aes(phaseNo, value, colour = name))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  scale_color_manual(values = c('black', 'yellow','orange', 'red'),
                     name = 'Scenario',
                     labels = c('Observed', 'Simulated step 1',
                                'Simulated step 2','Simulated step 3'),
                     guide = guide_legend(override.aes = list(size = c(3,4,3,2)-0.5,
                                                              alpha = 1)))+
  facet_wrap(~species, ncol = 1)+
  scale_size_manual(values = c(3,4,3,2)-0.5, guide = 'none')+
  scale_alpha_manual(values = c(0.1, 0.2,0.2, 0.3), guide = 'none')+
  geom_point(aes(size = name, alpha = name), position = position_dodge(0.5))+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.81,0.30),
        panel.grid = element_blank(),
        legend.text.align = 0,
         legend.background = element_rect(colour = 'grey50'),
        strip.text.x = element_text(face = 'italic', size = 12),
        strip.background = element_blank()
        # legend.key.size = unit(0.25, units = 'cm'),
        # legend.box.spacing = unit(0.1, units = 'cm'),
        # legend.text = element_text(size = 6),
        # legend.title = element_text(size = 8)
  )+
  ylim(0,0.04)+
  ylab(expression(italic('F')[ST]))+
  xlab('Sampling phase #') -> figb;figb

ggsave('./figures/figS1_fst_simulation_steps.png',
       figb, units = 'cm', height = 14, width = 12, dpi = 600)
