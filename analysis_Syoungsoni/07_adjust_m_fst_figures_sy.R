# fst confidence intervals
phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])
list.files('/data/scratch/emily/simulations/sy/fst_mAdjust/')
dfs <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

dfs_N50 <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_N200_mEquation/',
                         full.names = T), read.csv) %>% 
  do.call('bind_rows', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))

  dfs_N30 <- lapply(list.files('/data/scratch/emily/simulations/sy/meta_N200_mAdjust/',
                             full.names = T), read.csv) %>% 
  do.call('rbind', .) %>% 
  mutate(phaseNo = factor(phaseNo, 
                          levels = paste0(rep(c('L','I', 'D'), 3), 
                                          rep(1:3, each = 3))[-1]))


dfs %>% head
dfs$N %>% table
dfs %>% group_by(phaseNo, migration) %>% 
  summarise(N = mean(N))
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]


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
dfs_phaseData$name %>% table()
#filter(!duplicated(name) | name == 'fstphase') %>% 
dfs_phaseData %>% 
  ggplot(aes(phaseNo, value, colour = name))+
 # geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
#  geom_hline(yintercept = 0.034, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  
  scale_color_manual(values = c('coral3', 'yellow','orange', 'red'),
                     name = '',
                     labels = c('Observed', 'Simulated N 100', 'Simulated N 200','Simulated m adjusted (200)'),
                     guide = guide_legend(override.aes = list(size = c(3,4,3,2),
                                                              alpha = 1)))+
  scale_size_manual(values = c(3,4,3,2), guide = 'none')+
  scale_alpha_manual(values = c(0.1, 0.2,0.2, 0.3), guide = 'none')+
  geom_point(aes(size = name, alpha = name), position = position_dodge(0.5))+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = 'bottom', #c(0.8,0.9),
        panel.grid = element_blank()
        # legend.background = element_rect(colour = 'grey50'),
        # legend.key.size = unit(0.25, units = 'cm'),
        # legend.box.spacing = unit(0.1, units = 'cm'),
        # legend.text = element_text(size = 6),
        # legend.title = element_text(size = 8)
  )+
  #ylim(0,0.04)+
  ylab('Fst')+
  xlab('Sampling phase #') -> figb;figb
  
  dfs_N30 %>% 
   #filter(!phaseNo %in% c('I1')) %>%
    mutate(migration = round(migration, 3)) %>% 
  ggplot(aes(migration, fstphase, colour = phase))+
    geom_hline(aes(yintercept = fstch2), colour = 'red', lty=2)+
    scale_colour_manual(values = phaseCol)+
    geom_line()+
    theme_bw()+
    facet_grid(rep~phase, scale = 'free')

  dfs_N30 %>% 
    group_by(phaseNo) %>% 
    mutate(close = abs(fstch2 - fstphase),
           minclose = close - min(close),
           select = ifelse(minclose == 0, TRUE, FALSE)) %>% 
    filter(select) %>% 
    select(phaseNo,phase, migration, fstphase, fstch2, close, minclose) %>%
    unique() %>% 
    arrange(phaseNo)
  