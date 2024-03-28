
em.alf_npp_model <- function(a1, iteration){
  #a1 <- alfdf %>% filter(loci == alfdf$loci[2])
  
  m <- lm(alf2 ~ npp.log, a1)
  m_summary <- summary(m)
  
  
  data.frame(loci = a1$loci[1], uid = a1$uid[1], i = iteration, R2 = m_summary$r.squared,
             intercept = coef(m_summary)[1,1], slope = coef(m_summary)[2,1], 
             df = m_summary$df[2], tvalue =  coef(m_summary)[2,3],
             pvalue =  coef(m_summary)[2,4])
  
}

# npp random -----------
meta <- ph@other$ind.metrics %>% 
  group_by(period, phaseNo, phase, npp) %>% 
  summarise(n = n(),
            ngrids = length(unique(gridId))) %>% 
  rename(nppTrue = npp)
random <- meta[,c('phaseNo', 'nppTrue')]  %>% unique()

set.seed(444)
random$npp <- sample(random$nppTrue, 9); cor(random$nppTrue, random$npp)
plot(log(npp)~log(nppTrue), data = random)
random
meta <- left_join(meta, random)

# real -------
alfReal <- em.alf_df(ph_phase = ph_phase, meta)
alfReal
alfdf <- filter(alfReal, phaseNo != 'L1')

alfloci <- split(alfdf, alfdf$loci)
system.time(m.alfnpp <- lapply(alfloci, em.alf_npp_model, iteration = 'real') %>% 
  do.call('rbind', .))
rownames(m.alfnpp) <-  NULL
m.alfnpp$type <- 'Real'
nrow(m.alfnpp)

# cores- -----
ncores <- 20
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

dir.create('/data/scratch/emily/simulations/sy/models')

# sims m -------
models_nSim <- foreach(i = 1:20) %dopar% {
  library(dartR)
  library(tidyverse)
  print(paste('sim', i, '....Go!!'))
  fileload <- list.files('/data/scratch/emily/simulations/sy/genlights',
                         full.names = T)[i]
  simx <- readRDS(fileload)[12:19]
  nSim_list <- list()
  for(j in 1:25) { # can change to more
    nsim<- em.sample_sim(samplex, simx)
    nSimalf <- em.alf_df(nsim, meta)
    
    alfloci <- split(nSimalf, nSimalf$loci)
    m.alfnpp <- lapply(alfloci, em.alf_npp_model, iteration = j/100) %>% 
      do.call('rbind', .)
    nSim_list[[j]] <- m.alfnpp 
  }
  
  to_save_nSim <- do.call('bind_rows', nSim_list)
  to_save_nSim$i <- to_save_nSim$i + i
  to_save_nSim$type <- 'nSim'
  rownames(to_save_nSim) <- NULL
  mfilename <- sub('\\.rds', '_alf_models.rds', fileload)
  mfilename <- sub('genlights', 'models', mfilename)
  saveRDS(to_save_nSim, mfilename) 
  to_save_nSim
}

stopCluster(cl)

m.alfsim <- do.call('rbind', models_nSim)
nrow(m.alfsim)
m.alfnpp$i <- 1

m.alf <- bind_rows(m.alfnpp, m.alfsim) %>% 
  filter(pvalue < 0.05) %>% 
  pivot_longer(cols = c(R2,slope, pvalue)) %>% 
  mutate(name = factor(name, levels = c('slope', 'R2', 'pvalue')))
names(m.alf)
m.alf$type %>% table
m.alfsim$type
g1 <- ggplot(filter(m.alf, name != 'R2'), aes(x = name, y = value))+
 # geom_histogram(colour = 'black', fill = '#FF9999')+
  geom_boxplot(aes(colour = type))+
  #facet_grid(type~name, scale = 'free')+
  facet_grid(~name, scale = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab('Frequency')+
  xlab('Estimate')
g1  

dfslopes <- bind_rows(m.alfnpp, m.alfsim) %>% 
  filter(pvalue < 0.05) %>% 
  mutate(direction = ifelse(slope < 0, 'neg', 'pos')) %>% 
  group_by(direction, type, i) %>% 
  summarise(p = mean(pvalue, na.rm = T),
            plower = p - sd(pvalue),
            pupper = p + sd(pvalue),
            s = mean(slope, na.rm = T),
            slower = s - sd(slope),
            supper = s + sd(slope)) 

ggplot(dfslopes, aes(x = type, y = p))+
  geom_errorbar(aes(ymin = plower, ymax = pupper), width = 0)+
  geom_point()+
  facet_grid(~direction, scale = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
m.alfnpp$i <- 21


# fig 5 ------------------
## save data -------
dfslopes.raw <- bind_rows(m.alfnpp, m.alfsim) %>% 
  filter(pvalue < 0.05,
         #grepl('01', i) | i == '51'
         #R2 > 0.8
  ) %>% 
  mutate(direction = ifelse(slope < 0, 'negative slope', 'positive slope'),
         dslope = slope,
         slope = abs(slope),
         sim_i = round(i))

saveRDS(list(m.alfnpp, m.alfsim), './output/sy_vis_random.rds')


dfslopes <- bind_rows(m.alfnpp, m.alfsim) %>% 
  filter(pvalue < 0.05,
         grepl('01', i) | i == '51'
         #R2 > 0.8
         ) %>% 
  mutate(direction = ifelse(slope < 0, 'negative slope', 'positive slope'),
         dslope = slope,
         slope = abs(slope),
         sim_i = round(i)) %>%
  group_by(type,  i) %>% 
  summarise(p = mean(pvalue, na.rm = T),
            plower = p - sd(pvalue)/sqrt(n())*2,
            pupper = p + sd(pvalue)/sqrt(n())*2,
            s = mean(slope, na.rm = T),
            slower = s - sd(slope)/sqrt(n())*2,
            supper = s + sd(slope)/sqrt(n())*2,
            r = mean(R2, na.rm = T),
            rlower = r - sd(R2)/sqrt(n())*2,
            rupper = r + sd(R2)/sqrt(n())*2) 
  # R2 is not a good way to find loci acting abnormally - ignore it. 
  # use slope and pvalue...?


ggplot(dfslopes, aes(x = i, y = s, colour = type))+
  geom_errorbar(aes(ymin = slower, ymax = supper), width = 0)+
  geom_point(size = 1)+
  #facet_grid(~direction, scale = 'free')+
  theme_bw()+
  scale_color_manual(values = c('grey40', 'coral3'),
                     name = 'Data',
                     labels = c('Simulated', 'Real'))+
  theme(panel.grid = element_blank(),
        plot.margin = margin(l = 1, r = 1, b = 1),
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  ylab(expression('Mean model | '* italic(slope) *' |'))+
  xlab('Simulation number') -> g5


ggplot(dfslopes.raw, aes(x = sim_i, y = slope, colour = type,
                         group = sim_i))+
  geom_boxplot(fill = 'grey90', width = 0.5, size = 0.3, outlier.size = 0.5)+
 # facet_grid(~direction, scale = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        #axis.title.y = element_blank(),
        plot.margin = margin(l = 5.5, r = 1))+
  xlab('Simulation number')+
  scale_color_manual(values = c('grey40', 'coral3'),
                     name = 'Data',
                     labels = c('Simulated', 'Real'))+
  ylab(expression('Model | '* italic(slope) *' |')) -> g4


yleft <- textGrob(expression('Model | '* italic('β'[2]) *' |'), 
                  rot = 90, gp = gpar(fontsize = 14))

slopeFig <- grid.arrange(g5, g4)

ggsave('./figures/fig5_model_slopes_npp.png',
       slopeFig, units = 'cm', height = 14, width = 18, dpi = 300)

sum(models_nSim[[1]]$pvalue < 0.05, na.rm = T)
names(m.alf)
m.alf <- bind_rows(m.alfnpp, m.alfsim)
m.alf.summary <- m.alf %>% group_by(type, i) %>% 
  summarise(p05 = sum(pvalue < 0.05, na.rm = T),
          #  R2strong = sum(R2 > 0.9, na.rm = T),
            steep = sum(slope > 0.04 | slope < -0.04, na.rm = T))



my_labeller <- as_labeller(c(p05="p < 0.05", R2strong="R^2 > 0.9", steep="slope > 0.04"),
                           default = label_parsed)
m.alf.summary  %>%
  pivot_longer(cols = p05:steep) %>% 
  filter(type == 'nSim') %>% 
  ggplot(aes(x = value))+
  geom_histogram(bins = 30, fill = 'grey', colour = 'grey50')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(colour = 'grey'),
        strip.background = element_blank())+
  geom_vline(data = data.frame(value = unlist(m.alf.summary[1,3:4]),
                               name = names(m.alf.summary[1,3:4]),
                               type = 'Real'),
             aes(xintercept = value, colour = 'observerd'), lwd = 1)+
  scale_color_manual(values = 'red', name = NULL)+
  facet_wrap(~name, scale = 'free_x',
             labeller = my_labeller)+
  ylab('Frequency')+
  xlab('Number of loci') ->g2;g2


slopeMax <- max(m.alfnpp$slope)
slopeMin <- min(m.alfnpp$slope)
df1 <- m.alf %>% group_by(type, i) %>% 
  filter(pvalue < 0.05) %>% 
  summarise(minp = min(pvalue),
           # maxR = max(R2),
            slopeabs = max(abs(slope))) 

my_labeller <- as_labeller(c(minp="minimum-p", maxR="maximum-R^2", 
                             slopeabs="steepest-slope"),
                           default = label_parsed)
df1 %>% 
  filter(type == 'nSim') %>% 
pivot_longer(cols = minp: slopeabs) %>% 
  #filter(type == 'nSim') %>% 
  mutate(name= factor(name, levels = names(df1[1,3:4]))) %>% 
  ggplot(aes(x = value))+
  geom_histogram(bins = 30, fill = 'grey', colour = 'grey50')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(colour = 'grey'),
        strip.background = element_blank())+
  geom_vline(data = data.frame(value = unlist(df1[1,3:4]),
                               name = names(df1[1,3:4]),
                               type = 'Real'),
             aes(xintercept = value, colour = 'observerd'), lwd = 1)+
  scale_color_manual(values = 'red', name = NULL)+
  facet_wrap(~name, scale = 'free_x', labeller = my_labeller)+
  scale_y_continuous(breaks = seq(0,80,10))+
  ylab('Frequency')+
  xlab('Estimate') ->g3;g3

grid.arrange(g2, g3, heights = c(2.2,2))

# real negative slopes  ----------
phaseCol <- c(low = terrain.colors(10)[7],
              increase = terrain.colors(10)[1],
              decrease = terrain.colors(10)[4])
#alfReal <- em.alf_df(ph_phase = ph_phase, meta)

#alfdf <- filter(alfReal, phaseNo != 'L1')
alfReal
alfdf

alflociL1 <- split(alfReal, alfReal$loci)
system.time(m.alfnppL1 <- lapply(alflociL1, 
                                 em.alf_npp_model, iteration = 'real') %>% 
              do.call('rbind', .))
rownames(m.alfnppL1) <-  NULL
m.alfnppL1$type <- 'Real'
nrow(m.alfnppL1)


locmeta_ph <- ph@other$loc.metrics[,c(1,24,5,6, 25)]
names(locmeta_ph) <- c('AlleleID', 'uid', 'tas_genome', 'tas_position','rdepth')  
locmeta_ph %>% head
locmeta <- locmeta_ph %>% 
  mutate(scaffold = sub('GL', '', tas_genome),
         scaffold = ifelse(scaffold == '', 0, scaffold),
         scaffold = as.numeric(scaffold)) %>% 
  group_by(scaffold) %>% 
  mutate(position = tas_position - min(tas_position))
locmeta %>% 
 # filter(scaffold < 25 & scaffold>0) %>% 
ggplot(aes(scaffold, position))+
  geom_point()+
  theme_classic()
  
# fig 6 ----
outliers_withoutL1 <- m.alfnpp %>% mutate(abslope = abs(slope)) %>% 
  filter(abslope > 0.02, pvalue < 0.05) %>% 
  dplyr::select(loci) %>% unlist

negslopes <- m.alfnppL1 %>%
  mutate(abslope = abs(slope),
         direction = ifelse(slope < 0, 'negative slope', 'positive slope')) %>%
  filter(abslope > 0.02, pvalue < 0.05, loci %in% outliers_withoutL1) %>% 
  left_join(locmeta) #%>% 
  # mutate(Pse_genome = factor(Pse_genome),
  #        mus_genome = factor(mus_genome)) %>% 
  # rename(pse_genome = Pse_genome)
nrow(negslopes)
scaftb <- negslopes$scaffold %>% table 
scaftb
names(scaftb[scaftb > 9])

table(negslopes$loci %in% outliers_withoutL1)

names(negslopes)
negslopes$R2 %>% sort(., decreasing = T)
names(alfdf)
dfplot <- alfReal %>%
  left_join(negslopes) %>% 
  filter(uid %in% negslopes$uid) 
# intLoci <- dfplot$loci[dfplot$phaseNo == 'I1' & dfplot$alf2 > 0.1 & dfplot$alf2 < 0.25]
# 
# dfplot <- dfplot %>% filter(#loci %in% intLoci, 
#                             scaffold %in% c(1,18,20),
#                             phase != 'decrease')
nrow(dfplot)
ggplot(dfplot, aes(x = npp.log, y = alf2, colour = loci))+
  geom_smooth(method = 'lm', se = T) +
  geom_point()+
  facet_wrap(~direction)+
  theme_classic()+
  theme(legend.position = 'none')
## g6 ------
filter(dfplot, phase != 'decrease2') %>% 
  mutate(alf = ifelse(direction == 'negative slope', alf1, alf2),
         phase = factor(phase, levels = c('low',
                                          'increase',
                                          'decrease'))) %>%
ggplot(aes(phaseNo, alf, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line(linewidth = 0.5, alpha = 0.5)+
  #geom_hline(yintercept = c(0,1), colour = 'grey', lty = 3)+
  # geom_boxplot(aes(group = phaseNo, fill = phase), 
  #              colour = 'grey50',
  #              alpha = 0.5, width = 0.1, size = 0.4)+#
  geom_point(size = 1, alpha = 0.5, aes(group = phaseNo))+
 # facet_wrap(~direction, ncol = 1, scale = 'free')+
  scale_colour_manual(values = phaseCol,#rep('grey',3),
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
  xlab('Population phase #') -> g6
ggsave('./figures/fig6_alf_change_real_1.png', g6,
       units = 'cm', height = 8, width = 16, dpi = 300)

## g8 lm ----
dfplot$pse_genome %>% table
filter(dfplot, phase != 'decrease2',
  #     pse_genome != ''
        #pse_genome == 'HiC_scaffold_18',
       # mus_genome %in% 
       #   c('NC_000070.6_C57BL/6J_chromosome_4_GRCm38.p3_C57BL/6J',           
       #     'NC_000071.6_C57BL/6J_chromosome_5_GRCm38.p3_C57BL/6J',         
       #     'NC_000072.6_C57BL/6J_chromosome_6_GRCm38.p3_C57BL/6J')
       ) %>% 
  mutate(alf = ifelse(direction == 'negative slope', alf1, alf2),
         phase = factor(phase, levels = c('low',
                                          'increase',
                                          'decrease'))) %>% 
 # filter(R2 > 0.75) %>% 
 # arrange(pse_genome) %>% 
  mutate(noloci = 1:nrow(.)) %>% 
  ggplot(aes(npp, alf, 
             group = loci))+
 # geom_point(size = 1, pch = 21, colour = 'grey50') +
  geom_jitter(aes(colour = phase),size = 0.5, 
              width = 0.01) +
  #geom_line()+
 # geom_vline(xintercept = 0.001)+
   scale_colour_manual(values = c(phaseCol, Models = 'grey'),
                       name = c('Population phase'),
                       guide = guide_legend(override.aes = list(
                         size = c(2,2,2,1),
                         linetype = c(rep("blank", 3), "solid"),
                         shape = c(rep(16, 3), NA))))+
  # scale_colour_manual(values = virid(5)[1:4],
  #                     labels = c('unassigned', 'scaffold 18', 'scaffold 22',
  #                                'scaffold 23'))+
  
  geom_smooth(method = 'lm', se = F, aes(colour = 'Models'), 
              alpha = 0.5, linewidth = 0.5)+
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
  ylab('Allele Frequency') -> g8; g8

g8g6 <- grid.arrange(g8, g6+theme(legend.position = 'none'), ncol =1)

ggsave('./figures/fig6_alf_change_real_2.png',g8g6,
       units = 'cm', height = 14, width = 12, dpi = 300)
## g7 -----
filter(dfplot, phase != 'decrease') %>% 
  mutate(alf = ifelse(direction == 'negative slope', alf1, alf2),
         phase = factor(phase, levels = c('low',
                                          'increase',
                                          'decrease'))) %>% 
  ggplot(aes(phaseNo, alf, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line(linewidth = 0.5, alpha = 0.5)+
  #geom_hline(yintercept = c(0,1), colour = 'grey', lty = 3)+
  # geom_boxplot(aes(group = phaseNo, fill = phase), 
  #              colour = 'grey50',
  #              alpha = 0.5, width = 0.1, size = 0.4)+#
  # geom_point(size = 1, alpha = 0.5, aes(group = phaseNo))+
  # facet_wrap(~direction, ncol = 1, scale = 'free')+
  scale_colour_manual(values =  c('grey','grey'),
                      labels = c('low to increase', 
                                 'increase to low'),
                      name = 'Transition')+
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
  xlab('Population phase #') -> g7

alfchangeFig <- grid.arrange(g6 +theme(legend.position = 'none',
                       legend.direction="vertical"),
             g7 +theme(legend.position = 'none',
                       legend.direction="vertical",
                       #plot.margin = margin(b = 22, 7,7,7)
                       ),
             widths = c(2,2))

ggsave('./figures/fig6_alf_change_real.png',
       alfchangeFig, units = 'cm', height = 8, width = 18, dpi = 300)

       