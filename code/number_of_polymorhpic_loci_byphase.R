
# real  ----
ph_phaseMono <- lapply(ph_phase, gl.filter.monomorphs)

sapply(ph_phaseMono, nLoc)/sapply(ph_phase, nLoc)
dfloc <- data.frame(total = sapply(ph_phase, nLoc), 
                    real = sapply(ph_phaseMono, nLoc))
# sim -----
ncores <- 20
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

# sims nloc -------
nSim_Loc <- foreach(i = 1:20) %dopar% {
  library(dartR)
  library(tidyverse)
  print(paste('sim', i, '....Go!!'))
  fileload <- list.files('/data/scratch/emily/simulations/genlights',
                         full.names = T)[i]
  ph_sim <- readRDS(fileload)[20:27]
 
  nsim <- em.sample_sim(samplex, ph_sim)
    
  mono <- lapply(nsim, gl.filter.monomorphs)
  
  sim <- data.frame(x = sapply(mono, nLoc))
  colnames(sim) <- paste0('sim', i)
  
  sim
}

stopCluster(cl)
dfsimloc <- do.call('cbind', nSim_Loc) %>% 
  mutate(phaseNo = rownames(.),
         phaseNo = factor(phaseNo, levels = names(ph_phase)))
dflocAll <- dfloc %>% 
  mutate(phaseNo = rownames(.),
         phaseNo = factor(phaseNo, levels = names(ph_phase)),
         species = 'ph') %>% 
  left_join(dfsimloc) %>% 
  filter(complete.cases(sim1)) %>%
  relocate(phaseNo, species) %>% 
  pivot_longer(cols = real:sim20) %>% 
  mutate(data = substr(name, 1,3))



#write.csv(dflocAll, './output/dataframes/SNPs_real_sim_ph.csv', row.names = F)
dflocAllsy <- read.csv('./output/dataframes/SNPs_real_sim_sy.csv')

dflocAll %>% 
  rbind(dflocAllsy) %>% 
  mutate(species = ifelse(species == 'sy', 'S. youngsoni',
                          'P. hermannsburgensis'),
         value = total - value) %>% 
ggplot(aes(phaseNo, value, colour = data))+
  geom_point(size = 2, aes(alpha = data))+
  scale_colour_manual(values = c('coral3', 'grey60'),
                      labels = c('Observed', 'Simulated'))+
  scale_alpha_manual(values = c(0.8, 0.5), guide = 'none')+
  theme_bw()+
  theme(legend.background = element_rect(colour = 'grey'),
        panel.grid = element_blank(),
        legend.position = c(0.85, 0.9),
        legend.key.size = unit(0.4, units = 'cm'),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'italic'))+
  facet_wrap(~species, ncol = 1, scale = 'free')+
  xlab('Population phase #') +
  ylab('Monomorphic loci') -> fignloc;fignloc

ggsave('./figures/figx_nloc.png',fignloc,
       units = 'cm', height = 14, width = 12, dpi = 300)

## the ---------
dir.create('/data/scratch/emily/simulations/models')
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

# sims m -------
models_nSim <- foreach(i = 1:20) %dopar% {
  library(dartR)
  library(tidyverse)
  print(paste('sim', i, '....Go!!'))
  fileload <- list.files('/data/scratch/emily/simulations/genlights',
                         full.names = T)[i]
  simx <- readRDS(fileload)[20:27]
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