
# conditons ----------
initialise_conditions <- data.frame(gen = 0, 
                                    phaseNo = 'Initial',
                                    N = 100,
                                    migration = 0.63,
                                    nInd = 4800,
                                    offspring = 4,
                                    pops = nPop(simStart),
                                    loci = nLoc(simStart))
initialise_conditions 

# run sim (generation 0) -----------
initialiseSim <- pblapply(1:20, function(x) em.simulate(simStart, 
                                                        mat, 
                                                        initialise_conditions,
                                                        ncores = 24))

# fst ------
simx <- initialiseSim
tofst <- 1:length(simx)

ncores <- length(simx)
cl <- parallel::makeCluster(ncores)
registerDoParallel(cl)

system.time({
  simFst_initiate <- foreach(k = tofst) %dopar% {
    x <- TRUE
    while(x){
      fst <- dartR::gl.fst.pop(simx[[k]][[1]], nboots = 1)
      if(is.matrix(fst)) x <- FALSE
    }
    fst
  }
}) # 10 minutes
stopCluster(cl)


data.frame(fst = sapply(simFst_initiate, mean, na.rm = T),
           gen = 1:length(simFst_initiate)) %>% 
  filter(complete.cases(fst)) %>% 
  mutate(f17 = abs(fst - 0.017),
         min = ifelse(min(f17)==f17, 'Select', ''),
         grp = 'cool') %>% 
  ggplot(aes(gen, fst, colour = min))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_point(size =4)+
  theme_classic()

nSel <- 16
initialisedSelected <- initialiseSim[[nSel]][[1]]
simFst_initiate[[nSel]] %>% mean(., na.rm = T)

alfselected <- gl.alf(initialisedSelected)


alfsim$type <- 'sim'
alfreal$type <- 'real'
alfselected$type <- 'selected'
hist(alfreal$alf2 - alfselected$alf2)

bind_rows(alfreal, alfselected) %>%
  #rbind(alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 30, position = position_dodge())+
  #geom_freqpoly(bins = 40, lwd = 2)+
  theme_classic()

initialisedSelected <- initialiseSim[[nSel]][[1]]
saveRDS(initialisedSelected, 
        './output/initialised_simulation_gen0.rds')
write.csv(initialise_conditions,
          './output/dataframes/initialise_conditions.csv', row.names = F)
