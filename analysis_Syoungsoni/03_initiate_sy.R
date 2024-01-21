em.mRate(0.003, 200)
#mean(gl.fst.pop(simStart, nboots = 1))
# conditons ----------
initialise_conditions <- data.frame(gen = 0, 
                                    phaseNo = 'Initial',
                                    N = 100,
                                    migration = 0.93,#seq(0.90, 0.95, 0.01),
                                    nInd = 4800,
                                    offspring = 6,
                                    pops = nPop(simStart),
                                    loci = nLoc(simStart))

initialise_conditions 

# run sim (generation 0) -----------
initialiseSim  <- pblapply(1:10, function(x) em.simulate(simStart, 
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
    #  fst <- dartR::gl.fst.pop(simx[[k]], nboots = 1)
     fst <- dartR::gl.fst.pop(simx[[k]][[1]], nboots = 1)
      if(is.matrix(fst)) x <- FALSE
    }
    fst
  }
}) # 1.5 minutes
stopCluster(cl)


init_df <- data.frame(fst = sapply(simFst_initiate, mean, na.rm = T),
           gen = 1:length(simFst_initiate),
           m = initialise_conditions$migration) %>% 
  filter(complete.cases(fst)) %>% 
  mutate(f17 = abs(fst - 0.003),
         min = ifelse(min(f17)==f17, 'Select', ''),
         grp = 'cool') 
  ggplot(init_df, aes(gen, fst, colour = m))+
  geom_hline(yintercept = 0.003, colour = 'grey', lty = 2)+
  geom_point(aes(shape = min),size =4)+
  theme_classic()

nSel <- which(init_df$min == 'Select')
initialisedSelected <- initialiseSim[[nSel]][[1]]
simFst_initiate[[nSel]] %>% mean(., na.rm = T)

alfselected <- gl.alf(initialisedSelected)


alfsim$type <- 'sim'
alfreal$type <- 'real'
alfselected$type <- 'selected'
hist(alfreal$alf2 - alfselected$alf2)
max(abs(alfreal$alf2 - alfselected$alf2))
bind_rows(alfreal, alfselected) %>%
  #rbind(alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 30, position = position_dodge())+
  #geom_freqpoly(bins = 40, lwd = 2)+
  theme_classic()

initialisedSelected <- initialiseSim[[nSel]][[1]]
saveRDS(initialisedSelected, 
        './output/initialised_simulation_gen0_sy.rds')
write.csv(initialise_conditions,
          './output/dataframes/initialise_conditions_sy.csv', row.names = F)
