alfReal %>% 
  mutate(fixed = ifelse(alf1 == 1 | alf1 == 0, TRUE, FALSE)) %>% 
  group_by(phaseNo) %>% 
  summarise(fixed = sum(fixed))

ph_phaseMono <- lapply(ph_phase, gl.filter.monomorphs)


sapply(ph_phaseMono, nLoc)/sapply(ph_phase, nLoc)


sim.gl <- lapply(list.files('/data/scratch/emily/simulations/sy/genlights/', full.names = T), readRDS)



sim_phaseMono <- lapply(, function(x) )
