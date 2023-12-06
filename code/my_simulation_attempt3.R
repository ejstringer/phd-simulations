
# setup --------
source('code/libraries.R')
source('../phd-analysis2/code/functions_genetic-data-creation.R')
source('code/functions_simulate_geneflow.R')

# what we know: low pop size = 25* with migration rate 0
#               2007 had an increase rate of ~4.17, round to 4*
#               number of grids to take alf from is 20
#               2007 Fst = 0.014 from Ch1 and 0.017* from Ch2
#               generation time every 3months during booms (increase decrease)
#               geneartion time every year during busts (low)
#               migration = 0 during decreases and lows

# step 1: data -----------
ph <- readRDS('../phd-analysis2/output/Pseudomys_hermannsburgensis_filtered_genotypes.rds')

rain <- em.rain_caps_phase() %>% 
  mutate(trip = factor(trip))
periodID <- data.frame(phaseNo = paste0(rep(c('L', 'I', 'D'), 3),
                                        rep(1:3, each = 3)),
                       period = paste0('period ', 1:9))


ph@other$ind.metrics <- ph@other$ind.metric %>% 
  left_join(rain[, c('trip', 'phase', 'period')]) %>% 
  left_join(periodID)

ph2 <- gl.keep.pop(ph, as.pop = 'phaseNo', pop.list = 'I1')

nSample <- table(ph@other$ind.metrics$phaseNo,
              factor(ph@other$ind.metrics$gridId)) %>% 
  data.frame()
gridsKeep <- nSample$Var2[nSample$Var1 == 'I1' & nSample$Freq>3] %>% unique()
gridsKeep %>% length

gridsxphase <- nSample[nSample$Freq>3,]$Var1 %>% table

pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)
glBaseMaf <- gl.filter.maf(glBase, threshold = 0.1)

# system.time(fstrealmaf <- gl.fst.pop(glBaseMaf, nboots = 1))
# mean(fstrealmaf, na.rm = T)
# median(fstrealmaf, na.rm = T)

popsReal <- seppop(glBase)

# step 2: define site gene flow probabilities ----------
xy <- matrix(nrow = length(gridsKeep),
             ncol = 2,
             data = 1:(length(gridsKeep)*2))
rownames(xy) <- gridsKeep

mat <- as.matrix(dist(xy))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal

colSums(mat)
nrow(mat)
mat[1:5, 1:5]

# step 3: simulate alf ---------
alfreal <- gl.alf(glBase)
hist(alfreal$alf2)

alfreal$minAlf <- ifelse(alfreal$alf1 < alfreal$alf2,
                         alfreal$alf1,
                         alfreal$alf2)
hist(alfreal$minAlf, breaks = 20)

popSim <- pblapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]))

popSim <- lapply(popSim, fxsex)


simStart <- reduce(popSim, em.gl.join)

# system.time(fststart <- gl.fst.pop(simStart, nboots = 1))
# mean(fststart, na.rm = T)
# median(fststart, na.rm = T)

alfsim <- gl.alf(simStart)
alfsim$minAlf <- ifelse(alfsim$alf1 < alfsim$alf2,
                         alfsim$alf1,
                         alfsim$alf2)

hist(alfreal$minAlf, breaks = 20)
alfsim$type <- 'sim'
alfreal$type <- 'real'

hist(alfreal$minAlf - alfsim$minAlf)
sd(alfreal$minAlf - alfsim$minAlf)
max(alfreal$minAlf - alfsim$minAlf) - min(alfreal$minAlf - alfsim$minAlf)



rbind(alfreal, alfsim) %>% 
  ggplot(aes(x = minAlf, fill = type))+
  geom_histogram(bins = 50, position = position_dodge())+
  #geom_freqpoly(bins = 30, lwd = 2)+
  theme_classic()

# step 4: initalise simulation ---------

ne <- nInd(popSim$FRN1)
fst <- 0.017
mEqu <- (1/fst-1)/(4*(ne)) # at equilibrium
m <- ceiling(mEqu*100)/100 # round up 

initialise_conditions <- data.frame(gen = 1, 
                                    migration = 0.64, N = ne, 
                                    phase = 'I', offspring = 4,
                                    npop = nPop(simStart)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))

initialise_conditions$migration[1] <- 0.64

initialiseSim <- lapply(1:4,function(x) em.simulate(simStart, mat,initialise_conditions))

sapply(initialiseSim, nPop)
sapply(initialiseSim, nInd)

fstinit <- pblapply(initialiseSim[1:3], function(x) gl.fst.pop(x[[1]], nboots = 1))

data.frame(fst = sapply(fstinit, mean, na.rm = T),
           gen = 1:length(fstinit)) %>% 
  mutate(min = ifelse(min(fst)==fst, 'Select', ''),
         grp = 'cool') %>% 
  ggplot(aes(gen, fst, colour = min, group = grp))+
 # geom_hline(yintercept = 0.03, colour = 'orange', lty = 2)+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_smooth(method = 'lm')+
  geom_point(size =4)+
  theme_classic()
saveRDS(initialiseSim, './output/intial_sim.rds')
saveRDS(fstinit, './output/intial_sim_fst.rds')
fstinit[[6]] %>% mean(., na.rm = T)
alfselected <- gl.alf(initialiseSim[[4]][[1]])

alfselected$minAlf <- ifelse(alfselected$alf1 < alfselected$alf2,
                        alfselected$alf1,
                        alfselected$alf2)

hist(alfreal$minAlf - alfselected$minAlf)
sd(alfreal$minAlf - alfselected$minAlf)
max(alfreal$minAlf - alfselected$minAlf) - min(alfreal$minAlf - alfselected$minAlf)


hist(alfreal$minAlf, breaks = 20)
alfsim$type <- 'sim'
alfreal$type <- 'real'
alfselected$type <- 'selected'
hist(alfreal$minAlf - alfselected$minAlf)

rbind(alfreal, alfselected) %>%
  #rbind(alfsim) %>% 
  ggplot(aes(x = minAlf, fill = type))+
  geom_histogram(bins = 40, position = position_dodge())+
  geom_freqpoly(bins = 40, lwd = 2)+
  theme_classic()

# step 5: define simulation conditions ---------
phaseNo <- paste0(rep(c('I', 'D', 'L'), 3), rep(1:3, each = 3))[-9]

phaselength <- c(2, 2, 2, 4, 2, 4, 2, 2)

conditions <- data.frame(gen = 1:sum(phaselength), 
                         period = rep(1:8, phaselength),
                         phaseNo = rep(phaseNo, phaselength), 
                         grids = rep(c(23, 9, 3, 20, 9, 12, 16,15), phaselength), 
                         replace = T,
                         year = 0.5) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         rep = str_sub(phaseNo, 2,2),
         migration = ifelse(phase == 'I', 0.15, 0),
         N = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 100,
           phase == 'D' ~ 100/2,
         ),
         N = ifelse(phase == 'D' & duplicated(phaseNo), 100/4, N),
         N = ceiling(N),
         offspring = 4,
         N = ifelse(phaseNo == 'I2', 250, N),
         migration = ifelse(phaseNo == 'I1', 0.15, migration))

for(i in 2:nrow(conditions)){
  x <- ifelse(conditions$phase[i] == 'L', 1, 0.5)
  conditions$year[i] <- conditions$year[i-1] + x
  
  
}

conditions


# step 6: run simulation ----------
fstinit[[6]] %>% mean(., na.rm = T)


sim <- em.simulate(initialiseSim[[6]], mat, conditions)
saveRDS(sim, './output/my_beautiful_sim.rds')


sapply(sim, nInd)
sapply(sim, nPop)
round(sapply(sim, nInd)/sapply(sim, nPop))

# step 7: fst ---------

system.time(simfst <- pblapply(sim, gl.fst.pop, nboots = 1))

conditions$fst <- sapply(simfst, mean, na.rm = T)
# step 8: plot --------

ggplot(conditions, aes(year, fst, colour = phase))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)

Inc <- rev(!duplicated(rev(conditions$phaseNo)) & rev(conditions$phase) == 'I')

conditionsM <- conditions %>% 
  rowwise() %>% 
  mutate(yearsSince = case_when(
    rep == '1' ~ year - 0,
    rep == '2' ~ year - 4,
    rep == '3' ~ year - 11,
  ))

m <- lm(log(fst) ~ yearsSince * rep, data = conditionsM[Inc | conditionsM$phase != 'I',])
m %>% summary

conditionsM[Inc | conditionsM$phase != 'I',] %>% 
  ggplot(aes(yearsSince, fst, colour = phase, group = rep))+
  geom_smooth(method = 'lm', colour = 'grey50', se = F) +
  geom_point(aes(size = phase))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_size_manual(values = c(4,8,4))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.05, colour = 'grey', lty = 2)+
  facet_grid(~rep, scales = 'free_x', space = 'free_x')

ggplot(conditions, aes(year, fst, colour = phase))+
  geom_point(size = 4)+
  theme_classic()+
  geom_hline(yintercept = 0.05, colour = 'orange',lwd = 1, lty = 2,
             alpha = 0.5)+
  geom_hline(yintercept = 0.017, colour = 'lightgreen',
             alpha = 0.5, lwd = 1, lty = 2)
