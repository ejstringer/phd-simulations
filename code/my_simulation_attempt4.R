
# setup --------
source('code/libraries.R')

source('./code/functions_simulate_geneflow.R')

# what we know: low pop size = 25* with migration rate 0
#               2007 had an increase rate of ~4.17, round to 4*
#               number of grids to take alf from is 20
#               2007 Fst = 0.014 from Ch1 and 0.017* from Ch2
#               generation time every 3months during booms (increase decrease)
#               geneartion time every year during busts (low)
#               migration = 0 during decreases and lows

# step 1: data -----------
ph <- readRDS('./output/pherm_filtered_genotypes_phases.rds')

ph2 <- gl.keep.pop(ph, as.pop = 'phaseNo', pop.list = 'I1')

nSample <- table(ph@other$ind.metrics$phaseNo,
              factor(ph@other$ind.metrics$gridId)) %>% 
  data.frame()
gridsKeep <- nSample$Var2[nSample$Var1 == 'I1' & nSample$Freq>3] %>% unique()
gridsKeep %>% length

gridsxphase <- nSample[nSample$Freq>3,]$Var1 %>% table

pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)
glBase <- gl.filter.monomorphs(glBase)

gl.report.monomorphs(glBase)

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

alfreal$minAlf <- ifelse(alfreal$alf1 < alfreal$alf2,
                         alfreal$alf1,
                         alfreal$alf2)

hist(alfreal$minAlf, breaks = 20)

popSim <- pblapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]))

popSim <- lapply(popSim, fxsex)


simStart <- reduce(popSim, em.gl.join)

# compare alf --------------------
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

alfdiff <- alfreal$alf1 - alfsim$alf1 
hist(alfdiff)
max(abs(alfdiff))

rbind(alfreal, alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 40, position = position_dodge())+
  #geom_freqpoly(bins = 30, lwd = 2)+
  theme_classic()

# step 4: initalise simulation ---------
em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
ne <- nInd(popSim$FRN1)
fst <- 0.017
mEqu <- (1/fst-1)/(4*(ne)) # at equilibrium
m <- ceiling(mEqu*100)/100 # round up 

initialise_conditions <- data.frame(gen = 1, 
                                    migration = 0.6, N = ne, 
                                    phase = 'I', offspring = 4,
                                    npop = nPop(simStart)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))
initialise_conditions

initialiseSim <- pblapply(1:10,function(x) em.simulate(simStart, mat,initialise_conditions))

fstinit <- pblapply(initialiseSim, function(x) gl.fst.pop(x[[1]], nboots = 1))

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

nSel <- 4
fstinit[[nSel]] %>% mean(., na.rm = T)
alfselected <- gl.alf(initialiseSim[[nSel]][[1]])

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
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 40, position = position_dodge())+
  #geom_freqpoly(bins = 40, lwd = 2)+
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
fstinit[[nSel]] %>% mean(., na.rm = T)


sim <- em.simulate(initialiseSim[[nSel]][[1]], mat, conditions)
saveRDS(sim, './output/my_beautiful_sim.rds')


sapply(sim, nInd)
sapply(sim, nPop)
round(sapply(sim, nInd)/sapply(sim, nPop))

# step 7: fst ---------

system.time(simfstextra <- mclapply(sim[c(6,6,6,6,6)], 
                                    function(x) gl.fst.pop(x, nboots = 1), mc.cores = 20))
#6 15
simfst[[6]] <- gl.fst.pop(sim[[6]], nboots = 1, nclusters = 20)
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

# Conditions 2 --------------


em.mRate(0.03, 50)
em.mRate(0.034, 50)
sapply(c,)
em.mRate(0.01, 100)
em.mRate(0.005, 250)

em.mRate(0.01, 50)
em.mRate(0.026, 25)

## step 5: define simulation conditions ---------

k <- c(23, 9, 3, 20, 9, 12, 16,15)
k*(k-1)/2

phaseNo <- paste0(rep(c('I', 'D', 'L'), 3), rep(1:3, each = 3))[-9]
fstch2 <- c(0.017, 0.034, 0.027, 0.005, 0.012, 0.02, 0.01, 0.03)


phaselength <- c(1, 2, 2, 4, 2, 4, 2, 2)

conditions2 <- data.frame(gen = 1:sum(phaselength), 
                         period = rep(1:8, phaselength),
                         fstch2 = rep(fstch2, phaselength),
                         phaseNo = rep(phaseNo, phaselength), 
                         grids = rep(c(23, 9, 3, 20, 9, 12, 16,15), phaselength), 
                         replace = T,
                         year = 0.5) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         rep = str_sub(phaseNo, 2,2),
         Ne = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 100,
           phase == 'D' ~ 100/2,
         ),
         Ne = ifelse(phaseNo == 'I2', 250, Ne),
         N = ifelse(phase == 'D' & duplicated(phaseNo), 100/4, Ne),
         N = ceiling(N),
         offspring = 4)

for(i in 2:nrow(conditions2)){
  x <- ifelse(conditions2$phase[i] == 'L', 1, 0.5)
  conditions2$year[i] <- conditions2$year[i-1] + x
}

conditions2 <- conditions2 %>% 
  rowwise() %>% 
  mutate(yearsSince = case_when(
    rep == '1' ~ year - 0,
    rep == '2' ~ year - 3.5,
    rep == '3' ~ year - 10.5,
  ),
  migration = em.mRate(fstch2, Ne))

changeM <- which(!duplicated(conditions2$phaseNo))[-1]
conditions2$migration[changeM] <- conditions2$migration[changeM-1]
conditions2

## step 6: run simulation ----------
fstinit[[nSel]] %>% mean(., na.rm = T)


sim <- em.simulate(initialiseSim[[nSel]][[1]], mat, conditions2)



sapply(sim, nInd)
sapply(sim, nPop)
round(sapply(sim, nInd)/sapply(sim, nPop))
saveRDS(sim, 'my_beautiful_second_sim_ch2.rds')
## step 7: fst ---------
# simfst <- list()
# for(i in 1:length(sim)){
#   simfst[[i]] <- gl.fst.pop(sim[[i]], nboots = 1, nclusters = 20)
#   cat(paste0("\033[0;", 35, "m", i,"\033[0m","\n"))
# }


system.time(simfst <- pblapply(sim, function(x) gl.fst.pop(x, nboots = 1)))
#6 15
#simfst[[6]] <- gl.fst.pop(sim[[6]], nboots = 1, nclusters = 20)
conditions2$fst <- sapply(simfst, mean, na.rm = T)
## step 8: plot --------

ggplot(conditions2, aes(year, fst, colour = phase))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)


Inc <- rev(!duplicated(rev(conditions2$phaseNo)) & rev(conditions2$phase) == 'I')


m <- lm(log(fst) ~ yearsSince * rep, data = conditions2[Inc | conditions2$phase != 'I',])
m %>% summary

conditions2[Inc | conditions2$phase != 'I',] %>% 
  ggplot(aes(yearsSince, fst, colour = phase, group = rep))+
  geom_smooth(method = 'lm', colour = 'grey50', se = F) +
  geom_point(aes(size = phase))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_size_manual(values = c(4,8,4))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.05, colour = 'grey', lty = 2)+
  facet_grid(~rep, scales = 'free_x', space = 'free_x')

ggplot(conditions2, aes(year, fst, colour = phase))+
  geom_point(size = 4)+
  theme_classic()+
  geom_hline(yintercept = 0.05, colour = 'orange',lwd = 1, lty = 2,
             alpha = 0.5)+
  geom_hline(yintercept = 0.017, colour = 'lightgreen',
             alpha = 0.5, lwd = 1, lty = 2)

# Conditions 3 --------------


em.mRate(0.03, 25)
em.mRate(0.034, 25)
sapply(c,)
em.mRate(0.01, 100)
em.mRate(0.005, 250)

em.mRate(0.01, 50)
em.mRate(0.026, 25)

## step 5: define simulation conditions ---------

k <- c(23, 9, 3, 20, 9, 12, 16,15)
k*(k-1)/2

phaseNo <- paste0(rep(c('I', 'D', 'L'), 3), rep(1:3, each = 3))[-9]
fstch2 <- c(0.017, 0.034, 0.027, 0.005, 0.012, 0.02, 0.01, 0.03)


phaselength <- c(1, 2, 2, 4, 2, 4, 2, 2)

conditions2 <- data.frame(gen = 1:sum(phaselength), 
                          period = rep(1:8, phaselength),
                          fstch2 = rep(fstch2, phaselength),
                          phaseNo = rep(phaseNo, phaselength), 
                          grids = rep(c(23, 9, 3, 20, 9, 12, 16,15), phaselength), 
                          replace = T,
                          year = 0.5) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         rep = str_sub(phaseNo, 2,2),
         Ne = case_when(
           phase == 'L' ~ 30,
           phase == 'I' ~ 120,
           phase == 'D' ~ 30,
         ),
         Ne = ifelse(phaseNo == 'I2', 300, Ne),
         N = Ne,
        # N = ifelse(phase == 'D' & duplicated(phaseNo), 100/4, Ne),
        # N = ceiling(N),
         offspring = 4)

for(i in 2:nrow(conditions2)){
  x <- ifelse(conditions2$phase[i] == 'L', 1, 0.5)
  conditions2$year[i] <- conditions2$year[i-1] + x
}

conditions2 <- conditions2 %>% 
  rowwise() %>% 
  mutate(yearsSince = case_when(
    rep == '1' ~ year - 0,
    rep == '2' ~ year - 3.5,
    rep == '3' ~ year - 10.5,
  ),
  migration = em.mRate(fstch2, Ne))

conditions2$migration <- rep(c(0.15, 0.05, 0.37, 0.5,0.3, 0.37, 0.5, 0.1), phaselength)
conditions2


## step 6: run simulation ----------
fstinit[[nSel]] %>% mean(., na.rm = T)


sim <- em.simulate(initialiseSim[[nSel]][[1]], mat, conditions2)



sapply(sim, nInd)
sapply(sim, nPop)
round(sapply(sim, nInd)/sapply(sim, nPop))
#saveRDS(sim, 'my_beautiful_third_sim.rds')
## step 7: fst ---------
# simfst <- list()
# for(i in 1:length(sim)){
#   simfst[[i]] <- gl.fst.pop(sim[[i]], nboots = 1, nclusters = 20)
#   cat(paste0("\033[0;", 35, "m", i,"\033[0m","\n"))
# }


system.time(simfst <- pblapply(sim, function(x) gl.fst.pop(x, nboots = 1)))
#6 15
#simfst[[6]] <- gl.fst.pop(sim[[6]], nboots = 1, nclusters = 20)
conditions2$fst <- sapply(simfst, mean, na.rm = T)

simJoin <- reduce(sim, em.gl.join)
simJoin@other$ind.metrics$phaseNo <- factor(simJoin@other$ind.metrics$phaseNo,
                                            levels = phaseNo)
pop(simJoin) <- simJoin@other$ind.metrics$phaseNo

simPhase <- seppop(simJoin)

simPhase <- lapply(simPhase, em.gl_assign_pop, define.pop = 'pop')

system.time(simfstphase <- pblapply(simPhase, function(x) gl.fst.pop(x, nboots = 1)))
conditions2$fstphase <- rep(sapply(simfstphase, mean, na.rm = T), phaselength)
conditions2$phaseNo <- factor(conditions2$phaseNo, levels = phaseNo)
## step 8: plot --------

ggplot(conditions2, aes(gen, fst, colour = phase))+
  geom_point()+
  geom_point(aes(y = fstch2), colour = 'grey')+
  theme_classic()+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)


Inc <- rev(!duplicated(rev(conditions2$phaseNo)) & rev(conditions2$phase) == 'I')


m <- lm(log(fst) ~ yearsSince * rep, data = conditions2[Inc | conditions2$phase != 'I',])
m %>% summary

conditions2[Inc | conditions2$phase != 'I',] 
  conditions2 %>% 
  ggplot(aes(yearsSince, fst, colour = phase, group = rep))+
  geom_smooth(method = 'lm', colour = 'grey50', se = F) +
  geom_point(aes(size = phase))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_size_manual(values = c(4,8,4))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)+
  facet_grid(~rep, scales = 'free_x', space = 'free_x')

ggplot(conditions2, aes(year, fst, colour = phase))+
  geom_point(size = 4)+
  theme_classic()+
  geom_hline(yintercept = 0.05, colour = 'orange',lwd = 1, lty = 2,
             alpha = 0.5)+
  geom_hline(yintercept = 0.017, colour = 'lightgreen',
             alpha = 0.5, lwd = 1, lty = 2)


# grids ---------

fstCH2 <- read.csv('../phd-analysis2/output/Pseudomys_hermannsburgensis_fst.csv')
ph@other$ind.metrics$phaseNo <- factor(ph@other$ind.metrics$phaseNo,
                                       levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                       rep(1:3, each = 3)))
tb <- table(ph@other$ind.metrics$gridId, 
      ph@other$ind.metrics$phaseNo)

tb[rowSums(tb>3) > 0,] %>% nrow

tb <- table(fstCH2$pairs, fstCH2$period)
tb[tb[,2] == 0,]


median(fst2007sep$fst)