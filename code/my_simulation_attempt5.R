# setup --------
source('libraries.R')

source('functions_simulate_geneflow.R')

# what we know: low pop size = 25* with migration rate 0
#               2007 had an increase rate of ~4.17, round to 4*
#               number of grids to take alf from is 20
#               2007 Fst = 0.014 from Ch1 and 0.017* from Ch2
#               generation time every 3months during booms (increase decrease)
#               geneartion time every year during busts (low)
#               migration = 0 during decreases and lows

# step 1: data -----------
ph <- readRDS('pherm_filtered_genotypes_phases.rds')

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

glBaseMAF <- gl.filter.maf(glBase, threshold = 0.25)

gl.report.monomorphs(glBase)

#glBase <- glBase[,1:1000]
popsReal <- seppop(glBase)
popsReal <- list(glBase)
popsReal <- popsReal[rep(1, 24)]
names(popsReal) <- gridsKeep
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

popSim <- mclapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]), mc.cores = 24)

popSim <- lapply(popSim, fxsex)


simStart <- reduce(popSim, em.gl.join)

# compare alf --------------------

alfsim <- gl.alf(simStart)

alfsim$type <- 'sim'
alfreal$type <- 'real'

alfdiff <- alfreal$alf1 - alfsim$alf1 
hist(alfdiff)
max(abs(alfdiff))

rbind(alfreal, alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 30, position = position_dodge())+
  #geom_freqpoly(bins = 30, lwd = 2)+
  theme_classic()

# step 4: initalise simulation ---------
em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
ne <- 100#nInd(popSim$FRN1)
fst <- 0.017
mEqu <- (1/fst-1)/(4*(ne)) # at equilibrium
m <- ceiling(mEqu*100)/100 # round up 

initialise_conditions <- data.frame(gen = 1:4, 
                                    migration = 0, N = ne, 
                                    phase = 'I', offspring = 4,
                                    npop = nPop(simStart)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))
initialise_conditions

initialiseSim <- pblapply(1,function(x) em.simulate(simStart, mat,initialise_conditions, ncores = 24))
initialiseSim <- em.simulate(simStart, mat,initialise_conditions, ncores = 24)
fstinit <- mclapply(initialiseSim, function(x) gl.fst.pop(x[[1]], nboots = 1), mc.cores = 4)




data.frame(fst = sapply(fstinit, mean, na.rm = T),
           gen = 1:length(fstinit)) %>% 
  filter(complete.cases(fst)) %>% 
  mutate(f17 = abs(fst - 0.017),
         min = ifelse(min(f17)==f17, 'Select', ''),
         grp = 'cool') %>% 
  ggplot(aes(gen, fst, colour = min, group = grp))+
  # geom_hline(yintercept = 0.03, colour = 'orange', lty = 2)+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_smooth(method = 'lm')+
  geom_point(size =4)+
  theme_classic()
saveRDS(initialiseSim, './output/intial_sim.rds')
saveRDS(fstinit, './output/intial_sim_fst.rds')

nSel <- 2
fstinit[[nSel]] %>% mean(., na.rm = T)

alfselected <- gl.alf(initialiseSim[[nSel]][[1]])




alfsim$type <- 'sim'
alfreal$type <- 'real'
alfselected$type <- 'selected'
hist(alfreal$alf2 - alfselected$alf2)

bind_rows(alfreal, alfselected) %>%
  #rbind(alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 40, position = position_dodge())+
  #geom_freqpoly(bins = 40, lwd = 2)+
  theme_classic()

initialisedSelected <- initialiseSim[[nSel]][[1]]
saveRDS(initialisedSelected, 'initialised_simulation_gen0.rds')


# Conditions 3 --------------
initialisedSelected <- readRDS('initialised_simulation_gen0.rds')

phaseCol <- c(L = terrain.colors(10)[7],
              I = terrain.colors(10)[1],
              D = terrain.colors(10)[4])
em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)
mPhase <- c(L = em.mRate(mean(c(0.032,0.027,0.02)), 25), # low
            D = em.mRate(mean(c(0.034,0.03,0.012)), 25), # decrease
            I = em.mRate(mean(c(0.01, 0.005, 0.017)), 50))# increase
# # remove "outliers"
# mPhase <- c(L = em.mRate(mean(c(0.032,0.027)), 25), # low
#             D = em.mRate(mean(c(0.034,0.03)), 25), # decrease
#             I = em.mRate(mean(c(0.01,0.017)), 50))# increase
mPhase


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
           phase == 'D' ~ 25,
         ),
         Ne = ifelse(phaseNo == 'I2', 250, Ne),
         N = Ne,
         # N = ifelse(phase == 'D' & duplicated(phaseNo), 100/4, Ne),
         # N = ceiling(N),
         offspring = 4,
         migration = rep(rep(rev(mPhase),3)[-9], phaselength),
         migration = round(migration, 4))

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
  ))

#conditions2$migration <- rep(c(0.15, 0.05, 0.37, 0.5,0.3, 0.37, 0.5, 0.1), phaselength)
conditions2$migration <- rep(c(0.15, 0.09, 0.42, 0.66,0.3, 0.42, 0.66, 0.09), phaselength)
conditions2


## step 6: run simulation ----------
fstinit[[nSel]] %>% mean(., na.rm = T)


sim <- em.simulate(initialisedSelected, mat, conditions2, ncores = 20)



conditions2$nInd <- sapply(sim, nInd)
conditions2$perPop <- round(sapply(sim, nInd)/sapply(sim, nPop))
conditions2$nInd
saveRDS(sim, 'my_beautiful_fourth_sim.rds')
## step 7: fst ---------
# simfst <- list()
# for(i in 1:length(sim)){
#   simfst[[i]] <- gl.fst.pop(sim[[i]], nboots = 1, nclusters = 20)
#   cat(paste0("\033[0;", 35, "m", i,"\033[0m","\n"))
# }


system.time(simfst <- mclapply(sim, function(x) gl.fst.pop(x, nboots = 1), mc.cores = 19))

missing <- which(!sapply(simfst, is.matrix))

simfstextra <- mclapply(sim[rep(missing,3)], gl.fst.pop, nboots = 1, mc.cores = length(missing)*3)

names(simfstextra) <- paste0('miss', rep(letters[1:3], each =length(missing)), '_',missing)
missing
simfstextra2 <- simfstextra[which(sapply(simfstextra, is.matrix))]

saveRDS(simfst, 'fst_backup.rds')

for (i in 1:length(missing)) {
  missingx <- missing[i]
  fstx <- which(as.numeric(str_sub(names(simfstextra2),7,8)) == missingx)[1]
  simfst[[missingx]] <- simfstextra2[[fstx]]
  
}

sapply(simfst, is.matrix)






#6 15
#simfst[[6]] <- gl.fst.pop(sim[[6]], nboots = 1, nclusters = 20)
conditions2$fst <- sapply(simfst, mean, na.rm = T)

simJoin <- reduce(sim, em.gl.join)
simJoin@other$ind.metrics$phaseNo <- factor(simJoin@other$ind.metrics$phaseNo,
                                            levels = phaseNo)
pop(simJoin) <- simJoin@other$ind.metrics$phaseNo

simPhase <- seppop(simJoin)

em.gl_assign_pop <- function (glx, define.pop = "gridId") 
{
  if (length(define.pop) == 1) {
    if (!define.pop %in% names(glx@other$ind.metrics)) {
      stop(paste(define.pop, "not in ind.metrics; names(glx@other$ind.metrics)"))
    }
    pop(glx) <- glx@other$ind.metrics[, define.pop]
  }
  else {
    pop(glx) <- define.pop
  }
  return(glx)
}

simPhase <- lapply(simPhase, em.gl_assign_pop, define.pop = 'pop')

system.time(simfstphase <- mclapply(simPhase, function(x) gl.fst.pop(x, nboots = 1), mc.cores = 8))

missing <- which(!sapply(simfstphase, is.matrix))

simfstextra <- mclapply(simPhase[rep(missing,3)], gl.fst.pop, nboots = 1, mc.cores = length(missing)*3)
simfstextra2 <- simfstextra[which(sapply(simfstextra, is.matrix))]

simfstphase$D1 <- simfstextra2$D1
simfstphase$L1 <- simfstextra2$L1


sapply(simfstphase, is.matrix)




conditions2$fstphase <- rep(sapply(simfstphase, mean, na.rm = T), phaselength)
conditions2$phaseNo <- factor(conditions2$phaseNo, levels = phaseNo)
## step 8: plot --------

ggplot(conditions2, aes(gen, fst, colour = phase, fill = phase))+
  geom_line(aes(group = replace), linewidth = 1)+
  geom_point(size = 3, pch = 21, colour = 'black')+
  #  geom_point(aes(y = fstch2), colour = 'grey60', alpha = 0.5)+
  theme_classic()+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  theme(legend.position = 'none')+
  scale_x_continuous(breaks = seq(1,19,2))+
  scale_fill_manual(values = phaseCol)+
  scale_color_manual(values = phaseCol)+
  ylim(0,0.04)+
  ylab('Fst')+
  xlab('generation') +
  
  ggplot(conditions2, aes(phaseNo, fstphase, colour = phase))+
  geom_hline(yintercept = 0.017, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.03, colour = 'grey', lty = 2)+
  geom_hline(yintercept = 0.005, colour = 'grey', lty = 2)+
  scale_color_manual(values = c(obs_FST = 'grey40', phaseCol),
                     name = 'Population phase',
                     labels = c('Decrease', 'Increase', 'low', 'observed FST'))+
  geom_point(size = 5)+
  geom_point(aes(y = fstch2, colour = 'obs_FST'),size = 4, alpha = 0.5)+
  theme_classic()+
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.8,0.9),
        legend.background = element_rect(colour = 'grey50'),
        legend.key.size = unit(0.25, units = 'cm'),
        legend.box.spacing = unit(0.1, units = 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))+
  ylim(0,0.04)+
  ylab('Fst')+
  xlab('population phase #') 

ggsave('fst_simulation.png', units = 'cm',
       height = 9, width = 16, dpi = 600)

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

conditions2$pops <- 24
conditions2$loci <- '-'
conditions2$loci[1] <- nLoc(glBase)
names(conditions2)
conditions2 %>% 
  dplyr::select(gen, phaseNo, phase, year, yearsSince,
                loci, pops, offspring, N, migration, nInd,
                fst, fstphase) %>% 
write.csv('conditions_fourth_sim.csv', row.names = F)

### checks --------

simD3 <- sim[[19]]
simD3 %>% gl.report.monomorphs()

phD3 <- gl.keep.pop(ph, pop.list = 'D2', as.pop = 'phaseNo')
phD3 <- gl.keep.loc(phD3, loc.list = simD3$loc.names)
phD3 %>% gl.report.monomorphs()

simD3poly <- gl.filter.monomorphs(simD3)
sfsReal <- gl.sfs(glBase, singlepop = TRUE)
sfsSim  <- gl.sfs(simD3poly, singlepop = TRUE)

sfsSim  <- gl.sfs(simStartx, singlepop = TRUE)

data.frame(type = names(sfsReal),
           real = sfsReal) %>%
  full_join(data.frame(type = names(sfsSim),
                       sim = sfsSim)) %>% 
  mutate(real = ifelse(is.na(real), 0, real),
         bin = as.numeric(sub('d', '', type))) %>% 
  ggplot(aes(bin, real))+
  geom_bar(stat = 'identity')+
  geom_bar(aes(y = sim), stat = 'identity', fill = 'orange', alpha = 0.5)+
  theme_classic()

phobs <- gl.keep.loc(ph, loc.list = glBase@loc.names)
pop(phobs) <- phobs@other$ind.metrics$phaseNo 




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