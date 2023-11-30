source('code/libraries.R')
source('code/functions_simulate_geneflow.R')

# zero ------------
ph <- readRDS('../phd-analysis2/output/Pseudomys_hermannsburgensis_filtered_genotypes.rds')

fst2007sep <- read.csv('../boom-bust-genetics/data_processed/pherm_fst_min4.csv') %>% 
  filter(trip == '2007-09-01')
median(fst2007sep$fst)

ph2 <- gl.keep.pop(ph, as.pop = 'trip', pop.list = '2007-09-01')
meta <- ph@other$ind.metrics
nSample <- table(ph@other$ind.metrics$trip,
                 factor(ph@other$ind.metrics$gridId))
gridsKeep <- names(nSample[22,][nSample[22,]>3])
gridsKeep %>% length


pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)

popsReal <- seppop(glBase)


gen0 <- pblapply(popsReal, 
                 function(x) em.gl.sim.ind(x, n = 100, popname = pop(x)[1]))
sapply(gen0, nInd)-100
gen0x <- list()
for (i in 1:20) {
  bb <- sapply(popsReal, nInd)[i]
  bx <- nInd(gen0[[i]])
  gen0x[[i]] <- gen0[[i]][sample(1:bx, bb),] 
  
}
names(gen0x) <- names(gen0)
gen0all <- reduce(gen0x, em.gl.join)
pop(gen0all) %>% table

sapply(popsReal, nInd)

alfrealgrid <- lapply(popsReal, gl.alf)

# zero ------------
fstx <- 0.02
ne <- 250
(1/fstx-1)/(4*(ne))
p <- pblapply(pops, function(x) gen.zero.alf(x, x@pop[1], 500))

zero <- reduce(p, em.gl.join)
table(pop(zero)) 
nPop(zero)

gen0 <- gene.flow(zero, mat, 0.05)
pops <- seppop(gen0)
gen0Repro <- lapply(pops, next.gen, n.offspring = 1)
zeroEqu <- reduce(gen0Repro, em.gl.join)
table(pop(zeroEqu)) 
nPop(zeroEqu)

# fst wrights ----
popsReal <- seppop(glBase)
names(gen0x) <- names(popsReal)
popx <- gen0
popx <- popsReal
he <- sapply(popx, function(x) mean(gl.He(x), na.rm = T))
he

xx <- t(combn(1:20, 2))
FST <- list()
for(i in 1:nrow(xx)){
  a <- xx[i, 1]
  b <- xx[i, 2]
  
  total <- em.gl.join(popx[[a]], popx[[b]])
  Hs <- mean(he[c(a,b)], na.rm = T)
  Ht <- mean(gl.He(total), na.rm = T)
  
  FST[[i]] <- (Ht-Hs)/Ht
  if(i %% 10 == 0) print(i)
}
reduce(FST, c) %>% mean(., na.rm = T)
simfstx
# fst start -------
zero2 <- zero[seq(1,nInd(zero), 4),alfsim$alf2>0.2]
gen0all@ind.names
table(pop(zero2)) 
test@pop %>% table 
test <- gen0all[c(grep('sim7',gen0all@ind.names),grep('sim6',gen0all@ind.names)),]
system.time(fstreal <- gl.fst.pop(glBase, nboots = 1))
system.time(fststart <- gl.fst.pop(test, nboots = 1))

data.frame(fst = colSums(as.matrix(as.dist(fstreal))), pop = colnames(fstreal)) %>% 
  left_join(data.frame(n = sapply(popsReal, nInd), pop = names(popsReal))) %>% 
  ggplot(aes(n, fst))+
  geom_point()

alfsim <- gl.alf(zero)
alfreal <- gl.alf(glBase)

alfsim$alf1
table(rownames(alfsim)==rownames(alfreal))
dfalf <- data.frame(loci = rownames(alfsim),
                    sim = alfsim$alf1, real = alfreal$alf1) %>% 
  mutate(diff = sim-real)  %>% 
  pivot_longer(cols = sim:real) %>% 
  mutate(diff = ifelse(1:(nrow(alfreal)*2) %% 2 == 0, NA, diff))

hist(dfalf$diff)
ggplot(dfalf,aes(x = value, colour = name))+
  geom_freqpoly()+
  theme_classic()

mean(fststart, na.rm = T)
median(fststart, na.rm = T)

mean(fstreal, na.rm = T)
median(fstreal, na.rm = T)
median(fst2007sep$fst)


# mat -------------
m <- round((1/mean(c(0.014))-1)/(4*(100)), 2)

xy <- matrix(nrow = length(gridsKeep),
             ncol = 2,
             data = 1:(length(gridsKeep)*2))
rownames(xy) <- gridsKeep

mat <- as.matrix(dist(xy))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal

colSums(mat)
nrow(mat)
mat

# conditions ------

phaseNo <- paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3))


phaselength <- c(2, 2, 2, 4, 2, 4, 2, 2)

conditions <- data.frame(gen = 1:sum(phaselength), period = rep(1:8, phaselength),
                         phaseNo = rep(phaseNo[-1], phaselength), 
                         grids = sample(4:10, sum(phaselength), replace = T)
) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         rep = str_sub(phaseNo, 2,2),
         grids = ifelse(phase == 'I', grids + 10, grids),
         grids = ifelse(phase == 'D', grids + 3, grids),
         migration = ifelse(phase == 'I', 0.25, 0),
         N = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 250,
           phase == 'D' ~ 250/2,
         ),
         N = ifelse(phase == 'D' & duplicated(phaseNo), 250/4, N),
         N = ceiling(N),
         offspring = 4) %>% 
  mutate(rep = ifelse(phaseNo == 'L3', 2, rep),
         rep = ifelse(phaseNo == 'L2', 1, rep),
         genx = c(1:6, 1:10, 1:4),
         years = c(1,1.25,2,2.2,3,4,
                   1.1,1.2,1.3,1.4, 2, 2.2, 3,4,5,6,
                   1,1.2,2,2.2))

conditions