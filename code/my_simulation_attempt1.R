source('code/libraries.R')
source('code/functions_simulate_geneflow.R')

# genetics ----------
ph2 <- readRDS('../data/simulation_base_phgenlight_2007sep.rds')
fst2007sep <- read.csv('../boom-bust-genetics/data_processed/pherm_fst_min4.csv') %>% 
  filter(trip == '2007-09-01')

median(fst2007sep$fst)
pop(ph2) <- ph2@other$ind.metrics$gridId
ph2 <- gl.filter.callrate(ph2, threshold = 1)
ph2 <- gl.filter.monomorphs(ph2)
glBase <- gl.keep.pop(ph2, as.pop = "sex", pop.list = c("m", "f"))
table(glBase@other$ind.metrics$sex, exclude = F)

# spatial ----------
grids <- read.csv("../data/em_gridcoortable.csv") 
grids[grids$coor== "site mean",c("lat", "lon")] <- NA
grids <- grids %>% 
  filter(complete.cases(lon))
xyGrids <- rgdal::project(as.matrix(grids[,c("lon", "lat")]),
                          proj = '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

rownames(xyGrids) <- grids$gridId
#xy <- xyGrids[rownames(xyGrids) %in% c("KSE2", "MC4",  "SS3",  "FRN1", "FRS2", "MC2"),]
#xy <- xyGrids[rownames(xyGrids) %in% c("MC2", "MC4"),]

xy <- xyGrids[1:8,]
nrow(xy)

#xy[,1] <- rep(1:3, 2)
#xy[,2] <- rep(1:2, each = 3)
#plot(xy, col = rainbow(6), pch = 16)
lcost <- xy
mat <- as.matrix(dist(lcost))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal


#mat[mat > -100] <- 1/(nrow(mat))
colSums(mat)
mat
## zero (gen1) ----------

sites <- rownames(xy) # populations
p <- lapply(sites, function(x) gen.zero(glBase, x, 4))
zero <- reduce(p, em.gl.join)
table(pop(zero)) 
nPop(zero)

fstList <- list()
indkeep <- lapply(1:30, function(x) sample(1:nInd(zero), 6*9))
lapply(indkeep, function(x) zero[x,])
fst <- pbapply::pblapply(indkeep, function(x) gl.fst.pop(zero[x,])$Fsts)
sapply(fst, function(x) log(mean(x[x>0], na.rm = T))) %>% mean %>% exp
zero@other$ind.metrics
n.zero <- zero[indkeep,]
system.time(fst <- gl.fst.pop(n.zero))

lapply(indkeep, )

fst$Fsts[fst$Fsts>0] %>% mean(.,na.rm = T) %>% log

# my system ---------
phaseNo <- paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3))
100/10
882/231
1076/167
738/353

25*10

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
# sim -----
# migration, N, phase, offspring
conditions <- data.frame(gen = 1:10, 
                         migration = 0.1, N = rep(c(50, 190), c(6,4)), 
                         phase = rep(c('L', 'I'), c(6,4)), offspring = 4,
                         npop = nPop(zero)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))


head(conditions)
mat

em.simulate <- function(zero, mat, conditions){
sim <- list()
gen <- zero#sim[[11]]

system.time({
for(i in 1:nrow(conditions)){
  # setup
  icon <- conditions[i,]
  m.ind <- icon$migration
  max.N <- icon$N 
  pha <- icon$phase
  off <- icon$offspring
  
  # mice dynamics
  gen <- gene.flow(gen, mat, m.ind)
  pops <- seppop(gen)
  genRepro <- lapply(pops, next.gen, n.offspring = off)
  genSurvive <- lapply(genRepro, the.survivors, max.N = max.N)
  gen <- reduce(genSurvive, em.gl.join)
  gen@other$ind.metrics$iteration <- i
  gen@other$ind.metrics$nInd <- nInd(gen)
  gen@other$ind.metrics <- cbind(gen@other$ind.metrics, icon)
  sim[[i]] <- gen
  txt <- paste("  Generation:", i, pha)
  cat(paste0("\033[0;", 36, "m", txt,"\033[0m","\n"))

}
})

return(sim)
}

em.simulate(zero, mat, conditions[1:3,])
# fst ---------
sapply(sim, nInd)/nPop(zero)

sapply(sim, nInd)*0.1
#sim2 <- lapply(sim, function(x)  x[sample(1:nInd(x), ceiling(nInd(x)*0.2))])
#system.time(fst <- lapply(sim2, gl.fst.pop))
system.time(fst <- pblapply(sim, function(x) gl.fst.pop(x, nboots = 1)))

fstAve <- sapply(fst, function(x) mean(x, na.rm = T))
#118 = 2mins for sim2 20% of samples
#142 = 2mins+ for sim  100% of samples
df1 <- conditions %>% 
  mutate(years = gen, 
         rep = 1)
df1$fst <- fstAve


ggplot(df1, aes(years, fst, colour = phase)) +
  geom_hline(yintercept = c(median(fst2007sep$fst), 0.005, 0.032),
             colour = 'grey', lty = 2)+
  facet_grid(~rep, scale = 'free_x', space = 'free_x')+
  theme_classic()+
  geom_smooth(colour = 'grey', fill = 'grey80')+
  geom_point()+
  geom_text(aes(label = gen), nudge_y = 0.01, colour = 'black')
  scale_y_log10(breaks = c(0.003, 0.005, 0.010, 0.018))
sim[[15]] %>% gl.filter.monomorphs()
lm(log(fst) ~ years * rep, data = df1) %>% summary


ggplot()

simgens <- reduce(sim, em.gl.join)
pop(simgens) <- simgens@other$ind.metrics$phaseNo

glPhase <- seppop(simgens)
