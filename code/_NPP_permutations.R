
# histogram of sd models

em.alf_npp_model <- function(a1, iteration){
  #a1 <- alfdf %>% filter(loci == alfdf$loci[2])
  
  m <- lm(alf2 ~ npp.log, a1)
  m_summary <- summary(m)
  
  
  data.frame(loci = a1$loci[1], uid = a1$uid[1], i = iteration, R2 = m_summary$r.squared,
             intercept = coef(m_summary)[1,1], slope = coef(m_summary)[2,1], 
             df = m_summary$df[2], tvalue =  coef(m_summary)[2,3],
             pvalue =  coef(m_summary)[2,4])
  
}


em.alf_pvalue_n <- function(A1, xnpp) {
  
  nppdf <- data.frame(period = paste('period', 2:9), npp = xnpp)
  
  A1 <- A1 %>% 
    dplyr::select(-npp, -npp.log) %>% 
    left_join(nppdf) %>% 
    mutate(npp.log = log(npp))
  
  alfloci <- split(A1, A1$loci)
  system.time(m.alfnpp <- lapply(alfloci[1], em.alf_npp_model, iteration = 'unknown') %>% 
                do.call('rbind', .))
  rownames(m.alfnpp) <-  NULL
  n <- sum(m.alfnpp$pvalue < 0.05, na.rm = T)
  
}

em.sd_npvalue <- function(pp){
  xbar <- mean(pp[-1])
  s <- sd(pp[-1])
  x <- pp[1]
  
  z <- (x-xbar)/s
  z
}


em.z_sd_npp <- function(alf_n, xnpp){
  system.time(p_alf_sim <- lapply(alf_n,
                                  function(x) em.alf_pvalue_n(x, xnpp)))
  
  pp <- do.call('c', p_alf_sim)
  s <- em.sd_npvalue(pp)
  s
}



source('./code/libraries_dungog.R')
#source('./code/functions_simulate_geneflow.R')
#ph <- readRDS('./output/pherm_filtered_genotypes_phases.rds')
npp <- read.csv('./output/dataframes/npp_means.csv')[-1,]
npp
# npp random -----------


library(gtools)
library(parallel)
library(foreach)

nppAll <- permutations(8, 8, npp$npp , repeats.allowed=F)
dim(nppAll)
head(nppAll)

fx <- function(x,y) x*y
reduce(8:1, fx)

nppcor <- apply(nppAll, MARGIN = 1, function(x) cor(x, npp$npp))
pabs <- which(abs(nppcor)<0.70)

nppAll_r70 <-  nppAll[pabs,]
dim(nppAll_r70)

#nppAll_n <- nppAll[sample(1:nrow(nppAll_r70), 30),]
#nppAll_n <- nppAll[1:10000,]
#nppAll_n <- nppAll[10001:nrow(nppAll),]
nppAll_n <- nppAll

dim(nppAll)

npp.list <- setNames(split(nppAll_n, seq(nrow(nppAll_n))), NULL)


# load real alf -----
system.time(real <- readRDS('/data/scratch/emily/simulations/alf/real_alf.rds') %>% 
              filter(phaseNo != 'L1') %>% 
              mutate(npp = nppTrue,
                     npp.log = log(npp)))

system.time(p <- em.alf_pvalue_n(real, xnpp = npp$npp))


system.time(realsy <- readRDS('/data/scratch/emily/simulations/sy/alf/real_alf.rds') %>% 
              filter(phaseNo != 'L1') %>% 
              mutate(npp = nppTrue,
                     npp.log = log(npp)))

system.time(psy <- em.alf_pvalue_n(realsy, xnpp = npp$npp))

# all npp for real --------


system.time(pvalues_n_sy <- mclapply(npp.list, 
                                     function(x) em.alf_pvalue_n(real, xnpp = x),
                                     mc.cores = 30))

pvalues_n <- readRDS('./output/alf_random_npp/ph_real_40320x_npp_pvalues_n.rds')
pvalues_n_sy <- readRDS('./output/alf_random_npp/ph_real_40320x_npp_pvalues_n_sy.rds')

pp <- do.call('c', pvalues_n)
ppsy <- do.call('c', pvalues_n_sy)

df2 <- data.frame(type = 'sim', n = pp, spp = 'P. hermannsburgensis', real = p) %>% 
  rbind(data.frame(type = 'sim', n = ppsy, spp = 'S. youngsoni', real = psy))
#rbind(data.frame(type = 'real', n = NA))
dim(df2)

## probability -------------
sum(pvalues_n_sy > psy)/length(pvalues_n_sy)
sum(pvalues_n > p)/length(pvalues_n)


## figure ----------
ggplot(df2, aes(x = n))+
  geom_histogram(colour = 'black',fill = 'grey70') + 
  facet_grid(~spp, scale = 'free')+
  theme_bw()+
  theme(strip.background = element_blank(),
        #  strip.placement = 'outside',
        strip.text.x = element_text(face = 'italic', size = 12),
        legend.position.inside = c(0.9,0.85), 
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(0.5, units = 'cm'),
        # strip.text.y = element_text(face = 'bold', size = 10),
        legend.background = element_rect(colour = 'grey'),
        panel.grid = element_line(colour = 'grey97'))+
  geom_vline(aes(xintercept = real) , colour = 'coral2', lwd = 1)+
  xlab('Number of pvalues < 0.05')+
  ylab('Count') -> gRealandNPP;gRealandNPP

ggsave('./figures/fig_real_npvalues_40320npp.png', plot = gRealandNPP, dpi = 300, units = 'cm',
       height = 8, width = 16)

# 

ggplot(df2, aes(x = 'sim', y = n))+
  geom_boxplot(colour = 'black',fill = 'grey70') + 
  theme_classic()+
  geom_hline(yintercept = p, colour = 'coral2', lwd = 2.5)+
  scale_y_continuous(breaks = c(1275, seq(0, 4000,500)))





# load all ----

dirAlf <- '/data/scratch/emily/simulations/alf/'
alfFiles <- list.files(dirAlf, full.names = TRUE)
head(alfFiles)
alf <- mclapply(alfFiles, readRDS, mc.cores = 30)

# SIM random -------------

set.seed(93)
simSample <- sample(2:length(alf), 30)
hist(simSample, breaks = 30)
sub('/data/scratch/emily/simulations/', '', alfFiles[simSample])

alf_n <- c(alf[1], alf[simSample])
length(alf_n)

# NPP random -------
set.seed(24)
nppSample <- sample(2:length(npp.list), 1000)
npp.list_n <- npp.list[nppSample[301:1000]]

colSums(sapply(npp.list_n, function(x)  x==npp$npp))

npp.list_n_real <- c(1, npp.list_n)
npp.list_n_real[[1]] <- npp$npp
length(npp.list_n_real)
# real npp -----------
sd_value <- NA
for (npp_test in 1:length(npp.list_n_real)) {
  xnpp <- npp.list_n_real[[npp_test]]
  system.time(p_alf_sim <- mclapply(alf_n,
                                    function(x) em.alf_pvalue_n(x, xnpp),
                                    mc.cores = 30))
  
  pp <- do.call('c', p_alf_sim)
  s <- em.sd_npvalue(pp)
  sd_value[npp_test] <- s
  
  print(paste0(npp_test, ': ', round(s,4)))
  
  df2 <- data.frame(type = 'sim', n = pp)
  df2$type[1] <- 'real'
  
  g <- ggplot(df2, aes(n, fill = type, colour = type))+
    geom_histogram(bins = 30)+
    theme_classic()+
    ggtitle(paste0('NPP order: ', paste(round(xnpp), collapse = '; ')))+
    scale_colour_manual(values = c('coral2', 'black'))+
    scale_fill_manual(values = c('coral2', 'grey'))
  ggsave(filename = paste0('/data/scratch/emily/simulations/sy/hist_figures/histogram_',
                           npp_test, '.png'), plot = g, height = 4, width = 6)
  
}



## probability -------------
sum(sd_value[-1] > sd_value[1])/length(sd_value)

## figure ---------------
sd_valueph <- readRDS('./output/alf_random_npp/ph_zvalues_sim30_npp30.rds')

df3 <- data.frame(NPP = 'random', z = sd_value, spp = 'S. youngsoni')
df3$NPP[1] <- ' real'

df4<- data.frame(NPP = 'random', z = sd_valueph, spp = 'P. hermannsburgensis')
df4$NPP[1] <- 'real'


rbind(df3, df4) %>% 
  mutate(NPP = factor(NPP, levels = c('real', 'random'))) %>% 
  ggplot(aes(z, fill = NPP, colour = NPP))+
  geom_histogram(bins = 30)+
  facet_grid(~spp)+
  theme_bw()+
  theme(strip.background = element_blank(),
        #  strip.placement = 'outside',
        strip.text.x = element_text(face = 'italic', size = 12),
        legend.position.inside = c(0.9,0.85), 
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(0.5, units = 'cm'),
        # strip.text.y = element_text(face = 'bold', size = 10),
        legend.background = element_rect(colour = 'grey'),
        panel.grid = element_line(colour = 'grey97'))+
  scale_colour_manual(values = c(random = 'black',real = 'coral2'))+
  scale_fill_manual(values = c(random = 'grey',real = 'coral2'))

# 
# same but mclapply not forloop -----------


system.time(ph_z_values700 <- mclapply(npp.list_n,
                                       function(x) em.z_sd_npp(alf_n, xnpp = x),
                                       mc.cores = 30))

xxx <- do.call('c', ph_z_values700)
#saveRDS(xxx, './output/alf_random_npp/ph_zvalues_sim30_npp700.rds')

## dunnart -----------
# load all ----

dirAlf <- '/data/scratch/emily/simulations/sy/alf/'
alfFiles <- list.files(dirAlf, full.names = TRUE)
head(alfFiles)
alf <- mclapply(alfFiles, readRDS, mc.cores = 30)

## randomly select 30 sims

set.seed(93)
simSample <- sample(2:length(alf), 30)
hist(simSample, breaks = 30)
sub('/data/scratch/emily/simulations/', '', alfFiles[simSample])

alf_nSy <- c(alf[1], alf[simSample])
length(alf_nSy)

system.time(sy_z_values700 <- mclapply(npp.list_n,
                                       function(x) em.z_sd_npp(alf_nSy, xnpp = x),
                                       mc.cores = 30))

xxx <- do.call('c', sy_z_values700)
saveRDS(xxx, './output/alf_random_npp/sy_zvalues_sim30_npp700.rds')

## probability -------------
sd_valueph <- readRDS('./output/alf_random_npp/ph_zvalues_sim30_npp300.rds') %>% 
  c(readRDS('./output/alf_random_npp/ph_zvalues_sim30_npp700.rds'))
sd_value <- readRDS('./output/alf_random_npp/sy_zvalues_sim30_npp300.rds') %>% 
  c(readRDS('./output/alf_random_npp/sy_zvalues_sim30_npp700.rds'))
sum(sd_value[-1] > sd_value[1])/length(sd_value)
sum(sd_valueph[-1] > sd_valueph[1])/length(sd_valueph)

## figure ---------------


df3 <- data.frame(NPP = 'random', z = sd_value, spp = 'S. youngsoni')
df3$NPP[1] <- 'real'

df4<- data.frame(NPP = 'random', z = sd_valueph, spp = 'P. hermannsburgensis')
df4$NPP[1] <- 'real'


rbind(df3, df4) %>% 
  mutate(NPP = factor(NPP, levels = c('real', 'random'))) %>% 
  ggplot(aes(z, fill = NPP, colour = NPP))+
  geom_vline(data = rbind(df3[1,], df4[1,]), linewidth = 1, alpha = 0.25,
             aes(xintercept = z, colour = NPP))+
  geom_histogram(binwidth = 0.25, linewidth = 0.2)+
  facet_grid(~spp, scale = 'free_x', space = 'free_x')+
  theme_bw()+
  theme(strip.background = element_blank(),
        #  strip.placement = 'outside',
        strip.text.x = element_text(face = 'italic', size = 11),
        legend.position = c(0.11,0.75), 
        plot.margin = margin(10,10,10,10),
        legend.key.size = unit(0.5, units = 'cm'),
        # strip.text.y = element_text(face = 'bold', size = 10),
        legend.background = element_rect(colour = 'grey'),
        panel.grid = element_line(colour = 'grey97'))+
  scale_colour_manual(values = c(random = 'black',real = 'coral2'))+
  scale_fill_manual(values = c(random = 'grey',real = 'coral2'))+
  xlab('standard deviations away from the mean') -> gSim; gSim

ggsave('./figures/fig_sim_zvalues.png', plot = gSim, dpi = 300, units = 'cm',
       height = 8, width = 16)

# 