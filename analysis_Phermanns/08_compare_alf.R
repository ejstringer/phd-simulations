
phases <- readRDS('./output/phase_df.rds')
ph <- readRDS('./output/pherm_filtered_genotypes_phases.rds')
#ph@other$ind.metrics$trip <- ymd(ph@other$ind.metrics$trip)
#ph@other$ind.metrics<- ph@other$ind.metrics %>% left_join(phases)

ph@other$ind.metrics$phaseNo <- factor(ph@other$ind.metrics$phaseNo,
                                       levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                       rep(1:3, each = 3)))

npp <- read.csv('./output/dataframes/npp_means.csv')
npp

ph@other$ind.metrics <- ph@other$ind.metrics %>% left_join(npp)
# sample sizes
phase_nGrids<- rowSums(table(ph@other$ind.metrics$phaseNo, ph@other$ind.metrics$gridId) > 0)
phase_nGrids
# grids 

grids_more2n <- names(which(colSums(table(ph@other$ind.metrics$phaseNo, 
                                          ph@other$ind.metrics$gridId)) > 2))
grids_2n <- names(which(colSums(table(ph@other$ind.metrics$phaseNo, 
                                      ph@other$ind.metrics$gridId)) == 2))
gridNames <- c(grids_more2n, grids_2n[3])
gridNames %>% length # the grids to simulate and to be used in the analysis 

# base of simulation ---------

ph2 <- gl.keep.pop(ph, as.pop = 'phaseNo', pop.list = 'I1')

nSample <- table(ph@other$ind.metrics$phaseNo,
                 factor(ph@other$ind.metrics$gridId)) %>% 
  data.frame()

gridsKeep <- nSample$Var2[nSample$Var1 == 'I1' & nSample$Freq>3] %>% unique()
gridsKeep %>% length # 24 alf grids

pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)
glBase <- gl.filter.monomorphs(glBase)

gl.report.monomorphs(glBase)

# filter real data ----------------------
ph 
gridNames %>% length
glBase@loc.names

ph <- gl.keep.loc(ph, loc.list = glBase@loc.names)
ph_real <- gl.keep.pop(ph, pop.list = gridNames, as.pop = 'gridId')
nInd(ph)-nInd(ph_real)

# sep phase --------------

pop(ph_real) <- ph_real@other$ind.metrics$phaseNo
ph_phase <- seppop(ph_real) 
ph_phase <- lapply(ph_phase, em.gl_assign_pop)

# phase summary
meta <- ph@other$ind.metrics %>% 
  group_by(period, phaseNo, phase, npp) %>% 
  summarise(n = n(),
            ngrids = length(unique(gridId)))

# all freq ----------
em.alf_df <- function(ph_phase, meta){
alfphase <- lapply(ph_phase, gl.alf)

for(i in 1:length(alfphase)){
  
  a <- alfphase[[i]]
  alfphase[[i]]$loci <- rownames(a)
  alfphase[[i]]$phaseNo <- names(alfphase)[i]
  
}

alfdf <- do.call('bind_rows', alfphase) %>% 
  left_join(meta) %>% 
  mutate(npp.log = log(npp),
         uid = sub('-', '_', loci),
         uid = sub('-.*', '', uid),
         uid = sub('_', '-', uid))
rownames(alfdf) <- NULL
table(alfdf$phaseNo)
alfdf$phaseNo <- factor(alfdf$phaseNo,
                        levels = levels(ph@other$ind.metrics$phaseNo))

return(alfdf)
}

em.alf_change <- function(alfdf, type = 'real', iteration = 1){

alfdiff <- alfdf %>% dplyr::select(phaseNo, loci, alf2) %>% 
  pivot_wider(names_from = phaseNo, values_from = alf2) %>% 
  dplyr::select(-loci) %>% 
  as.matrix()

col1 <- ncol(alfdiff)
col2 <- col1-1
changedf <- as.data.frame((alfdiff[,1:col2] - alfdiff[,2:col1]))
names(changedf) <- paste(colnames(alfdiff)[1:col2], 
                         colnames(alfdiff)[2:col1], sep = '_')
head(changedf)

changedf_long <- changedf %>% pivot_longer(cols = everything()) %>% 
  mutate(name = factor(name, 
                       levels =paste(colnames(alfdiff)[1:col2], 
                                     colnames(alfdiff)[2:col1], sep = '_')),
         type = type,
         i = iteration) 
changedf_long
return(changedf_long)
}



 # subset -------
samplex <- table(ph_real@other$ind.metrics$phaseNo, 
                 factor(ph_real@other$ind.metrics$gridId))[-1,]
ncol(samplex) 

em.sample_sim <- function(samplex, simxx){
 samplelist <- list()
  
for(i in 1:nrow(samplex)){
  simxx[[i]]@ind.names <- paste0(simxx[[i]]@ind.names, '_g',
                                 simxx[[i]]@other$ind.metrics$gen)
   metaInd <- simxx[[i]]@other$ind.metrics
   nChosenOnes <- c()
   for (j in 1:ncol(samplex)) {
     n <- samplex[i,j]
     
     if(n > 0){
       p <- rownames(samplex)[i]
       g <- colnames(samplex)[j]
       nChosen <- sample(which(metaInd$phaseNo == p & metaInd$pop == g), n)
       nChosenOnes <- c(nChosenOnes, simxx[[i]]@ind.names[nChosen])
     }
   }
   samplelist[[i]] <- nChosenOnes
 }
 samplex[1:2,1:2]
 rowSums(samplex)
 sapply(samplelist, length)
 
 newSimsub <- list()
 for(x in 1:length(simxx)){
   
   newSimsub[[x]] <- simxx[[x]][which(simxx[[x]]@ind.names %in% samplelist[[x]]),]
   
 }
 newSimsub
 sapply(newSimsub, nInd)
 names(newSimsub) <- names(simxx)
 return(newSimsub)
 }


# alf ----------------
alfReal <- em.alf_df(ph_phase = ph_phase, meta)
changeReal <- em.alf_change(alfReal, type = 'Real', iteration = 1)%>% 
  group_by(name, type, i) %>% 
  summarise(meandiff = mean(value),
            vardiff = var(value))



# setup cores -----------------
ncores <- 25
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

changeSim <- foreach(i = 1:50) %dopar% {
  library(dartR)
  library(tidyverse)
  print(paste('sim', i, '....Go!!'))
  fileload <- list.files('/data/scratch/emily/simulations/genlights/',
                         full.names = T)[i]
  simx <- readRDS(fileload)[20:27]
  simalf <- em.alf_df(ph_phase = simx, meta)
  change_list <- em.alf_change(simalf, type = 'Sim', iteration = i)
  change_list %>% 
    group_by(name, type, i) %>% 
    summarise(meandiff = mean(value),
              vardiff = var(value))
}

change_nSim <- foreach(i = 1:50) %dopar% {
  library(dartR)
  library(tidyverse)
  print(paste('sim', i, '....Go!!'))
  fileload <- list.files('/data/scratch/emily/simulations/genlights/',
                         full.names = T)[i]
  simx <- readRDS(fileload)[20:27]
  nSim_list <- list()
 for(j in 1:20) {
    nsim<- em.sample_sim(samplex, simx)
    nSimalf <- em.alf_df(nsim, meta)
    changenSim <- em.alf_change(nSimalf, type = 'nSim', iteration = j/100)
    nSim_list[[j]] <- changenSim %>% 
      group_by(name, type, i) %>% 
      summarise(meandiff = mean(value),
                vardiff = var(value))
  }
  
  to_save_nSim <- do.call('bind_rows', nSim_list)
  to_save_nSim$i <- to_save_nSim$i + i
  to_save_nSim
}

stopCluster(cl)

changeSim2 <- do.call('bind_rows', changeSim)
change_nSim2 <- do.call('bind_rows', change_nSim)
table(change_nSim2$i)

# save data --------

changeReal$species <- 'ph'
change_nSim2$species <- 'ph'
changeSim2$species <- 'ph'

saveRDS(list(real = changeReal, 
             nsim = change_nSim2,
             sim = changeSim2), './output/ph_vis_alf.rds')


#### t.test --------

rtest <- t.test(changeReal$meandiff, mu = 0)
nsimtest <- t.test(change_nSim2$meandiff, mu = 0)
simtest <- t.test(changeSim2$meandiff, mu = 0)
hist(changeSim2$meandiff)
ttest <- data.frame(species = 'ph',
                       data = c('Real', 'n Simulated', 'Simulated'),
                       mean = c(rtest$estimate, nsimtest$estimate, simtest$estimate),
                       df = c(rtest$parameter, nsimtest$parameter, simtest$parameter),
                       t = c(rtest$statistic, nsimtest$statistic, simtest$statistic),
                       pvalue = c(rtest$p.value, nsimtest$p.value, simtest$p.value))
ttest
write.csv(ttest, './output/dataframes/ph_t_tests.csv', row.names = F)
 # summary stats --------
 change_nSim2 %>% 
  bind_rows(changeReal) %>%
  bind_rows(changeSim2) %>%
  mutate(type = factor(type, levels = c('Real', 'nSim', 'Sim')),
         transition = gsub('[[:digit:]]+', '', name)) %>% 
  group_by(type, i) %>% 
  summarise(meandiff = sum(abs(meandiff)),
            vardiff = sum(vardiff)) %>% 
  pivot_longer(cols = meandiff:vardiff) %>% 
  mutate(name = ifelse(name == 'meandiff', 'Mean', 'Variance')) -> df_sum_meanchage
  # group_by(type) %>% 
  # summarise(mean = mean(meandiff),
  #           sd = sd(meandiff),
  #           lower = mean - sd*2,
  #           upper = mean + sd*2) %>% 
df_sum_meanchage %>% 
  ggplot(aes(type, y = value, fill = type, colour = type))+
 # ggplot(aes(type, y = mean, colour = type))+
  facet_wrap(~name, scale = 'free')+
  geom_boxplot(aes(size = type), outlier.size = 1)+
  #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
  #geom_point(size = 3)+
  theme_bw()+
  scale_fill_manual(values = viridis_pal()(10)[c(1,4,7)],
                      labels = c( 'Real',
                                  expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                 'Simulated'),
                      name = 'Data')+
  scale_colour_manual(values = c(viridis_pal()(15)[c(1)], 'black', 'black'),
                      labels = c( 'Real',
                                  expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                  'Simulated'),
                      name = 'Data',)+
  scale_size_manual(values = c(1,0.5, 0.5))+
  guides(colour = 'none', size = 'none') +
  scale_x_discrete(labels=c("Real" = "Real", "nSim" = expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                            "Sim" = "Simulated"))+
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = 'bold', size = 12),
     #   legend.position = c(0.8,0.8),
        legend.position = 'none',
       # legend.background = element_rect(colour = 'grey'),
        legend.text.align = 0,
        axis.title.y = element_text(size = 11)) +
  ylab('Sum of Change in Allele Frequency') -> fig_sumMeans; fig_sumMeans
 
ggsave('./figures/fig4c_alf_meandiff__sum.png',
       fig_sumMeans, units = 'cm', height = 8, width = 16, dpi = 300)

sum(df_sum_meanchage[df_sum_meanchage$type == 'nSim' & df_sum_meanchage$name == 'Mean',]$value > 0.00223)


  fig_mean <- ggplot(change_nSim2, 
          aes(name, meandiff, colour = type, group = i))+
   geom_line(linewidth = 0.3, alpha = 0.2)+
   geom_line(data = changeSim2, linewidth = 0.2)+
   geom_line(data = changeReal[-1,], linewidth = 1)+
   theme_classic()+
    theme(legend.position = 'none',
          axis.title.x = element_blank())+
    scale_colour_manual(values = viridis_pal()(10)[c(4,1,7)],
                        labels = c('sampled simulation',
                                   'real',
                                   'simulated'),
                        name = 'Data')+
  scale_x_discrete(labels = sub('_', ' - ', changeReal$name[-1]))+
    ylab('Mean Change in Allele Frquency')
  
   
  fig_var <- ggplot(change_nSim2,aes(name, vardiff, colour = type, group = i))+
     geom_line(linewidth = 0.3, alpha = 0.2)+
     geom_line(data = changeSim2, linewidth = 0.2)+
     geom_line(data = changeReal[-1,], linewidth = 1)+
     theme_classic()+
    theme(legend.position = c(0.85,0.8),
          legend.text.align = 0,
          legend.background = element_rect(colour = 'grey'))+
     scale_colour_manual(values = viridis_pal()(10)[c(4,1,7)],
                         labels = c(expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                    'Real',
                                    'Simulated'),
                         name = 'Data')+
    scale_x_discrete(labels = sub('_', ' - ', changeReal$name[-1]))+
    labs(x ="Population Phases")+
    ylab('Variation of Change in Allele Frquency')
 fig_var

  fig_meanvar <- grid.arrange(fig_mean, fig_var)
  ggsave('/home/stringer2/simulations/fig4_alf_diff_simulations.png',
         fig_meanvar, units = 'cm', height = 20, width = 16, dpi = 300)
  
  
  df3 %>% group_by(type, name) %>% 
   summarise(mean = mean(value),
             var = var(value)) %>% 
   ggplot(aes(type, mean, colour = name, group = name))+
   geom_line()+
   geom_point()+
   theme_classic()
 ph_real@other$ind.metrics$phaseNo %>% table
 # npp sim alf-----
 
 # separate by phase ----------
 type.lab <- expression(Real = "Real", nSim = {
 } * italic(n) ~ textstyle(group("", Simulated, "")), Sim = "Simulated")
 typeLabs <- c('Real' = 'Real', 'nSim' = 'n Simulated', 'Sim' = 'Simulated')
 change_nSim2 %>% 
   bind_rows(changeReal) %>%
   bind_rows(changeSim2) %>%
   mutate(type = factor(type, levels = c('Real', 'nSim', 'Sim')),
          transition = gsub('[[:digit:]]+', '', name),
          transition = case_when(
            transition == 'I_D' ~ 'Increase - Decrease',
            transition == 'D_L' ~ 'Decrease - Low',
            transition == 'L_I' ~ 'Low - Increase',
          ),
          transition = factor(transition,
                              levels = c('Increase - Decrease',
                                         'Decrease - Low',
                                         'Low - Increase'))) %>% 
   group_by(type, transition, i) %>% 
   summarise(meandiff = mean(abs(meandiff)),
             vardiff = mean(vardiff)) %>% 
   pivot_longer(cols = meandiff:vardiff) %>% 
   mutate(name = ifelse(name == 'meandiff', 'Mean', 'Variance')) -> df_sum_meanchage
 # group_by(type) %>% 
 # summarise(mean = mean(meandiff),
 #           sd = sd(meandiff),
 #           lower = mean - sd*2,
 #           upper = mean + sd*2) %>% 
 df_sum_meanchage %>% 
   filter(name == 'Variance',
          type != 'Sim'
          ) %>% 
   ggplot(aes(transition, y = value, fill = type, colour = type))+
   # ggplot(aes(type, y = mean, colour = type))+
   #facet_wrap(~type, scale = 'free_y', labeller = labeller(type = typeLabs))+
   geom_boxplot(aes(size = type), outlier.size = 1)+
   #geom_errorbar(aes(ymin = lower, ymax = upper), width = 0)+
   #geom_point(size = 3)+
   theme_bw()+
   scale_fill_manual(values = viridis_pal()(10)[c(1,4,7)],
                     labels = c( 'Real',
                                 expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                 'Simulated'),
                     name = 'Data')+
   scale_colour_manual(values = c(viridis_pal()(15)[c(1)], 'black', 'black'),
                       labels = c( 'Real',
                                   expression({}*italic(n)~textstyle(group("", Simulated, ""))),
                                   'Simulated'),
                       name = 'Data',)+
   scale_size_manual(values = c(1,0.5, 0.5))+
   guides(colour = 'none', size = 'none') +
   # scale_x_discrete(labels=c("Real" = "Real", "nSim" = expression({}*italic(n)~textstyle(group("", Simulated, ""))),
   #                           "Sim" = "Simulated"))+
   # scale_x_discrete(labels=c('Increase - Decrease'= 'Increase',
   #                           'Decrease - Low' = '   -   Decrease     -',
   #                           'Low - Increase' = 'Low'))+
   theme(axis.title.x = element_blank(),
          #axis.text.x = element_text(angle = 10, vjust = 0.75),
         # axis.ticks.x = element_blank(),
         panel.grid = element_blank(),
         strip.background = element_blank(),
         strip.text.x = element_text(size = 12),
            legend.position = c(0.9,0.9),
         #legend.position = 'none',
          legend.background = element_rect(colour = 'grey'),
         legend.text.align = 0,
         axis.title.y = element_text(size = 11)) +
   ylab('Mean Variation Change in Allele Frequency') -> fig_sumMeans; fig_sumMeans
 