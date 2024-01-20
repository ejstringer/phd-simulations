## data creastion -------
source('./code/functions_genetic-data-creation.R')
source('./code/libraries.R')

ph <- em.filtering('ph')
sy <- em.filtering('sy')

ph@other$ind.metrics %>% names()
rain <- em.rain_caps_phase() %>% 
  mutate(trip = factor(trip)) %>% 
  left_join(data.frame(period = paste('period', 1:9),
            phaseNo = paste0(rep(c('L', 'I', 'D'), 3),
                             rep(1:3, each = 3))))

ph@other$ind.metrics <- ph@other$ind.metrics %>% 
  left_join(rain[, c('trip', 'phase', 'period', 'phaseNo')])

sy@other$ind.metrics <- sy@other$ind.metrics %>% 
  left_join(rain[, c('trip', 'phase', 'period', 'phaseNo')])



saveRDS(ph, './output/pherm_filtered_genotypes_phases.rds')
saveRDS(sy, './output/syoung_filtered_genotypes_phases.rds')

saveRDS(rain, './output/phase_df.rds')
