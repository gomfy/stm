library(ggplot2)
# library(doSNOW)
# library(foreach)
# cl<-makeCluster(4) #change to your number of CPU cores
# registerDoSNOW(cl)
# 
# r <- foreach(K=3:30, .packages=c('stm')) %dopar% {
#   deter <- stm(poliblog5k.docs, poliblog5k.voc, K=K, #prevalence=~rating, 
#                data=poliblog5k.meta, max.em.its = 100)
#   
#   rand <- stm(poliblog5k.docs, poliblog5k.voc, K=K, #prevalence=~rating, 
#               data=poliblog5k.meta, seed=13, 
#               max.em.its = 100, control=list(randomize=TRUE))
#   
#   
#   results = list(deter=deter, rand=rand)
# } 
rel_diff <- function(a, b) { return(mean(abs(a-b))/mean(abs(a))) }

N = length(r)
metrics <- data.frame(K=rep(NA, N), beta_diff=rep(NA, N), topics_equal=rep(NA, N))
for(ii in 1:length(r))
{
  deter = r[[ii]][[1]]
  rand = r[[ii]][[2]]
  K = deter$settings$dim$K
  cat('Num Topics: ', K, '\n')
  deter.sum <- summary(deter)
  rand.sum <- summary(rand)
  metrics$K[ii] = K
  metrics$beta_diff[ii] <- rel_diff(unlist(deter$beta[[1]]), unlist(rand$beta[[1]]))
  metrics$topics_equal[ii] <- paste(all.equal(deter.sum, rand.sum), collapse=', ')
  
  #all.equal(deter, rand)
  cat('\t Summary Equal?: ',  all.equal(deter.sum, rand.sum), '\n')
}

metrics$topics_equal_bool <- metrics$topics_equal == 'TRUE'
metrics$num_its_deter <- sapply(r, function(x) { x[[1]]$convergence$its })
metrics$num_its_rand <- sapply(r, function(x) { x[[2]]$convergence$its })


win.graph(); ggplot(subset(metrics, beta_diff < 0.003), aes(x=K, y=beta_diff)) + 
  geom_point(aes(colour=topics_equal_bool)) + geom_line() +
  xlab('Num Topics') + 
  ylab('Beta (Mean Relative Error Between Random and Deterministic)')
