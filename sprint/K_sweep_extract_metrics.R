library(stm)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) 
    stop("Must pass the directory with the RDS files from runs.", call.=FALSE)

directory = args[1]

ALL_FILES = list.files(directory, '*.rds', full.names=T)
matches = regmatches(ALL_FILES, gregexpr("[[:digit:]]+", ALL_FILES))
nums = as.numeric(unlist(matches))
ALL_FILES <- ALL_FILES[order(nums)]

rel_diff <- function(a, b) { return(mean(abs(a-b))/mean(abs(a))) }

metrics <- ldply(ALL_FILES, 
       function(x) {
         cat(sprintf('Loading %s\n', x))
         r = readRDS(x)
         deter = r[[1]]
         deter.sum = labelTopics(deter)
         rand = r[[2]]
         N = length(rand)
         metrics <- data.frame(K=deter$settings$dim$K, 
                               beta_diff=rep(NA, N), 
                               topics_equal=rep(NA, N))
         for(ii in 1:length(rand))
         {
            rand_i = rand[[ii]]
            rand_i.sum = labelTopics(rand_i)
            metrics$beta_diff[ii] = rel_diff(unlist(deter$beta[[1]]), unlist(rand_i$beta[[1]])) 
            metrics$topics_equal[ii] <- paste(all.equal(deter.sum, rand_i.sum), collapse=', ')
            metrics$num_its_deter[[ii]] <- deter$convergence$its
            metrics$num_its[[ii]] <- rand_i$convergence$its
         }
         
         return(metrics)
       })

write.csv(metrics, sprintf('%s_metrics.csv', gsub('/$', '', directory)), row.names=F)
