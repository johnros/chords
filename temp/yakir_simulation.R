rm(list=ls())
rds.sample3 <- read.csv('temp/binary2000.csv')
names(rds.sample3)

sink(file = 'temp/yakir_simulation.output.txt')
for(k in unique(rds.sample3$repeat.)){
  #   k <- 10
  .data <- filter(rds.sample3, repeat.==k)
  rds.object3 <- chords:::rdsObjectConstructor()
  rds.object3$I.t <- seq_len(nrow(.data))
  rds.object3$degree.in <- .data$degree
  rds.object3$degree.out <- rep(0, length.out = nrow(.data))
  rds.object3$rds.sample$NS1 <- .data$degree
  rds.object3$rds.sample$interviewDt <- .data$time
  rds.object3$estimates <- estimate.b.k(rds.object3)
  cat(rep('_',40), "\n")
  cat('repeat.=',k,"\n")
  cat(rep('_',40), "\n")
  cat("N=",sum(rds.object3$estimates$Nk.estimates),"\n")
  cat("N[k]=\n",rds.object3$estimates$Nk.estimates,"\n")
  cat("Imputation report=\n",rds.object3$estimates$convergence,"\n")
  cat("Theta=",try(getTheta(rds.object3)$theta, silent=TRUE),"\n")
}
sink()


