setwd("Y:/Nova_Analysis/PolyG_Data")

pG1113HiR1 <- read.table('1113Hi.r1.dist.txt', sep='\t', header=FALSE)
colnames(pG1113HiR1) <- c('length', 'count')
pG1113HiR1 <- pG1113HiR1[pG1113HiR1$length > 19,]
pG1113HiR1$pdf <- pG1113HiR1$count / sum(pG1113HiR1$count)
#pG1113HiR1$cdf <- lapply(pG1113HiR1$pdf, function(x) x/sum(x))
pG1113HiR1$cdf <- cumsum(pG1113HiR1$pdf)
plot(pG1113HiR1$length, pG1113HiR1$cdf)

pG1113NovaR1 <- read.table('1113Nova.r1.dist.txt', sep='\t', header=FALSE)
colnames(pG1113NovaR1) <- c('length', 'count')
pG1113NovaR1 <- pG1113NovaR1[pG1113NovaR1$length > 19,]
pG1113NovaR1$pdf <- pG1113NovaR1$count / sum(pG1113NovaR1$count)
pG1113NovaR1$cdf <- cumsum(pG1113NovaR1$pdf)
plot(pG1113NovaR1$length, pG1113NovaR1$pdf)
points(pG1113HiR1$length, pG1113HiR1$pdf)
plot(pG1113NovaR1$length, pG1113NovaR1$cdf)
points(pG1113HiR1$length, pG1113HiR1$cdf)

