library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(plyr)

## function to handle data import and shaping
f1 <- function(file1, plat){
  
  ## add quality levels, seq platform, etc
  ## reshape data
  wrangle_data <- function(dat) {
    # assign alt allele
    colnames(dat) <- c('A', 'T', 'G', 'C', 'N')
    dat$qual <- NA
    dat$ref <- NA
    #assign quality levels
    #is there not an easier way to do this???
    for(i in 1:nrow(dat)){
      if(i < 5){dat[i,6] <- 1}
      else if(i > 4 & i < 9){dat[i,6] <- 2}
      else{dat[i,6] <- 3}
      # assign ref allele
      if(i %% 4 == 1){dat[i,7] <- 'A'}
      else if(i %% 4 == 2){dat[i,7] <- 'T'}
      else if(i %% 4 == 3){dat[i,7] <- 'G'}
      else{dat[i,7] <- 'C'}
    }
    dat <- melt(dat, id=c('ref','qual'))
    dat$pair <- paste(as.character(dat$ref), as.character(dat$variable), sep='-')
    dat <- dat[dat$value > 0,]# & dat$variable != 'N',]  get rid of 0 cells
    return(dat)
  }
  
  ## transform from counts to proportions
  calc_props <- function(dat){
    combos <- unique(dat$pair)
    for (combo in combos){
      dat[dat$pair == combo,]$value <- dat[dat$pair == combo,]$value/sum(dat[dat$pair == combo,]$value)
    }
    return(dat)
  }
  
  calc_props2 <- function(dat){
    dat_collapse <- data.frame(Platform = rep(plat, 3), Qual=c(1, 2, 3), 
                               Count = rep(0,3), Prop=rep(0,3))
    for (i in 1:nrow(dat_collapse)){
      dat_collapse[i,'Count'] <- sum(dat[dat$qual == dat_collapse[i,'Qual'],'value'])
    }
    combos <- unique(dat_collapse$Qual)
    for (i in 1:nrow(dat_collapse)){
      dat_collapse[i,"Prop"] <- dat_collapse[i,"Count"]/sum(dat_collapse[,"Count"]) 
    }
    #for (combo in combos){
    #  dat_collapse[dat_collapse$Qual == combo,]$Prop <- 
    #    dat_collapse[dat_collapse$Qual == combo,]$Count/sum(dat_collapse[dat_collapse$Qual == combo,]$Count)
    #}
    return(dat_collapse)
  }
  
  # load and wrangle data
  errdat <- read.table(file1, sep='\t', header=FALSE)
  errdat2 <- wrangle_data(errdat)
  errdat2$platform <- plat
  
  # deal with aggregate (quality-free) data
  coll <- data.frame(Platform = rep(plat, 16), 
                     Pair = unique(errdat2$pair), Count = rep(0,16), Prop = rep(0,16))
  # would like to do this in vectorized manner instead of loop...
  for(i in 1:nrow(coll)){
    coll[i,"Count"] <- sum(errdat2[errdat2$pair == coll[i,"Pair"] & 
                                  errdat2$platform == coll[i,"Platform"],'value'])
  }
  for(i in 1:nrow(coll)){
    coll[i,"Prop"] <- coll[i,"Count"]/sum(coll[coll$Platform==coll[i,"Platform"],"Count"])
  }
  # deal with quality data
  ###PUT THIS BACK IN LATER
  errdat3 <- calc_props(errdat2)
  #errdat3 <- calc_props2(errdat2)
  return(list(errdat3, coll)) 
}

perbase_calc <- function(dat, tot){
  dat$perbase <- dat$Count/tot
  return(dat)
}

summ_stats <- function(dat){
  dat$mean <- apply(dat[,3:ncol(dat)], 1, mean)
  dat$std <- apply(dat[,3:(ncol(dat) - 1)], 1, sd)
  dat$sem <- dat$std/sqrt(ncol(dat) - 4)
  return(dat)
}


#------------------------------------------------------------------------------------------------------------------
## ch 21 cap per-base (error_finder7.py - counting total bases only at well-covered and non-variant sites)
setwd("C:/Users/cw3026/Desktop/small/nfs/seqscratch11/cw3026/Dan")
df13h<-f1('1113.hi.errs.21.v7.txt', 'HiSeq'); df13h<-perbase_calc(data.frame(df13h[2]),60528910)
df13n<-f1('1113.nova.errs.21.v7.txt', 'NovaSeq'); df13n<-perbase_calc(data.frame(df13n[2]),188004829)
df14h<-f1('1114.hi.errs.21.v7.txt', 'HiSeq'); df14h<-perbase_calc(data.frame(df14h[2]),59905282)
df14n<-f1('1114.nova.errs.21.v7.txt', 'NovaSeq'); df14n<-perbase_calc(data.frame(df14n[2]),186723366)
df15h<-f1('1115.hi.errs.21.v7.txt', 'HiSeq'); df15h<-perbase_calc(data.frame(df15h[2]),62994353)
df15n<-f1('1115.nova.errs.21.v7.txt', 'NovaSeq'); df15n<-perbase_calc(data.frame(df15n[2]),189976239)
df38h<-f1('1138.hi.errs.21.v7.txt', 'HiSeq'); df38h<-perbase_calc(data.frame(df38h[2]),76063651)
df38n<-f1('1138.nova.errs.21.v7.txt', 'NovaSeq'); df38n<-perbase_calc(data.frame(df38n[2]),207359331)
df39h<-f1('1139.hi.errs.21.v7.txt', 'HiSeq'); df39h<-perbase_calc(data.frame(df39h[2]),71402768)
df39n<-f1('1139.nova.errs.21.v7.txt', 'NovaSeq'); df39n<-perbase_calc(data.frame(df39n[2]),222644447)
df40h<-f1('1140.hi.errs.21.v7.txt', 'HiSeq'); df40h<-perbase_calc(data.frame(df40h[2]),71087837)
df40n<-f1('1140.nova.errs.21.v7.txt', 'NovaSeq'); df40n<-perbase_calc(data.frame(df40n[2]),223851203)

# join all per-base data for within-platform props
dfh <- join_all(list(df13h[,c(-3,-4)],df14h[,c(-3,-4)],df15h[,c(-3,-4)],
                      df38h[,c(-3,-4)],df39h[,c(-3,-4)],df40h[,c(-3,-4)]), by=c("Platform", "Pair"))
dfn <- join_all(list(df13n[,c(-3,-4)],df14n[,c(-3,-4)],df15n[,c(-3,-4)],
                     df38n[,c(-3,-4)],df39n[,c(-3,-4)],df40n[,c(-3,-4)]), by=c("Platform", "Pair"))
dfh <- summ_stats(dfh)
dfn <- summ_stats(dfn)
dfavg <- rbind(dfh, dfn)

# plot the per-base by error type (averages)
ggplot(data=dfavg, aes(x=dfavg$Pair, y=dfavg$mean, fill=dfavg$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Proportion of errors (within platform)", x="Error type (ref-alt)",
       title="Distribution of per-base err. rates on HiSeq and NovaSeq (Chr21)") + 
  theme(legend.title=element_blank()) +
  geom_errorbar(aes(ymin=dfavg$mean-dfavg$sem, ymax=dfavg$mean+dfavg$sem),
                width=0.2,position=position_dodge(0.9))

samph <- data.frame(Platform = c(rep('HiSeq 2500', 6)),
                    Samp = c('1113', '1114','1115','1138','1139','1140'),
                    Rate = c(sum(df13h$perbase),sum(df14h$perbase),sum(df15h$perbase),
                             sum(df38h$perbase),sum(df39h$perbase),sum(df40h$perbase)))
sampn <- data.frame(Platform = c(rep('NovaSeq', 6)),
                    Samp = c('1113', '1114','1115','1138','1139','1140'),
                    Rate = c(sum(df13n$perbase),sum(df14n$perbase),sum(df15n$perbase),
                             sum(df38n$perbase),sum(df39n$perbase),sum(df40n$perbase)))
dfsamp <- rbind(samph, sampn)

# plot overall per-base by sample
ggplot(data=dfsamp, aes(x=dfsamp$Samp, y=dfsamp$Rate, fill=dfsamp$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Error Rate (errs/sequenced base)", x="Sample",
       title="Error rate comparison of HiSeq and NovaSeq (Ch21)") + 
  theme(legend.title=element_blank())

# plot overall per-base by sample for just the three "bad" hiseq samples
ggplot(data=dfsamp[dfsamp$Samp=='1113' | dfsamp$Samp=='1114' | dfsamp$Samp=='1115',], 
       aes(x=dfsamp[dfsamp$Samp=='1113' | dfsamp$Samp=='1114' | dfsamp$Samp=='1115',]$Samp, 
          y=dfsamp[dfsamp$Samp=='1113' | dfsamp$Samp=='1114' | dfsamp$Samp=='1115',]$Rate, 
          fill=dfsamp[dfsamp$Samp=='1113' | dfsamp$Samp=='1114' | dfsamp$Samp=='1115',]$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Error rate (per sequenced base)", #x="Sample",
       title="Overall error rate comparison of sequencing platforms") + 
  theme(legend.title=element_blank(),
        axis.title.x=element_blank())

#------------------------------------------------------------------------------------------------------------------
## ch 21 cap per-base (first 100 base only, error_finder_100base.py)
setwd("Y:/Nova_Analysis/Error_Data")
df13h<-f1('1113.hi.errs.21.100base.txt', 'HiSeq'); df13h<-perbase_calc(data.frame(df13h[2]),47224978)
df13n<-f1('1113.nova.errs.21.100base.txt', 'NovaSeq'); df13n<-perbase_calc(data.frame(df13n[2]),124410565)
df14h<-f1('1114.hi.errs.21.100base.txt', 'HiSeq'); df14h<-perbase_calc(data.frame(df14h[2]),46689496)
df14n<-f1('1114.nova.errs.21.100base.txt', 'NovaSeq'); df14n<-perbase_calc(data.frame(df14n[2]),123575086)
df15h<-f1('1115.hi.errs.21.100base.txt', 'HiSeq'); df15h<-perbase_calc(data.frame(df15h[2]),49175961)
df15n<-f1('1115.nova.errs.21.100base.txt', 'NovaSeq'); df15n<-perbase_calc(data.frame(df15n[2]),125679131)
df38h<-f1('1138.hi.errs.21.100base.txt', 'HiSeq'); df38h<-perbase_calc(data.frame(df38h[2]),59737034)
df38n<-f1('1138.nova.errs.21.100base.txt', 'NovaSeq'); df38n<-perbase_calc(data.frame(df38n[2]),137277697)
df39h<-f1('1139.hi.errs.21.100base.txt', 'HiSeq'); df39h<-perbase_calc(data.frame(df39h[2]),56007891)
df39n<-f1('1139.nova.errs.21.100base.txt', 'NovaSeq'); df39n<-perbase_calc(data.frame(df39n[2]),147467486)
df40h<-f1('1140.hi.errs.21.100base.txt', 'HiSeq'); df40h<-perbase_calc(data.frame(df40h[2]),55740836)
df40n<-f1('1140.nova.errs.21.100base.txt', 'NovaSeq'); df40n<-perbase_calc(data.frame(df40n[2]),148284341)

# join all per-base data for within-platform props
dfh <- join_all(list(df13h[,c(-3,-4)],df14h[,c(-3,-4)],df15h[,c(-3,-4)],
                     df38h[,c(-3,-4)],df39h[,c(-3,-4)],df40h[,c(-3,-4)]), by=c("Platform", "Pair"))
dfn <- join_all(list(df13n[,c(-3,-4)],df14n[,c(-3,-4)],df15n[,c(-3,-4)],
                     df38n[,c(-3,-4)],df39n[,c(-3,-4)],df40n[,c(-3,-4)]), by=c("Platform", "Pair"))
dfh <- summ_stats(dfh)
dfn <- summ_stats(dfn)
dfavg <- rbind(dfh, dfn)

# plot the per-base by error type (averages)
ggplot(data=dfavg, aes(x=dfavg$Pair, y=dfavg$mean, fill=dfavg$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Proportion of errors (within platform)", x="Error type (ref-alt)",
       title="Distribution of per-base err. rates on HiSeq and NovaSeq (Chr21, first 100 bases)") + 
  theme(legend.title=element_blank()) +
  geom_errorbar(aes(ymin=dfavg$mean-dfavg$sem, ymax=dfavg$mean+dfavg$sem),
                width=0.2,position=position_dodge(0.9))

samph <- data.frame(Platform = c(rep('HiSeq', 6)),
                    Samp = c('1113', '1114','1115','1138','1139','1140'),
                    Rate = c(sum(df13h$perbase),sum(df14h$perbase),sum(df15h$perbase),
                             sum(df38h$perbase),sum(df39h$perbase),sum(df40h$perbase)))
sampn <- data.frame(Platform = c(rep('NovaSeq', 6)),
                    Samp = c('1113', '1114','1115','1138','1139','1140'),
                    Rate = c(sum(df13n$perbase),sum(df14n$perbase),sum(df15n$perbase),
                             sum(df38n$perbase),sum(df39n$perbase),sum(df40n$perbase)))
dfsamp <- rbind(samph, sampn)

# plot overall per-base by sample
ggplot(data=dfsamp, aes(x=dfsamp$Samp, y=dfsamp$Rate, fill=dfsamp$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Error Rate (errs/sequenced base", x="Sample",
       title="Error rate comparison of HiSeq and NovaSeq (Ch21, first 100 bases)") + 
  theme(legend.title=element_blank())

#------------------------------------------------------------------------------------------------------------------
## ch 21 cap per-base (error_finder7.py - counting total bases only at well-covered and non-variant sites)
# compare different trio (f478) on diff flowcell to the original two trios (no NovaSeq data for these)
setwd("Z:/Nova_Analysis/Error_Data")
df13h<-f1('1113.hi.errs.21.v7.txt', 'HiSeq'); df13h<-perbase_calc(data.frame(df13h[2]),60528910)
df38h<-f1('1138.hi.errs.21.v7.txt', 'HiSeq'); df38h<-perbase_calc(data.frame(df38h[2]),76063651)
df29h<-f1('1329.hi.errs.21.txt', 'HiSeq'); df29h<-perbase_calc(data.frame(df29h[2]),85699881)
df62h<-f1('562.hi.errs.21.txt', 'HiSeq'); df62h<-perbase_calc(data.frame(df62h[2]),63540706)
df75h<-f1('875.hi.errs.21.txt', 'HiSeq'); df75h<-perbase_calc(data.frame(df75h[2]),60883682)
df56h<-f1('1456.hi.errs.21.txt', 'HiSeq'); df56h<-perbase_calc(data.frame(df56h[2]),79486114)
df25h<-f1('1525.hi.errs.21.txt', 'HiSeq'); df25h<-perbase_calc(data.frame(df25h[2]),96035111)
df10h<-f1('1210.hi.errs.21.txt', 'HiSeq'); df10h<-perbase_calc(data.frame(df10h[2]),83619211)
df37h<-f1('1237.hi.errs.21.txt', 'HiSeq'); df37h<-perbase_calc(data.frame(df37h[2]),72201687)
df82h<-f1('1282.hi.errs.21.txt', 'HiSeq'); df82h<-perbase_calc(data.frame(df82h[2]),66239859)

# join all per-base data for within-platform props
dfh <- join_all(list(df13h[,c(-3,-4)],#df14h[,c(-3,-4)],df15h[,c(-3,-4)],
                     df38h[,c(-3,-4)],#df39h[,c(-3,-4)],df40h[,c(-3,-4)],
                     df29h[,c(-3,-4)],#df30h[,c(-3,-4)],df31h[,c(-3,-4)],
                     df62h[,c(-3,-4)],df75h[,c(-3,-4)],df56h[,c(-3,-4)],
                     df25h[,c(-3,-4)],df10h[,c(-3,-4)],
                     df37h[,c(-3,-4)],df82h[,c(-3,-4)]), by=c("Platform", "Pair"))
dfh <- summ_stats(dfh)

#df87n <- f1('IGMcg1011087_v2_ch21cap_errs.txt', 'NovaSeq'); df87n <- perbase_calc(data.frame(df87n[2]),86799903)
df87n <- f1('IGMcg1011087_v2_ch21cap_errs.txt', 'NovaSeq'); df87n <- perbase_calc(data.frame(df87n[2]),85791515)
df17n <- f1('IGMcg1011117_ch21cap_errs.txt', 'NovaSeq'); df17n <- perbase_calc(data.frame(df17n[2]),91457835)
df97n <- f1('IGMcg1011197_ch21cap_errs.txt', 'NovaSeq'); df97n <- perbase_calc(data.frame(df97n[2]),110746760)
df99n <- f1('IGMcg1011199_ch21cap_errs.txt', 'NovaSeq'); df99n <- perbase_calc(data.frame(df99n[2]),92989288)
df76n <- f1('IGMcg1011176_ch21cap_errs.txt', 'NovaSeq'); df76n <- perbase_calc(data.frame(df76n[2]),151557610)
df13n<-f1('1113.nova.errs.21.v7.txt', 'NovaSeq'); df13n<-perbase_calc(data.frame(df13n[2]),188004829)
df38n<-f1('1138.nova.errs.21.v7.txt', 'NovaSeq'); df38n<-perbase_calc(data.frame(df38n[2]),207359331)
df96n <- f1('IGMcg1011396_ch21cap_errs.txt', 'NovaSeq'); df96n <- perbase_calc(data.frame(df96n[2]),134644301)
df36n <- f1('IGMcg1012236_ch21cap_errs.txt', 'NovaSeq'); df36n <- perbase_calc(data.frame(df36n[2]),114095587)
df587n <- f1('IGMcg1011587_ch21cap_errs.txt', 'NovaSeq'); df587n <- perbase_calc(data.frame(df587n[2]),165376121)
df599n <- f1('IGMcg1011599_ch21cap_errs.txt', 'NovaSeq'); df599n <- perbase_calc(data.frame(df599n[2]),100381341)
df81n <- f1('IGMcg1012081_ch21cap_errs.txt', 'NovaSeq'); df81n <- perbase_calc(data.frame(df81n[2]),189995215)
df34n <- f1('IGMcg1011934_ch21cap_errs.txt', 'NovaSeq'); df34n <- perbase_calc(data.frame(df34n[2]),175759185)

dfn <- join_all(list(df87n[,c(-3,-4)],
                     df17n[,c(-3,-4)],
                     df97n[,c(-3,-4)],
                     df99n[,c(-3,-4)],
                     df76n[,c(-3,-4)],
                     df13n[,c(-3,-4)],df38n[,c(-3,-4)],
                     df96n[,c(-3,-4)],
                     df36n[,c(-3,-4)],
                     df587n[,c(-3,-4)],df599n[,c(-3,-4)],
                     df81n[,c(-3,-4)],df34n[,c(-3,-4)]), by=c("Platform", "Pair"))
dfn <- summ_stats(dfn)

df <- rbind(dfh[,c(1,2,13:15)],dfn[,c(1,2,16:18)])

# plot the per-base by error type (averages)
ggplot(data=df, aes(x=df$Pair, y=df$mean, fill=df$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Error rate per sequenced base", x="Error type (ref-alt)",
       title="Distribution of per-base err. rates on HiSeq and NovaSeq (Chr21)") + 
  theme(legend.title=element_blank()) +
  geom_errorbar(aes(ymin=df$mean-df$sem, ymax=df$mean+df$sem),
                width=0.2,position=position_dodge(0.9))

samph <- data.frame(Platform = c(rep('HiSeq', 10),rep('NovaSeq',13)),
                    Samp = c('1113', #'1114','1115',
                             '1138',#'1139','1140',
                             '1329',#'1330','1331',
                             '562','875','1456','1525', '1210', '1237', '1282',
                             '87','17','97','99','76', '13', '38','96','36','587','599','81','34'),
                    Rate = c(sum(df13h$perbase),#sum(df14h$perbase),sum(df15h$perbase),
                             sum(df38h$perbase),#sum(df39h$perbase),sum(df40h$perbase),
                             sum(df29h$perbase),#sum(df30h$perbase),sum(df31h$perbase),
                             sum(df62h$perbase),sum(df75h$perbase),sum(df56h$perbase),
                             sum(df25h$perbase),sum(df10h$perbase),
                             sum(df37h$perbase),sum(df82h$perbase),
                             sum(df87n$perbase),sum(df17n$perbase),
                             sum(df97n$perbase),sum(df99n$perbase),
                             sum(df76n$perbase),
                             sum(df13n$perbase),sum(df38n$perbase),
                             sum(df96n$perbase),sum(df36n$perbase),
                             sum(df587n$perbase),sum(df599n$perbase),
                             sum(df81n$perbase),sum(df34n$perbase)))


# samph$flow <- rep(NA, nrow(samph))
# 
# for (i in 1:nrow(samph)){
#   if (samph[i,'Samp'] == 1113 || samph[i,'Samp'] == 1114 || samph[i,'Samp'] == 1115){
#     samph[i,'flow'] <- 'C9JN5ANXX'
#   }
#   else if (samph[i,'Samp'] == 1138 || samph[i,'Samp'] == 1139 || samph[i,'Samp'] == 1140){
#     samph[i,'flow'] <- 'C9TCVANXX'
#   }
#   else if (samph[i,'Samp'] == 1329 || samph[i,'Samp'] == 1330 || samph[i,'Samp'] == 1331){
#     samph[i,'flow'] <- 'CA30JANXX'
#   }
#   else if (samph[i,'Samp'] == 1210){
#     samph[i,'flow'] <- 'C9WVBANXX'
#   }
#   else if (samph[i,'Samp'] == 1456){
#     samph[i,'flow'] <- 'CA324ANXX'
#   }
#   else if (samph[i,'Samp'] == 1525){
#     samph[i,'flow'] <- 'CA8ULANXX'
#   }
#   else if (samph[i,'Samp'] == 875){
#     samph[i,'flow'] <- 'C9E02ANXX'
#   }
#   else if (samph[i,'Samp'] == 562){
#     samph[i,'flow'] <- 'C8B7RANXX'
#   }
#   else if (samph[i,'Samp'] == 1237){
#     samph[i,'flow'] <- 'C9NVHANXX'
#   }
#   else if (samph[i,'Samp'] == 1282){
#     samph[i,'flow'] <- 'C9NDFANXX'
#   }
# }

# plot overall per-base by sample
ggplot(data=samph, aes(x=samph$Samp, y=samph$Rate, fill=samph$Platform)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(y="Error Rate (errs/sequenced base)", x="Sample",
       title="Error rate comparison of HiSeq and NovaSeq samples") + 
  theme(legend.title=element_blank())





#------------------------------------------------------------------------------------------------------------------
# error by quality bins

df13h<-f1('1113.hi.errs.21.v7.txt', 'HiSeq 2500'); df13h<-data.frame(df13h[1])
df13n<-f1('1113.nova.errs.21.v7.txt', 'NovaSeq'); df13n<-data.frame(df13n[1])
df14h<-f1('1114.hi.errs.21.v7.txt', 'HiSeq 2500'); df14h<-data.frame(df14h[1])
df14n<-f1('1114.nova.errs.21.v7.txt', 'NovaSeq'); df14n<-data.frame(df14n[1])
df15h<-f1('1115.hi.errs.21.v7.txt', 'HiSeq 2500'); df15h<-data.frame(df15h[1])
df15n<-f1('1115.nova.errs.21.v7.txt', 'NovaSeq'); df15n<-data.frame(df15n[1])
df38h<-f1('1138.hi.errs.21.v7.txt', 'HiSeq 2500'); df38h<-data.frame(df38h[1])
df38n<-f1('1138.nova.errs.21.v7.txt', 'NovaSeq'); df38n<-data.frame(df38n[1])
df39h<-f1('1139.hi.errs.21.v7.txt', 'HiSeq 2500'); df39h<-data.frame(df39h[1])
df39n<-f1('1139.nova.errs.21.v7.txt', 'NovaSeq'); df39n<-data.frame(df39n[1])
df40h<-f1('1140.hi.errs.21.v7.txt', 'HiSeq 2500'); df40h<-data.frame(df40h[1])
df40n<-f1('1140.nova.errs.21.v7.txt', 'NovaSeq'); df40n<-data.frame(df40n[1])

df13 <- rbind(df13h, df13n)
df14 <- rbind(df14h, df14n)
df15 <- rbind(df15h, df15n)
df38 <- rbind(df38h, df38n)
df39 <- rbind(df39h, df39n)
df40 <- rbind(df40h, df40n)

# join all quality data for within-error type props by quality
dfq <- join_all(list(df13[,c(-1,-3)],df14[,c(-1,-3)],df15[,c(-1,-3)],
                     df38[,c(-1,-3)],df39[,c(-1,-3)],df40[,c(-1,-3)]), by=c("platform", "pair", "qual"))
dfq <-dfq[c(1,3,4,2,5,6,7,8,9)]
colnames(dfq) <- c("Qual", "Pair", "Platform", "prop13", "prop14", "prop15", "prop38", "prop39", "prop40")
dfq$mean <- apply(dfq[,4:ncol(dfq)], 1, mean)
dfq$std <- apply(dfq[,4:9], 1, sd)
dfq$sem <- dfq$std/sqrt(6)

#ggplot(data=dfq, aes(x=dfq$Qual, y=dfq$mean, fill=dfq$Platform)) + 
#  geom_bar(stat="identity", position="dodge") +
#  labs(y="Average proportion of errors (within platform)", x="Base quality bin",
#       title="Distribution of error rates by quality bin") + 
#  theme(legend.title=element_blank()) +
#  geom_errorbar(aes(ymin=dfq$mean-dfq$sem, ymax=dfq$mean+dfq$sem),width=0.3,position=position_dodge(0.9))

# join all non-qual data for within-platform props
# <- join_all(list(df13nq[,-3],df14nq[,-3],df15nq[,-3],
#                       df38nq[,-3],df39nq[,-3],df40nq[,-3]), by=c("Platform", "Pair"))
#colnames(dfnq) <- c("Platform", "Pair", "prop13", "prop14", "prop15", "prop38", "prop39", "prop40")
#dfnq$mean <- apply(dfnq[,3:ncol(dfnq)], 1, mean)
#dfnq$std <- apply(dfnq[,3:8], 1, sd)
#dfnq$sem <- dfnq$std/sqrt(6)

## plot the quality-less data
#ggplot(data=dfnq, aes(x=dfnq$Pair, y=dfnq$mean, fill=dfnq$Platform)) + 
#  geom_bar(stat="identity", position="dodge") +
#  labs(y="Proportion of errors (within platform)", x="Error type (ref-alt)",
#       title="Distribution of error types on HiSeq and NovaSeq - Chr15") + 
#  theme(legend.title=element_blank()) +
#  geom_errorbar(aes(ymin=dfnq$mean-dfnq$sem, ymax=dfnq$mean+dfnq$sem),width=0.2,position=position_dodge(0.9))


## function for plotting the quality data
plot1 <- function(dat, combo){
  p <- ggplot(data = dat[dat$Pair == combo,], aes(x=as.factor(dat[dat$Pair == combo,]$Qual), 
              y=dat[dat$Pair == combo,]$mean, fill=dat[dat$Pair == combo,]$Platform)) + 
    geom_bar(stat="identity", position="dodge") +
    labs(y="Err. prop. (within platform)", x="Quality Bin", 
         title=combo) + 
    theme(legend.title=element_blank()) +
    geom_errorbar(aes(ymin=dat[dat$Pair == combo,]$mean-dat[dat$Pair == combo,]$sem, 
                      ymax=dat[dat$Pair == combo,]$mean+dat[dat$Pair == combo,]$sem),
                  width=0.1,position=position_dodge(0.9))
  return(p)
}

## function to create multi plot with shared legend 
# https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# plot the quality data
p1 <- plot1(dfq, 'A-T')
p2 <- plot1(dfq, 'A-G')
p3 <- plot1(dfq, 'A-C')
p4 <- plot1(dfq, 'T-A')
p5 <- plot1(dfq, 'T-G')
p6 <- plot1(dfq, 'T-C')
p7 <- plot1(dfq, 'G-A')
p8 <- plot1(dfq, 'G-T')
p9 <- plot1(dfq, 'G-C')
p10 <- plot1(dfq, 'C-A')
p11 <- plot1(dfq, 'C-T')
p12 <- plot1(dfq, 'C-G')
grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
p1


#----------------------------------------------------------------------------------------------------------------
# positional analysis

pos_analyze <- function(fh, fn, ti){
  dfh <- read.table(fh, sep='\t', header=FALSE); colnames(dfh) <- c('pos', 'err', 'tot')
  dfn <- read.table(fn, sep='\t', header=FALSE); colnames(dfn) <- c('pos', 'err', 'tot')
  dfh$rate <- dfh$err/dfh$tot
  dfh$plat <- 'HiSeq'
  dfn$plat <- 'NovaSeq'
  dfn$rate <- dfn$err/dfn$tot
  #plot(dfn[dfn$rate<0.02,]$pos, dfn[dfn$rate<0.02,]$rate, col='blue')
  #points(dfh[dfh$rate<0.02,]$pos, dfh[dfh$rate<0.02,]$rate, col='red')
  #plot(dfn$pos, dfn$rate, col='blue', ylim=c(0, 0.04), xlab='position', ylab='err. rate')
  #points(dfh$pos, dfh$rate, col='red')
  both <- rbind(dfh, dfn)
  p1 <- ggplot(data=both, aes(x=both$pos, y=both$rate)) + 
    geom_point(stat='identity', aes(colour=factor(both$plat))) +
    labs(y="Error Rate", x="Position on read (0-based)",
         title=ti) + 
    theme(legend.title=element_blank()) +
    coord_cartesian(xlim=c(0, 151), ylim=c(0, 0.04))
  #return(list(dfh, dfn))
  return(p1)
}

setwd("Y:/Nova_Analysis/Error_Data/Pos_Data")

pos13 <- pos_analyze('1113.hi.errpos.21.txt', '1113.nova.errpos.21.txt', 'diagseq1113f404')
#hpos13 <- data.frame(pos13[1]); npos13 <- data.frame(pos13[2])
pos14 <- pos_analyze('1114.hi.errpos.21.txt', '1114.nova.errpos.21.txt', 'diagseq1114f404')
#hpos14 <- data.frame(pos14[1]); npos14 <- data.frame(pos14[2])
pos15 <- pos_analyze('1115.hi.errpos.21.txt', '1115.nova.errpos.21.txt', 'diagseq1115f404')
#hpos15 <- data.frame(pos15[1]); npos15 <- data.frame(pos15[2])
pos38 <- pos_analyze('1138.hi.errpos.21.txt', '1138.nova.errpos.21.txt', 'diagseq1138f413')
#hpos38 <- data.frame(pos38[1]); npos38 <- data.frame(pos38[2])
pos39 <- pos_analyze('1139.hi.errpos.21.txt', '1139.nova.errpos.21.txt', 'diagseq1139f413')
#hpos39 <- data.frame(pos39[1]); npos39 <- data.frame(pos39[2])
pos40 <- pos_analyze('1140.hi.errpos.21.txt', '1140.nova.errpos.21.txt', 'diagseq1140f413')
#hpos40 <- data.frame(pos40[1]); npos40 <- data.frame(pos40[2])

grid_arrange_shared_legend(pos13, pos14, pos15, pos38, pos39, pos40)

pos_analyze1 <- function(fh, fn, ti){
  dfh <- read.table(fh, sep='\t', header=FALSE); colnames(dfh) <- c('pos', 'err', 'tot')
  dfn <- read.table(fn, sep='\t', header=FALSE); colnames(dfn) <- c('pos', 'err', 'tot')
  dfh$rate <- dfh$err/dfh$tot
  dfh$plat <- 'HiSeq 2500'
  dfn$plat <- 'NovaSeq'
  dfn$rate <- dfn$err/dfn$tot
  dfh <- dfh[dfh$pos != 125,]
  dfn <- dfn[dfn$pos != 150,]
  #plot(dfn[dfn$rate<0.02,]$pos, dfn[dfn$rate<0.02,]$rate, col='blue')
  #points(dfh[dfh$rate<0.02,]$pos, dfh[dfh$rate<0.02,]$rate, col='red')
  #plot(dfn$pos, dfn$rate, col='blue', ylim=c(0, 0.04), xlab='position', ylab='err. rate')
  #points(dfh$pos, dfh$rate, col='red')
  both <- rbind(dfh, dfn)
  p1 <- ggplot(data=both, aes(x=both$pos, y=both$rate)) + 
    geom_point(stat='identity', aes(colour=factor(both$plat))) +
    labs(y="Error rate (per sequenced base)", x="Position on read (0-based)",
         title="Error rate by position on read") + 
    theme(legend.title=element_blank()) +
    coord_cartesian(xlim=c(0, 150), ylim=c(0, 0.04))
  #return(list(dfh, dfn))
  return(p1)
}
pos13 <- pos_analyze1('1113.hi.errpos.21.txt', '1113.nova.errpos.21.txt', 'diagseq1113f404')
pos13



#hpos13n <- read.table('1113.hi.errpos.21.NoN.txt', sep='\t', header=FALSE) 
# colnames(hpos13n) <- c('pos', 'err', 'tot')
#hpos13n$rate <- hpos13n$err/hpos13n$tot
#plot(hpos13$pos, hpos13$rate)
#points(hpos13n$pos, hpos13n$rate, col='blue')




