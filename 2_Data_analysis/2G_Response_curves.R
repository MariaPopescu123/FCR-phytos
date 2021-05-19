abundance <- phytos5
sampleScores <- Q$points[,3]
speciesScores <- Q$species[,3]

responseCurveWA <- function(abundance, sampleScores,
                            speciesScores) {
  
  numSpecies <- length(speciesScores)
  PE <- speciesScores
  ET <- vector(length=numSpecies)
  PA <- vector(length=numSpecies)
  
  for (species in 1:numSpecies){
    # get sample scores of all samples containing the taxon
    scoresWithSpecies <-
      sampleScores[abundance[ ,species]>0]
    # calculate ET
    ET[species] <- sd(scoresWithSpecies) 
    # find all samples within 1 ET of PE
    samplesNearOptimum <-
      (abs(sampleScores - PE[species]) <= ET[species])
    # count the number of these samples
    total <- length(samplesNearOptimum
                    [samplesNearOptimum==TRUE])
    # find abundances of taxon within those samples
    abundancesNearOptimum <-
      abundance[samplesNearOptimum, species]
    
    # extract those samples that actually contain the taxon
    present <- length(abundancesNearOptimum
                      [abundancesNearOptimum>0])
    
    # calculate PA, with correction factor from average PA
    # to peak PA. Note that this correction factor will not
    # be accurate for large PA
    PA[species] <- present/total * 100 * 1.168739
  }
  speciesResponse <- data.frame(PE, ET, PA)
  rownames(speciesResponse) <- colnames(abundance)
  speciesResponse
}

rc <- responseCurveWA(abundance, sampleScores, speciesScores)
no_mix <- c("Oocystis","Rhodomonas","Staurastrum","Spondylosium","Monomastix","Staurodesmus","Selenastrum")
genera <- row.names(rc)

rc1 <- rc %>%
  mutate(genus = genera) %>%
  filter(complete.cases(.)) %>%
  mutate(mix = ifelse(genus %in% no_mix, 0, 1))
for (i in 1:length(rc1$PE)){
  rc1$pdf[i] <- max(dnorm(seq(-2.5,2.5,0.1),mean = rc1$PE[i], sd = rc1$ET[i]))
}
rc1 <- arrange(rc1, ET)

rc2 <- rc %>%
  mutate(genus = genera) %>%
  filter(complete.cases(.)) %>%
  mutate(mix = ifelse(genus %in% no_mix, 0, 1))
for (i in 1:length(rc2$PE)){
  rc2$pdf[i] <- max(dnorm(seq(-2.5,2.5,0.1),mean = rc2$PE[i], sd = rc2$ET[i]))
}
rc2 <- arrange(rc2, ET)

rc3 <- rc %>%
  mutate(genus = genera) %>%
  filter(complete.cases(.)) %>%
  mutate(mix = ifelse(genus %in% no_mix, 0, 1))
for (i in 1:length(rc3$PE)){
  rc3$pdf[i] <- max(dnorm(seq(-2.5,2.5,0.1),mean = rc3$PE[i], sd = rc3$ET[i]))
}
rc3 <- arrange(rc3, ET)

plot(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc1$PE[1], sd = rc1$ET[1]), type = "n", main = "Response curves NMDS 1",xlab = "",ylab = "",ylim = c(0,10))
for (i in 1:length(rc1$PE)){
  if(rc1$mix[i] == 1){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc1$PE[i], sd = rc1$ET[i]),col = "black")}
  if(rc1$mix[i] == 0){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc1$PE[i], sd = rc1$ET[i]),col = "red")}
}

plot(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc2$PE[1], sd = rc2$ET[1]), type = "n", main = "Response curves NMDS 2",xlab = "",ylab = "",ylim = c(0,10))
for (i in 1:length(rc2$PE)){
  if(rc2$mix[i] == 1){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc2$PE[i], sd = rc2$ET[i]),col = "black")}
  if(rc2$mix[i] == 0){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc2$PE[i], sd = rc2$ET[i]),col = "red")}
}

plot(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc3$PE[1], sd = rc3$ET[1]), type = "n", main = "Response curves NMDS 3",xlab = "",ylab = "",ylim = c(0,10))
for (i in 1:length(rc3$PE)){
  if(rc3$mix[i] == 1){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc3$PE[i], sd = rc3$ET[i]),col = "black")}
  if(rc3$mix[i] == 0){lines(seq(-2.5,2.5,0.1),dnorm(seq(-2.5,2.5,0.1),mean = rc3$PE[i], sd = rc3$ET[i]),col = "red")}
}
