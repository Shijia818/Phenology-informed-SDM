
#################
library(lme4)
library(lmerTest)
library(dplyr)
library(car)
library(caper)
library(performance)
library(ape)
library(phyr)
library(phytools)
library(picante)
library(stringr)
library(phyr)
library(ROCR)
library(pROC)
library(reshape)
library(ggplot2)
library(sp)
library(maptools)
library(magrittr)
library(ggpubr)
library(landscapetools)
library(landscapemetrics)
library(raster)
library(foreach)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(reshape2)
library(showtext)
library(ggsci)
library(RColorBrewer)
library(corrplot)

#####################################
Growth_form <- read.csv("Growth_form.csv")
data.phenology <- read.csv("data.phenology.csv")
current.climate <- read.csv("current_climate.csv")
future.climate.AC <- read.csv("ACCESS-CM2.csv")
future.climate.CM <- read.csv("CMCC-ESM2.csv")
future.climate.GI <- read.csv("GISS-E2-1-G.csv")
future.climate.HD <- read.csv("HadGEM3.csv")
future.climate.IN <- read.csv("INM-CM4-8.csv")
future.climate.MI <- read.csv("MIROC6.csv")
data.soil <- read.csv("data.soil.csv",row.names=1)
data.altitude <- read.csv("Altitude.csv",row.names=1)
grid40km <-read.csv("grid40km.csv")
biomod.spdis <- read.csv("distribution_data.csv")

##########################
data.phenology$date <- as.Date(paste(data.phenology$year,"-",data.phenology$month,"-",data.phenology$day,sep=""),format="%Y-%m-%d")
data.phenology$julian <- as.numeric(format(data.phenology$date,format="%j"))


###########################
PDflr <- data.phenology[which(data.phenology$pfl >= 0.5),] 
PDfrt <- data.phenology[which(data.phenology$pfr >= 0.5),] 
PDbud <- data.phenology[which(data.phenology$pbd >= 0.5),] 

################################ construct phenology-relationship relationship ###########

d <- mutate_at(PDflr, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))

pheno.model <- lmerTest::lmer(julian ~ bio1 + bio12 + bio4 + bio15 + Growth_form + Status + year + bio1:Growth_form + 
                                bio12:Growth_form + bio4:Growth_form + bio15:Growth_form + bio1:Status + bio4:Status + bio12:Status + bio15:Status + (1|latitude:longitude) +
                                (1 | species) +
                                (0 + bio1 | species) + 
                                (0 + bio12 | species) + 
                                (0 + bio4 | species) + 
                                (0 + bio15 | species), data = d, REML =T,
                              control = lmerControl(optimizer = "bobyqa", 
                                                    optCtrl = list(maxfun = 2e5)))
summary(pheno.model)
model_performance(pheno.model)

##### check residuals of models to test spatial auto-correlations ####

dists <- dist(PDflr[,6:7])
data.dists <- as.matrix(dist(cbind(data$longitude,data$latitude)))
data.dists.inv <- 1/data.dists
diag(data.dists.inv) <- 0
data.dists.inv[!is.finite(data.dists.inv)] = 0
residual <- summary(pheno.model)$residuals
Moran.I(residual,data.dists.inv)

#### extract random components #####################
random_effects <- ranef(pheno.model)$species
random_intercept <-coef(pheno.model)$species


#####  check the models with phylogenetic LMMs that account for phylogenetic relationship ###

tip.label.1 <- str_extract(tree$tip.label,pattern = "_.*_") %>% str_replace_all(pattern = "_",replacement = " ") %>% str_squish() %>% str_replace(pattern = " ",replacement = "_")
spnames <- unique(PDflr$species) 
tree$tip.label <- tip.label.1
postemp <- which(!tree$tip.label %in% spnames)
tree.final <- drop.tip (tree, tree$tip.label[postemp])

d.pm <- mutate_at(PDflr, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))
pm_flr <- pglmm(julian ~ bio1 + bio12 + bio4 + bio15 + Growth_form  + year + bio1:Growth_form + 
                  bio12:Growth_form + bio4:Growth_form + bio15:Growth_form + (1 | latitude:longitude) +
                  (1 | species) +
                  (0 + bio1 | species) + 
                  (0 + bio12 | species) + 
                  (0 + bio4 | species) + 
                  (0 + bio15 | species), data = d.pm,
                family = "Gaussian",REML = TRUE,
                cov_ranef = list(sp = tree.final))


###### predict phenology at grid-cell level for each species ######

current.climate.1 <- current.climate[,c(1,2,5,13,16)]
current.climate.2 <- left_join(grid40km, current.climate.1,by = "Gridcode")

Species <- unique(PDflr$species)
test <- rbind()
for (i in 1:length(Species)){
  test <- rbind(test,current.climate.2)}

Species.rev <- rep(Species,nrow(current.climate)) %>% sort() %>% cbind(test)
colnames(Species.rev) [1] <- "species"
Species.rev <- Species.rev[,-2]

Growth_form <- Growth_form[which(Growth_form$Species %in% Species.rev$species),]
Growth_form.1 <- Growth_form[match(Species.rev$species,Growth_form$Species),]
Species.rev$Growth_form <- Growth_form.1$Growth_Form
Species.rev$Growth_form <- Growth_form.1$Growth_Form
Species.rev$Status <- Growth_form.1$Status
Species.rev$year <- "1970"
Species.rev$year <- as.integer(Species.rev$year)

model.current <- mutate_at(Species.rev, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))
predicted.current <- predict(pheno.model, model.current, type="response",allow.new.levels=T)
model.current$predict <- predicted.current

###### predict phenology at grid-cell level for each species under future conditions ######
### we use one GCM as example ###########

future.climate.1 <- future.climate.IN[,c(1,2,5,13,16)]
future.climate.2 <- left_join(grid40km, future.climate.1,by = "Gridcode")

Species <- unique(data$species)
test <- rbind()
for (i in 1:length(Species)){
  test <- rbind(test,future.climate.2)}

Species.rev <- rep(Species,nrow(future.climate)) %>% sort() %>% cbind(test)
colnames(Species.rev) [1] <- "species"
Species.rev <- Species.rev[,-2]

Growth_form <- Growth_form[which(Growth_form$Species %in% Species.rev$species),]
Growth_form.1 <- Growth_form[match(Species.rev$species,Growth_form$Species),]
Species.rev$Growth_form <- Growth_form.1$Growth_Form
Species.rev$Growth_form <- Growth_form.1$Growth_Form
Species.rev$Status <- Growth_form.1$Status
Species.rev$year <- "2070"
Species.rev$year <- as.integer(Species.rev$year)

model.future <- mutate_at(Species.rev, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))
predicted.future <- predict(pheno.model, model.future, type="response",allow.new.levels=T)
model.future$predict <- predicted.future

#### Using phenology to predict species distributions #### 

###Note that your predictors have very different ranges,so if we want to directly compare the effects or test for interactions ##########
##we need to first rescale them or conduct PCA ###########

data.current.T <- current.climate[,c(1,5,6,8,9,10,11)]
data.current.P <- current.climate[,c(12,13,14,16,17,18,19)]
data.current.V <- current.climate[,c(2,3,4,7,15)]

results.T <- prcomp(data.current.T,scale = TRUE)
results.P <- prcomp(data.current.P,scale = TRUE)
results.V <- prcomp(data.current.V,scale = TRUE)
results.soil <- prcomp(data.soil,scale = TRUE)



data.future.T <- future.climate.IN[,c(1,5,6,8,9,10,11)]
data.future.P <- future.climate.IN[,c(12,13,14,16,17,18,19)]
data.future.V <- future.climate.IN[,c(2,3,4,7,15)]

results.T <- prcomp(data.future.T,scale = TRUE)
results.P <- prcomp(data.future.P,scale = TRUE)
results.V <- prcomp(data.future.V,scale = TRUE)


### make grid-cell climate database ##

grid.climate.current <- as.data.frame(matrix(NA,nrow(current.climate),8))
rownames(grid.climate.current) <- rownames(current.climate)
colnames(grid.climate.current) <- c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","Altitude")
grid.climate.current[,1] <- (-1*results.T$x)[,1]
grid.climate.current[,c(2,3)] <- (-1*results.P$x)[,c(1,2)]
grid.climate.current[,c(4,5)] <- (-1*results.V$x)[,c(1,2)]
grid.climate.current[,c(6,7)] <- (-1*results.soil$x)[,c(1,2)]
grid.climate.current[,8] <- (data.altitude$Altitude - mean(data.altitude$Altitude))/sd(data.altitude$Altitude)

##############################################

grid.climate.future <- as.data.frame(matrix(NA,nrow(future.climate.IN),8))
rownames(grid.climate.future) <- rownames(future.climate.IN)
colnames(grid.climate.future) <- c("PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","Altitude")
grid.climate.future[,1] <- (-1*results.T$x)[,1]
grid.climate.future[,c(2,3)] <- (-1*results.P$x)[,c(1,2)]
grid.climate.future[,c(4,5)] <- (-1*results.V$x)[,c(1,2)]
grid.climate.future[,c(6,7)] <- (-1*results.soil$x)[,c(1,2)]
grid.climate.future[,8] <- (data.altitude$Altitude - mean(data.altitude$Altitude))/sd(data.altitude$Altitude)


###################################

mydata <- melt(biomod.spdis, id.vars = "Gridcode")
colnames(mydata)[2:3] <- c("Species", "Present")
model.data <- merge(mydata, grid.climate.current, by="Gridcode", all.x = TRUE) 
res <- which(model.data$Species %in% model.current$Species)
model.data <- model.data[res,]
model.data.final <- merge(model.data, model.current, by.x = c("Gridcode", "Species"), all.x = TRUE)
colnames(model.data.final)[22] <- 'flr'
full.model <-glmer(as.numeric(Presence) ~ PCA1 + PCA2  + PCA3 + PCA4 + 
                     PCA5 + PCA6 + PCA7 + Altitude + PCA1*flr + PCA2*flr+
                     PCA3*flr + PCA4*flr + PCA5*flr + PCA6*flr + PCA7*flr  +
                     Altitude*flr + (1 + PCA1 + PCA2+ PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA8 + Altitude | Species),
                   data = model.data.final,family = binomial)


###### Models without phenology #####
full.model.c <- glmer(as.factor(Presence) ~ PCA1 + PCA2  + PCA3 + PCA4 + 
                        PCA5 + PCA6 + PCA7 + Altitude + (1 + PCA1 + PCA2+ PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + Altitude | Species),
                      data = model.data.final,family = binomial)

####### predict probability of occurrence of each species at each grid cell ###

predicted <- predict(full.model, model.data.final, type="response")
predicted.current <- model.data.final[,c("Gridcode","Species","Growth_form","Status","flr")]
predicted.current$Predicted <- predicted

############# probability of occurrence in the future ###################

mydata <- melt(biomod.spdis, id.vars = "Gridcode")
colnames(mydata)[2:3] <- c("Species", "Present")
model.data <- merge(mydata, grid.climate.future, by="Gridcode", all.x = TRUE) 
model.data.final <- merge(model.data, model.future, by.x = c("Gridcode", "Species"), all.x = TRUE)
colnames(model.data.final)[22] <- 'flr'
predicted <- predict(full.model, model.data.final, type="response")
predicted.future <- model.data.final[,c("Gridcode","Species","Growth_form","Status","flr")]
predicted.future$Predicted <- predicted

##### calculate threshold ###

spnames <- unique(predicted.current$Species)
data.threshold <- as.data.frame(matrix(NA,length(spnames),2))
colnames(data.threshold) <- c("Species","threshold")
data.threshold[,1] <- spnames
biomod.spdis.1 <- biomod.spdis[,match(spnames,colnames(biomod.spdis))]

for(i in 1:length(unique(spnames))){
  sp <- predicted.current$Species[i]
  dat <- predicted.current[which(predicted.current$Species==sp),]
  prediction <- dat$Predicted
  pred <- prediction(prediction,as.numeric(biomod.spdis.1[,i]))
  ss <- performance(pred, "sens", "spec")
  threshold <- ss@alpha.values[[1]][which.max(ss@x.values[[1]]+ss@y.values[[1]])]
  data.threshold[i,2] <- threshold
}

######################################################

### get distributions under current and future conditions ##############

dis.current <- as.data.frame(matrix(0,1158,360))
colnames(dis.current) <- data.threshold$Species
rownames(dis.current) <- rownames(biomod.spdis.1)
for(i in 1:length(unique(spnames))){
  sp <- spnames[i]
  dat <- predicted.current[which(predicted.current$Species==sp),]
  prediction <- dat$Predicted
  dif <- prediction - data.threshold$threshold[i]
  res <- which(dif > 0)
  dis.current[res,i] <- 1
}

dis.future <- as.data.frame(matrix(0,1158,360))
colnames(dis.future) <- data.threshold$Species
rownames(dis.future) <- rownames(biomod.spdis.1)
for(i in 1:length(unique(spnames))){
  sp <- spnames[i]
  dat <- predicted.future[which(predicted.future$Species==sp),]
  prediction <- dat$Predicted
  dif <- prediction - data.threshold$threshold[i]
  res <- which(dif > 0)
  dis.future[res,i] <- 1
}
