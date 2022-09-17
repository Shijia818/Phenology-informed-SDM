
#######download species distribution data from GBIF #######
setwd("D:/phenology/Data")
spnames <- read.csv("species_list.csv")
library(spocc)
res <- occ(query = as.character(spnames$Species),from=c("gbif"),limit=100000)

####how to deal with GBIF data#####
library(data.table)
setwd("D:/phenology/Data/distribution_data/Trillium grandiflorum")
data <- data.table::fread("GBIF.csv",sep='\t',quote = '',encoding = "UTF-8")
data.new <- data[,c("species","scientificName","countryCode","locality","stateProvince",
                    "decimalLongitude","decimalLatitude","basisOfRecord","year")]

data.table::fwrite(data.new,"GBIF.1.csv")

######
setwd("D:/phenology/Climate/Future Climate/Tiff/ssp585/INM-CM4-8")
library(raster)
r <-stack("D:/phenology/Climate/Future Climate/Tiff/ssp585/WorldClim/wc2.1_2.5m_bioc_INM-CM4-8_ssp585_2061-2080.tif")
nlayers(r)
for(i in 1:nlayers(r)){
  band<-r[[i]]
  writeRaster(band,paste('D:/phenology/Climate/Future Climate/Tiff/ssp585/INM-CM4-8/','band',i,'.tif', sep=''))
}
allfile = dir() 
txtfile <- grep("*.tif.aux.xml", allfile)
file.remove(allfile[txtfile])

######
setwd("D:/phenology/Analyses")
data.phenology <- read.csv("data.phenology.csv")
PDflr <- data.phenology[which(data.phenology$pfl >= 0.5),] ### peak flowering ####
PDfrt <- data.phenology[which(data.phenology$pfr >= 0.5),] ### peak fruiting ####
PDbud <- data.phenology[which(data.phenology$pbd >= 0.5),] ## peak budding ####
FDflr <- data.phenology[which(data.phenology$pfl < 0.5 & data.phenology$fruit == 0),] #### first flowering####
FDflr.final <- FDflr[which(FDflr$flower >= 1 | FDflr$bud >= 1),]
FDfrt <- data.phenology[which(data.phenology$pfr < 0.5 & data.phenology$fruit >= 1),]   #### first fruiting ####
FDbud <- data.phenology[which(data.phenology$pbd < 0.5 & data.phenology$flower == 0),]  #### first budding ####

#### match climate data ####

#### Spring total precipitation ####
setwd("D:/work/results/tables/pdbud/csv/Precipitation")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
length(Year)

data.pre <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.pre) <- c("Longitude","Latitude","Year","Precipitation")
data.pre$Year <- PDbud$year
data.pre$Longitude <- PDbud$longitude
data.pre$Latitude <- PDbud$latitude

##########
sequence <- seq(6,1482,12)

######
for(j in 1:length(sequence)){
      res <- which(data.pre$Year == Year[j])
      pos <- sequence[j]
      data1 <- read.csv(files[pos])
      data2 <- read.csv(files[pos+1])
      data3 <- read.csv(files[pos+2])
      pre.total <- data1$PRISM_ppt_ + data2$PRISM_ppt_ + data3$PRISM_ppt_ 
      data.pre[res,4] <- pre.total
}

write.csv(data.pre, file="D:/phenology/Climate/County/PDbud/data.precipitaion.csv")

##### Spring mean air temperature #####
setwd("D:/work/results/tables/pdbud/csv/Tmean")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
###length(Year)

data.Tmean <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.Tmean) <- c("Longitude","Latitude","Year","Tmean")
data.Tmean$Year <- PDbud$year
data.Tmean$Longitude <- PDbud$longitude
data.Tmean$Latitude <- PDbud$latitude

##########
sequence <- seq(6,1482,12)

######
for(j in 1:length(sequence)){
  res <- which(data.Tmean$Year == Year[j])
  pos <- sequence[j]
  data1 <- read.csv(files[pos])
  data2 <- read.csv(files[pos+1])
  data3 <- read.csv(files[pos+2])
  Tmean <- (data1$PRISM_tmea + data2$PRISM_tmea + data3$PRISM_tmea)/3
  data.Tmean[res,4] <- Tmean
}

write.csv(data.Tmean, file="D:/phenology/Climate/County/PDbud/data.Tmean.csv")


####### mean annual temperature ####
setwd("D:/work/results/tables/pdbud/csv/Tmean")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
data.bio1 <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.bio1) <- c("Longitude","Latitude","Year","bio1")
data.bio1$Year <- PDbud$year
data.bio1$Longitude <- PDbud$longitude
data.bio1$Latitude <- PDbud$latitude

####3
sequence <- seq(1,1488,12)

####
#############
for(j in 1:length(sequence)){
  res <- which(data.bio1$Year == Year[j])
  pos <- sequence[j]
  data1 <- read.csv(files[pos])
  data2 <- read.csv(files[pos+1])
  data3 <- read.csv(files[pos+2])
  data4 <- read.csv(files[pos+3])
  data5 <- read.csv(files[pos+4])
  data6 <- read.csv(files[pos+5])
  data7 <- read.csv(files[pos+6])
  data8 <- read.csv(files[pos+7])
  data9 <- read.csv(files[pos+8])
  data10 <- read.csv(files[pos+9])
  data11 <- read.csv(files[pos+10])
  data12 <- read.csv(files[pos+11])
  data.total <- cbind(data1$PRISM_tmea,data2$PRISM_tmea,data3$PRISM_tmea,data4$PRISM_tmea,data5$PRISM_tmea,data6$PRISM_tmea,data7$PRISM_tmea,
                      data8$PRISM_tmea,data9$PRISM_tmea,data10$PRISM_tmea,data11$PRISM_tmea,data12$PRISM_tmea)
  bio1 <- apply(data.total,1,mean)
  data.bio1[res,4] <- bio1
}

write.csv(data.bio1, file="D:/phenology/Climate/County/PDbud/data.bio1.csv")

##### temperature seasonality #####
setwd("D:/work/results/tables/pdbud/csv/Tmean")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
###length(Year)

data.bio4 <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.bio4) <- c("Longitude","Latitude","Year","bio4")
data.bio4$Year <- PDbud$year
data.bio4$Longitude <- PDbud$longitude
data.bio4$Latitude <- PDbud$latitude

####
sequence <- seq(1,1488,12)

#############

#############
for(j in 1:length(sequence)){
  res <- which(data.bio4$Year == Year[j])
  pos <- sequence[j]
  data1 <- read.csv(files[pos])
  data2 <- read.csv(files[pos+1])
  data3 <- read.csv(files[pos+2])
  data4 <- read.csv(files[pos+3])
  data5 <- read.csv(files[pos+4])
  data6 <- read.csv(files[pos+5])
  data7 <- read.csv(files[pos+6])
  data8 <- read.csv(files[pos+7])
  data9 <- read.csv(files[pos+8])
  data10 <- read.csv(files[pos+9])
  data11 <- read.csv(files[pos+10])
  data12 <- read.csv(files[pos+11])
  data.total <- cbind(data1$PRISM_tmea,data2$PRISM_tmea,data3$PRISM_tmea,data4$PRISM_tmea,data5$PRISM_tmea,data6$PRISM_tmea,data7$PRISM_tmea,
                      data8$PRISM_tmea,data9$PRISM_tmea,data10$PRISM_tmea,data11$PRISM_tmea,data12$PRISM_tmea)
  bio4 <- apply(data.total,1,sd)
  data.bio4[res,4] <- bio4 *100
}

write.csv(data.bio4, file="D:/phenology/Climate/County/PDbud/data.bio4.csv")

#########
##### precipitation seasonality #####
setwd("D:/work/results/tables/pdbud/csv/Precipitation")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
###length(Year)

data.bio15 <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.bio15) <- c("Longitude","Latitude","Year","bio15")
data.bio15$Year <- PDbud$year
data.bio15$Longitude <- PDbud$longitude
data.bio15$Latitude <- PDbud$latitude

sequence <- seq(1,1488,12)

#############
for(j in 1:length(sequence)){
  res <- which(data.bio15$Year == Year[j])
  pos <- sequence[j]
  data1 <- read.csv(files[pos])
  data2 <- read.csv(files[pos+1])
  data3 <- read.csv(files[pos+2])
  data4 <- read.csv(files[pos+3])
  data5 <- read.csv(files[pos+4])
  data6 <- read.csv(files[pos+5])
  data7 <- read.csv(files[pos+6])
  data8 <- read.csv(files[pos+7])
  data9 <- read.csv(files[pos+8])
  data10 <- read.csv(files[pos+9])
  data11 <- read.csv(files[pos+10])
  data12 <- read.csv(files[pos+11])
  data.total <- cbind(data1$PRISM_ppt_,data2$PRISM_ppt_,data3$PRISM_ppt_,data4$PRISM_ppt_,data5$PRISM_ppt_,data6$PRISM_ppt_,data7$PRISM_ppt_,
                      data8$PRISM_ppt_,data9$PRISM_ppt_,data10$PRISM_ppt_,data11$PRISM_ppt_,data12$PRISM_ppt_)
  mean <- apply(data.total,1,mean)
  sd <- apply(data.total,1,sd)
  data.bio15[res,4] <- sd/mean *100
}

write.csv(data.bio15, file="D:/phenology/Climate/County/PDbud/data.bio15.csv")

##### mean annual precipitation #####
setwd("D:/work/results/tables/pdbud/csv/Precipitation")
PDbud <- read.csv("D:/phenology/Analyses/PDbud.csv")
files <- list.files()
Year <- c(1895:2018)
###length(Year)

data.bio12 <- as.data.frame(matrix(NA,nrow(PDbud),4))
colnames(data.bio12) <- c("Longitude","Latitude","Year","bio12")
data.bio12$Year <- PDbud$year
data.bio12$Longitude <- PDbud$longitude
data.bio12$Latitude <- PDbud$latitude

sequence <- seq(1,1488,12)

#############
for(j in 1:length(sequence)){
  res <- which(data.bio12$Year == Year[j])
  pos <- sequence[j]
  data1 <- read.csv(files[pos])
  data2 <- read.csv(files[pos+1])
  data3 <- read.csv(files[pos+2])
  data4 <- read.csv(files[pos+3])
  data5 <- read.csv(files[pos+4])
  data6 <- read.csv(files[pos+5])
  data7 <- read.csv(files[pos+6])
  data8 <- read.csv(files[pos+7])
  data9 <- read.csv(files[pos+8])
  data10 <- read.csv(files[pos+9])
  data11 <- read.csv(files[pos+10])
  data12 <- read.csv(files[pos+11])
  data.total <- cbind(data1$PRISM_ppt_,data2$PRISM_ppt_,data3$PRISM_ppt_,data4$PRISM_ppt_,data5$PRISM_ppt_,data6$PRISM_ppt_,data7$PRISM_ppt_,
                      data8$PRISM_ppt_,data9$PRISM_ppt_,data10$PRISM_ppt_,data11$PRISM_ppt_,data12$PRISM_ppt_)
  bio12 <- apply(data.total,1,mean)
  data.bio12[res,4] <- bio12
}

write.csv(data.bio12, file="D:/phenology/Climate/County/PDbud/data.bio12.csv")


##### match growth_form #####
#### setwd("D:/phenology/Analyses")
#### PDbud <- read.csv("PDbud.csv")
#### data <- read.csv("D:/phenology/Data/datasets/Growth_form.csv")
#### PDbud.1 <- PDbud[which(PDbud$species %in% data$Species),]
#### data.1 <- data[which(data$Species %in% PDbud.1$species),]
### data.2 <-data.1[match(PDbud.1$species,data.1$Species),]


###### construct relationship between DOY and climate ####
#### Previous studies provides background variables known to affect plant phenology ####
##### We extract spring mean temperature, bio4, spring total precipitation, bio15 in the county and year that the specimen was collected ####
library(lme4)
library(lmerTest)
library(dplyr)
library(car)
library(performance)
library(ape)
library(phyr)

setwd("D:/phenology/Analyses")
data <- read.csv("PDbud.csv")
Growth_form <- read.csv("D:/phenology/Data/datasets/Growth_form.csv")

##### all variables are scaled to have mean 0 and standard deviation 1 ###

######
d <- mutate_at(data, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))

pheno.model <- lmer(julian ~ bio1 + bio12 + bio4 + bio15 + Growth_form  + year + bio1:Growth_form + 
                      bio12:Growth_form + bio4:Growth_form + bio15:Growth_form + (1 | latitude) + (1 | longitude)+
                      (1 | species) +
                      (1 + bio1 | species) + 
                      (1 + bio12 | species) + 
                      (1 + bio4 | species) + 
                      (1 + bio15 | species), data = d, REML =T,
                    control = lmerControl(optimizer = "bobyqa", 
                                          optCtrl = list(maxfun = 2e5)))
summary(pheno.model)

Anova(pheno.model)
saveRDS(pheno.model,file = "D:/phenology/Results/Phenology_model/PDbud.rds")

#### PDbuds_model <- readRDS("D:/phenology/Results/Phenology_model/PDbud.rds") #####


###### Compute the marginal and conditional r-squared value for mixed-model####
r2(pheno.model)
model_performance(pheno.model)


######## 

###### predict phenology at grid-cell level for each species ######
current.climate <- read.csv("D:/phenology/Climate/Current Climate/current_climate.csv")
current.climate.1 <- current.climate[,c(1,2,5,13,16)]
grid40km <-read.csv("D:/phenology/Data/datasets/grid40km.csv")
current.climate.2 <- left_join(grid40km, current.climate.1,by = "Gridcode")

#### head(current.climate.1)

### construct dataframe ###
Species <- unique(data$species)

test <- rbind()
for (i in 1:116){
  test <- rbind(test,current.climate.2)}

Species.rev <- rep(Species,1158) %>% sort() %>% cbind(test)
colnames(Species.rev) [1] <- "species"
Species.rev <- Species.rev[,-2]

Growth_form.1 <- Growth_form[match(Species.rev$species,Growth_form$Species),]
Species.rev$Growth_form <- Growth_form.1$Growth_form
Species.rev$year <- "1970"
Species.rev$year <- as.integer(Species.rev$year)


######
model.current <- mutate_at(Species.rev, .vars = vars(bio1, bio12, bio4, bio15),.funs = function(x)(x-mean(x))/sd(x))
predicted.current <- predict(pheno.model, model.current, type="response",allow.new.levels=T)
model.current$predict <- predicted.current
write.csv(model.current,file = "D:/phenology/Results/Phenology_model/PDbuds_current.csv")





















#####PCA example#####
data.current <- read.csv("D:/phenology/current.csv",row.names=1)
results <- prcomp(data.current,scale = TRUE)
summary(results)
results$x <- -1*results$x
head(results$x)

####Note that your predictors have very different ranges,so if we want to directly compare the effects or test for interactions
##we need to first rescale them or conduct PCA


#####model fitting#####
library(pROC)

####Species distributions####
setwd("D:/phenology")
sp_dis <- read.csv("./Data/distribution_data/distribution_data.csv",row.names=1)
sp_dis_1 <- sp_dis[,1,drop=FALSE]

####Climate variables####
data.current.climate <- read.csv("./Climate/Current Climate/current_climate.csv",row.names=1)
dat <- data.current.climate[,c(4,10,11,15,18,19)]

#####Anemone canadensis####
dataset <- cbind(sp_dis_1,data.current.climate)
m1 <- glm(Anemone.canadensis ~bio4 + bio10 + bio11 + bio15 + bio18 + bio19, family = "binomial",data = dataset)

#####Model evaluation###
##### how to calculate AUC###
predicted <- predict(m1, dat, type="response") ###predictions under current conditions#####

#####AUC value####
AUC <- auc(dataset$Anemone.canadensis,predicted)

####if we want to divide data into train and test, respectively#####
sample <- sample(c(TRUE, FALSE), nrow(dataset), replace=TRUE,prob = c(0.7, 0.3))
train <- dataset[sample,]
test <- dataset[!sample, ] 
model <- glm(Anemone.canadensis ~ bio4 + bio10 + bio11 + bio15 + bio18 + bio19, family = "binomial",data = train)
predicted <- predict(model, test, type="response")
auc(test$Anemone.canadensis, predicted)

#######predict future probability of occurrence######

####future climate####
data.future.climate <- read.csv("./Climate/Future Climate/csv/ssp585/ACCESS-CM2/ACCESS-CM2.csv",row.names=1)
dat.future <- data.future.climate[,c(4,10,11,15,18,19)]
predicted.future<- predict(m1, dat.future, type="response")

######Generalized linear mixed model#####
#####add random component#####
library(lme4) 
library(reshape)#for melt function#####
setwd("D:/phenology/Analyses") #for lmer function#####
data.dis <- read.csv("distribution_data.csv")
climate.current <- read.csv("current_climate.csv")

### climate.current$bio4 <- (climate.current$bio4 - mean(climate$bio4)) / sd(climate$bio4)


mydata <- melt(data.dis, id.vars = "Gridcode")
colnames(mydata)[2:3] <- c("Taxon", "Present")
model.data <- merge(mydata, climate.current, by="Gridcode", all.x = TRUE) 
###merge presence/absence with environmental data#####


##########fitting model with generalized linear mixed model####
###full.model <- bglmer(Present ~ bio4 + bio10 + bio11 + bio15 + bio18 + bio19 +
                       #(bio4 + bio10 + bio11 + bio15 + bio18 + bio19 | Taxon), family = binomial, fixef.prior = normal(sd = 1),
                     #cov.prior = wishart, data = model.data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6))

######
full.model <- glmer(Present ~ bio4 + bio10 + bio11 + bio15 + bio18 + bio19 + (1 + bio4 + bio10 + bio11 + bio15 + bio18 + bio19 | Taxon),
                    data = model.data, family = binomial)
###
model.cov <- model.data[, c(2, 7, 13, 14, 11, 18, 21, 22)]

####
predicted <- predict(full.model, model.cov, type="response")
AUC <- auc(model.data$Present,predicted)

######## basic phenology variable #####
####
library(stringr)
setwd("D:/phenology/Analyses")
data.dis <- read.csv("distribution_data.csv",row.names = 1)
spnames <- colnames(data.dis)
spnames.col <- gsub("\\.","_",spnames)

###
data.phenology <- read.csv("D:/phenology/Phenology/Park.csv")
pheno <- subset(data.phenology, species %in% spnames.col)
####length(unique(pheno$species))

##### Formats columns ####

pheno$date <- as.Date(paste(pheno$year,"-",pheno$month,"-",pheno$day,sep=""), format="%Y-%m-%d")
pheno$julian <- as.numeric(format(pheno$date,format="%j")) #### calculate day of year (DOY)###
pheno$state_county <- paste(pheno$state,pheno$county,sep="-")
#### Multiple States may share the same county  name ###
#head(pheno)
pheno <- pheno[, c(2, 4, 5, 6, 13, 17, 18)]


###### construct relationship between DOY and climate ####
#### Previous studies provides background variables known to affect plant phenology ####
##### We extract bio1, bio4, bio12,bio15 in the county and year that the specimen was collected ####

pheno.model <- lmer(DOY ~ bio1 + bio4 + bio12 + bio15 + growth_form + (1 | year) + (1 | county) + (1 | species) +
                       (1 + bio1 | species) + 
                       (1 + bio4 | species) + 
                       (1 + bio12 | species) + 
                       (1 + bio15 | species), data = data, REML =T,
                       control = lmerControl(optimizer = "bobyqa", 
                       optCtrl = list(maxfun = 2e5)))

##### Then we predict the phenology under current (1960 - 2010) and future conditions (2070s) #####




##### create data frame of species-county phenology #####
model.phenology <- as.data.frame(matrix(NA, length(unique(pheno$species)),5))
colnames(model.phenology) <- c("species", "location","buds", "flowers", "fruits")

###### 
for(i in 1:length(unique(pheno$species))){
  sp.names <- unique(pheno$species)[i]
     dat.species <- subset(pheno, species %in% sp.names)
        for(j in 1:length(unique(dat.species$state_county ){
          county.names <- unique(dat.species$state_county)[j]
          dat.county <- subset(dat.species, state_county %in% county.names)
                phe.min <- 
    

                  boxplot(as.numeric(data.1$pct_flo))
                > hist(as.numeric(data.1$pct_flo))
                > quantile(as.numeric(data.1$pct_flo), probs=c(.1,.5))























