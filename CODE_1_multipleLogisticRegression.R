
if (!require("plyr")) install.packages("plyr")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("doParallel")) install.packages("doParallel")
if (!require("data.table")) install.packages("data.table")

set.seed(42)
library(tidyverse)
library(plyr)
library(doParallel)
library(data.table)
registerDoParallel(cores = 1)  ###set to 24+ if possible.


##creating data 
test_data <- data.frame(marker_1 = rnorm(3222, mean = -0.6, sd = 1.02), marker_2 = rnorm(3222, mean = 0.8, sd = 1.26) , sex = as.factor(sample(c(1,1,2), replace = T, size = 3222)), 
                        phenotype_1 = sample(c(0,0,1), replace = T, size = 3222), phenotype_2 = sample(c(0,0,0,1,1), replace = T, size = 3222) )


logPermut <- function(data, phenotype, marker, numberOFmodels, ID){
  #subsetting data	
  output <- data.frame()
  phenotype_x <- which(colnames(data) == phenotype)
  marker_x <- which(colnames(data)== marker)
  femaleControls <- data[data$sex == 2 & data[,phenotype_x] == 0,] #number of female(2) controls (0) 
  femaleCases <- data[data$sex == 2 & data[,phenotype_x] == 1,] #number of female(2) cases(1)
  maleControls <- data[data$sex == 1 & data[,phenotype_x] == 0,] #number of male(1) controls (0)
  maleCases <- data[data$sex == 1 & data[,phenotype_x] == 1,] #number of male(1) cases(1)
  print(paste("Female controls", nrow(femaleControls), sep = ":"))
  print(paste("Female cases", nrow(femaleCases), sep = ":"))
  print(paste("Male controls", nrow(maleControls), sep = ":"))
  print(paste("Male cases", nrow(maleCases), sep = ":"))
  #calculating female/male ratio
  maleRatio <- nrow(maleControls) / nrow(maleCases)
  print(paste("Male ratio con/cas",maleRatio, sep = ":"))
  femaleRatio <- nrow(femaleControls) / nrow(femaleCases)
  print(paste("Female ratio con/cas", femaleRatio, sep = ":"))
  #calculating sample sizes
  sampleSize <- ifelse(maleRatio < femaleRatio, floor(nrow(femaleCases) * maleRatio), floor(nrow(maleCases) * femaleRatio))
  print(paste("Sample size", sampleSize, sep = ":"))
  #model calculation x times
    output <- foreach(i = 1:numberOFmodels, .combine = "rbind.fill" ) %dopar% {
    if(femaleRatio > maleRatio){
      sampleFemaleControls <- femaleControls[sample(1:nrow(femaleControls), sampleSize, replace = F),]
      modelSet <- rbind(sampleFemaleControls, maleControls, maleCases,femaleCases)
      logitModel <- glm(formul = as.factor(modelSet[, phenotype_x]) ~ modelSet[, marker_x]  + as.factor(sex), modelSet, family = binomial)
      sampling.sex <- paste("Female", nrow(sampleFemaleControls), sep = ":")
    }else{
      sampleMaleControls <- maleControls[sample(1:nrow(maleControls), sampleSize, replace = F),]
      modelSet <- rbind(sampleMaleControls, femaleControls, maleCases, femaleCases)
      logitModel <- glm(formul = as.factor(modelSet[, phenotype_x]) ~ modelSet[, marker_x]  + as.factor(sex), modelSet, family = binomial)
      sampling.sex <- paste("Male", nrow(sampleMaleControls), sep = ":")
    }
      #Extracting results
      coxSnell  = 1 - exp(-(logitModel$null.deviance - logitModel$deviance) / length(logitModel$fitted.values))
      collect <- data.frame(P.value.marker = as.numeric(summary(logitModel)$coefficients[2,4]) ,
                            P.value.CI2.5.marker = round(as.numeric(exp(confint(logitModel))[2,1]), digits = 2) ,
                            P.value.CI97.5.marker = round(as.numeric(exp(confint(logitModel))[2,2]), digits = 2) ,
                            OR.marker = as.numeric(exp(coef(logitModel)))[2] ,
                            #Hosmer-Lemeshow :
                            Hosmer.lemshow =  1- logitModel$deviance / logitModel$null.deviance ,
                            #Cox-Snell residuals
                            coxSnell  =   coxSnell ,
                            #Nagelkerke
                            Nagelkerke = coxSnell / (1-(exp(-logitModel$deviance / length(logitModel$fitted.values)))) ,
                            P.value.constant = as.numeric(summary(logitModel)$coefficients[1,4])  , 
                            P.value.sex = as.numeric(summary(logitModel)$coefficients[3,4]) ,
                            constantCoefficient = as.numeric(summary(logitModel)$coefficients[1,1]),
                            constantSE = as.numeric(summary(logitModel)$coefficients[1,2]) ,
                            markerCoefficient = as.numeric(summary(logitModel)$coefficients[2,1]) ,
                            markerSE = as.numeric(summary(logitModel)$coefficients[2,2]) , 
                            sexCoefficient = as.numeric(summary(logitModel)$coefficients[3,1]) ,
                            sexSE = as.numeric(summary(logitModel)$coefficients[3,2]), sampleing.sex = sampling.sex)    
    }
  
  report <-  data.frame(Phenotype = ID, Marker = marker, perms = numberOFmodels, P.mean = mean(output[,1]),
  OR.mean = mean(output$OR.marker), mean.CI2.5 = mean(output$P.value.CI2.5.marker), mean.CI97.5 = mean(output$P.value.CI97.5.marker), 
  P.mean.sex = mean(output$P.value.sex), Hosmer.Lemshow = mean(output$Hosmer.lemshow),
  coxSnell = mean(output$coxSnell), Nagelkerk = mean(output$Nagelkerk), P.mean.constant = mean(output$P.value.constant), 
  constant.Coefficient = mean(output$constantCoefficient), 
  constantSE = mean(output$constantSE), marker.Coefficient = mean(output$markerCoefficient), marker.SE = mean(output$markerSE), 
  sex.Coefficient = mean(output$sexCoefficient), sex.SE = mean(output$sexSE), gender.samples = output$sampleing.sex[1] )  
  return(report)
  }  	


#examples 

marker1_phenotype1 <- logPermut(test_data, "phenotype_1", "marker_1", 10, "Diagnosis_1")
marker2_phenotype1 <- logPermut(test_data, "phenotype_1", "marker_2", 10, "Diagnosis_1")
marker1_phenotype2 <- logPermut(test_data, "phenotype_2", "marker_1", 10, "Diagnosis_2")
marker2_phenotype2 <- logPermut(test_data, "phenotype_2", "marker_2", 10, "Diagnosis_2")

Results.output <- rbindlist(mget(ls(pattern = "marker*")))
rm(marker1_phenotype1, marker1_phenotype2, marker2_phenotype1, marker2_phenotype2)
