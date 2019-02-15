if (!require("factoextra")) install.packages("factoextra")
library(factoextra)
if (!require("ggpubr")) install.packages("ggpubr")
library("ggpubr")
if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")
if (!require("psych")) install.packages("psych")
library("psych")

#Creating data example:
test.data <- data.frame(ID = seq(1,3222,1), BDNF = rnorm(3222, mean = -0.6, sd = 1.02), sex = as.factor(sample(c(1,1,2), replace = T, size = 3222)), 
                        yearOFbirth = sample(c(1986,1987,1988,1989,1990,1991,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004), replace = T, size = 3222),
                        SNP.1 = sample(c(0,1,2), replace = T, size = 3222), SNP.2 = sample(c(0,0,2), replace = T, size = 3222), 
                        SNP.3 = sample(c(0,1,2,1), replace = T, size = 3222), SNP.4 = sample(c(1,2,2), replace = T, size = 3222))

#normalizing yearOFbirth
test_data$yearOFbirth <- (test_data$yearOFbirth - min(test_data$yearOFbirth)) / (max(test_data$yearOFbirth)-min(test_data$yearOFbirth))

#Function for finding all combinations of SNPs. sex and normalized yearOFbirth are in all models. 
lmAllcombi <- function(df){
  SNPs <- as.character(colnames(df[-c(1:4)]))
  combinations <- unlist(lapply(1:length(SNPs), function(x) combn(SNPs,x,FUN = paste, collapse = "+")))
  allCombinations <- data.frame(A = "BDNF", B = "sex", C = "yearOFbirth", SNP = combinations)
  allCombinations$formula <- paste0(allCombinations$A, "~", allCombinations$B, "+" ,allCombinations$C , "+", allCombinations$SNP)
  output <- data.frame()
  for(i in 1:nrow(allCombinations)){
    lmFormula <- as.formula(allCombinations$formul[i])
    lmModel <- lm(lmFormula, data = df, na.action = na.omit)
    coefficientValue <- data.frame(P = summary(lmModel)$coefficients[,4])
    lowest.SNP.P <- coefficientValue[4:nrow(coefficientValue),][which.min(coefficientValue$P)]
    SNPS.p <-  data.frame(P.SNPs = coefficientValue[4:nrow(coefficientValue),])
    sig.SNPs <- length(SNPS.p[SNPS.p$P.SNPs <= 0.01,])
    print(allCombinations$formul[i])
    collect <- data.frame( R.adj = summary(lmModel)$adj.r.squared, F.stats.p = pf(summary(lmModel)$fstatistic[1], summary(lmModel)$fstatistic[2], summary(lmModel)$fstatistic[3]),
                           AIC = AIC(lmModel), number.of.SNPs = nrow(coefficientValue) - 3  , lowest.SNP.P = lowest.SNP.P, sex.p.value =  summary(lmModel)$coefficients[2,4], 
                           number.sig.SNPs =  sig.SNPs, formula = allCombinations$formula[i])
    output <- rbind(output, collect)
    }
  return(output)
}


#example run
lmModels_output <- lmAllcombi(test.data)


#Creating data example  
test.data.genetics <- data.frame(ID = seq(1,1567,1), sex = as.factor(sample(c(1,1,2), replace = T, size = 1567)), 
                                 BDNF = rnorm(1567, mean = -0.6, sd = 1.02),
                                 SNP.1 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.2 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.3 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.4 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.5 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.6 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.7 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.8 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.9 = sample(c(rep(0,1189), rep(1,117), rep(2,261))),
                                 SNP.10 = sample(c(rep(0,1189), rep(1,17), rep(2,361))),
                                 SNP.11 = sample(c(rep(0,1189), rep(1,17), rep(2,361))))

#calculat PCs
PCA <- prcomp(test.data.genetics[,4:length(test.data.genetics)])
summary(PCA) 
df.PCA <- data.frame(PCA$x[,1:2])  
ggplot(df.PCA, aes(PC1, PC2)) + geom_point()  
#fining optimal k.  
fviz_nbclust(df.PCA, kmeans, method = "wss")
#clustering
kmeans.df.PCA <- kmeans(df.PCA, 4 , nstart = 100, iter.max = 1000)  
#generating plot:
plot.PCA.kmeans <- fviz_cluster(kmeans.df.PCA, data = df.PCA, ellipse.type = "euclid", labelsize = 0, ggtheme = theme_minimal())
ggpar(plot.PCA.kmeans)

#adding cluster ID variable  
test.data.genetics$cluster <-  as.factor(kmeans.df.PCA$cluster)
#and dummy codes
test.data.genetics.dummy <- cbind(test.data.genetics ,dummy.code(test.data.genetics$cluster) )
colnames(test.data.genetics.dummy)[names(test.data.genetics.dummy) %in% c("1","2","3","4")] <- c("cluster.1", "cluster.2","cluster.3","cluster.4")
#lm model
lmModel <- lm(formula = BDNF ~ cluster.1 + cluster.2 + cluster.3  + cluster.4 + as.factor(sex), data = test.data.genetics.dummy)
summary(lmModel)
lmfactorModel <- lm(formula = BDNF ~ as.factor(cluster) + as.factor(sex), data = test.data.genetics.dummy)
summary(lmfactorModel)

