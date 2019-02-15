if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")
if (!require("gridExtra")) install.packages("gridExtra")
library("gridExtra")



#Creating data example:
test.data <- data.frame(marker_1 = rlnorm(n=3222,log(600^2 / sqrt(850^2 + 600^2)) ,sqrt(log(1 + (850^2 / 600^2)))), marker_2 = rlnorm(n=3222,log(500^2 / sqrt(900^2 + 500^2)) ,sqrt(log(1 + (900^2 / 500^2)))), 
                        marker_3 =  rlnorm(n=3222,log(750^2 / sqrt(900^2 + 750^2)) ,sqrt(log(1 + (900^2 / 750^2)))),
                        sex = as.factor(sample(c(1,1,2), replace = T, size = 3222)), phenotype_1 = sample(c(0,0,1), replace = T, size = 3222), phenotype_2 = sample(c(0,0,0,1,1), replace = T, size = 3222),
                        yearOFbirth = sample(c(1986,1987,1988,1989,1990,1991,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004), replace = T, size = 3222))

Zstd_per_year <- function(df, phenotype, markers){
  output <- data.frame()
  phenotype_x <- which(colnames(df) == phenotype)
  levels.year <- unique(factor(df$yearOFbirth))
  for(leve in levels.year){
    year <- df[df$yearOFbirth == leve,]
    collect <- year
    controls.sub <- year[year[,phenotype_x] == "0",]
    for(k in markers){
      Zstd.marker <- data.frame( noName = (log(year[,k]) - mean(log(as.numeric(unlist(controls.sub[,k]))))) / sd(log(as.numeric(unlist(controls.sub[,k])))))
      colnames(Zstd.marker) <- paste("Zstd", k, sep = "_" )
      collect <- cbind(collect, Zstd.marker)
    }
    output <- rbind(output, collect)
  }
  
return(output)  
}

#run example
Zstd.markers.output <- Zstd_per_year(test.data, "phenotype_1", c("marker_1", "marker_2","marker_3"))


#plotting 
marker_1_org <- ggplot(Zstd.markers.output, aes(marker_1)) + geom_density()
marker_1_Zstd <- ggplot(Zstd.markers.output, aes(Zstd_marker_1)) + geom_density()
marker_2_org <- ggplot(Zstd.markers.output, aes(marker_2)) + geom_density()
marker_2_Zstd <- ggplot(Zstd.markers.output, aes(Zstd_marker_2)) + geom_density()
marker_3_org <- ggplot(Zstd.markers.output, aes(marker_3)) + geom_density()
marker_3_Zstd <- ggplot(Zstd.markers.output, aes(Zstd_marker_3)) + geom_density()
grid.arrange(marker_1_org, marker_1_Zstd, marker_2_org, marker_2_Zstd, marker_3_org, marker_4_Zstd, ncol = 2, nrow = 3)
rm()






