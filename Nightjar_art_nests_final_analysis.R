rm(list=ls())
library(dplyr)
library(sp)
library(rgdal)
library(nlme)
library(carData)
library(car)
library(effects)
library(lattice)


nightjar_art_nests <- read.csv(file.choose())
#data input: Nightjar_artificial_nests_final_data.csv

Nest_Code <- unique(nightjar_art_nests$Nest_Code)
nests_df <- data.frame(Nest_Code)
#Creates dataframe with a row for each nest


for (i in 1:nrow(nests_df)){
  nest <- nightjar_art_nests[which(nightjar_art_nests$Nest_Code == nests_df$Nest_Code[i]),]
  #subsets all data relating to single nest
  
  nests_df$total_dist_metric[i] <- sqrt(sum(nest$Daily_score))
  #calculates total disturbance metric for nest
  
  nest$removal <- ifelse(nest$Removed_YN == "Y", 1, 0)
  nests_df$removal[i] <- sum(nest$removal)
  #assigns 1 for any nest where removal occured, 0 if no removal occured
  
  events <- nest[which(nest$Daily_score>0),]
  #subsets for days when disturbance recorded at nest
  
  #calculates initial risk metric for each nest
  if (nrow(events)==0){
    if (nests_df$removal[i] == 1){
      #nest does not get initial risk metric if removed before any disturbance recorded
      nests_df$initial_risk_metric[i] <- NA
    }else{
      nests_df$initial_risk_metric[i] <- 0
    }
  }else{
    first <- events[which(events$Day_of_trial == min(events$Day_of_trial)),]
    nests_df$initial_risk_metric[i] <- first$Daily_score / first$Day_of_trial
  }
}


habitat_info <- unique(nightjar_art_nests[,c(2,6,7,9:15)])
#habitat information for each nest


full_data <- merge.data.frame(habitat_info, nests_df, by="Nest_Code", all.x=T, all.y=T)
#habitat information combined with metric values for each nest into full_data dataframe

full_data$Treatment <- factor(full_data$Treatment)


#Processing of latitude and longitude data for each nest
latitude <- sapply(strsplit(as.character(full_data$GPS), split = " "), "[[", 1)
longitude <- sapply(strsplit(as.character(full_data$GPS), split = " "), "[[", 2)
latitude <- as.numeric(sapply(strsplit(sapply(strsplit(latitude, split = "[()]"), "[[", 2), split = ","), "[[", 1))
longitude <- as.numeric(sapply(strsplit(longitude, split = "[()]"), "[[", 1))
full_data$latitude <- latitude
full_data$longitude <- longitude
coords <- data.frame(x = longitude, y = latitude)
latlong <- "+init=epsg:4326"
SPDF <- SpatialPointsDataFrame(coords = coords, data = full_data, proj4string = CRS(latlong))
BNG <- spTransform(SPDF, CRS("+init=epsg:27700"))
full_data$longBNG <- BNG@coords[,1]
full_data$latBNG <- BNG@coords[,2]


#spatial autocorrelation check
nests.dists <- as.matrix(dist(cbind(full_data$longitude, full_data$latitude)))
nests.dists.inv <- 1/nests.dists
diag(nests.dists.inv) <- 0
ape::Moran.I(full_data$total_dist_metric, nests.dists.inv, na.rm=T)
ape::Moran.I(full_data$initial_risk_metric, nests.dists.inv, na.rm=T)
#Both metrics show significant spatial autocorrelation so GLS models will be used



model1_data <- full_data[which(is.na(full_data$total_dist_metric)==F),]
#Subsets data for those which have a total disturbance metric value

#Total disturbance metric model
model1 <- gls(total_dist_metric ~ Distance_to_edge_of_heath + Round_of_experiment + Distance_to_path 
              + Distance_to_building + Number_of_trees_within_20m
              + Heath_degradation_Score + Treatment, correlation=corSpher(form=~longBNG + latBNG), data=model1_data)
model1.ml <- update(model1, . ~ ., method = "ML")
#updates model to maximum likelihood so models which differ in fixed effects can be compared


shapiro.test(model1$residuals)
plot(model1)
#Checks for approximate normality of residuals

summary(model1)
#Gives effect sizes for each variable


null1 <- gls(total_dist_metric ~ Round_of_experiment + Distance_to_path 
             + Distance_to_building + Number_of_trees_within_20m
             + Heath_degradation_Score + Treatment, correlation=corSpher(form=~longBNG + latBNG), data=model1_data)
null1.ml <- update(null1, . ~ ., method = "ML")
#Drop a single model 1 variable from null 1 to obtain p value for that variable:
anova(model1.ml, null1.ml)
#likelihood ratio test



model2_data <- full_data[which(is.na(full_data$initial_risk_metric)==F),]
#Subsets data for those which have an initial risk metric value


#Initial risk metric model
model2 <- gls(initial_risk_metric ~ Distance_to_edge_of_heath + Round_of_experiment + Distance_to_path 
              + Distance_to_building + Number_of_trees_within_20m
              + Heath_degradation_Score + Treatment, correlation=corSpher(form=~longBNG + latBNG), data=model2_data)
model2.ml <- update(model2, . ~ ., method = "ML")


shapiro.test(model2$residuals)
plot(model2)
#Checks for approximate normality of residuals

summary(model2)
#Gives effect sizes for each variable


null2 <- gls(initial_risk_metric ~ Round_of_experiment + Distance_to_path 
              + Distance_to_building + Number_of_trees_within_20m
              + Heath_degradation_Score + Treatment, correlation=corSpher(form=~longBNG + latBNG), data=model2_data)
null2.ml <- update(null2, . ~ ., method = "ML")
#Drop a single model 2 variable from null 2 to obtain p value:
anova(model2.ml, null2.ml)
#likelihood ratio test






#Plotting the significant effects below:

plot(allEffects(model1, confidence.level = 0.95), selection=1, ylim=c(-0.1,5.1), 
     main = "a)",
     family="Arial Unicode MS",
     xlab="Distance to the edge of the heath (m)", 
     ylab="Total Disturbance Metric",
     axes = list(grid=T, x=list(rug = FALSE, cex=1), y=list(cex=1)),
     lwd = 4,      colors = c("black")
)
trellis.focus("panel", 1, 1)
panel.points(model1_data$Distance_to_edge_of_heath, model1_data$total_dist_metric, pch=15, col = berryFunctions::addAlpha('dodgerblue',0.3))
panel.points(model1_data$Distance_to_edge_of_heath, model1_data$total_dist_metric,pch = 0, col = berryFunctions::addAlpha('grey15',0.25))
trellis.unfocus()


plot(allEffects(model2, confidence.level = 0.95), selection=1, ylim=c(-0.1,5.1), 
     family="Arial Unicode MS",
     xlab="Distance to the edge of the heath (m)", 
     main = "b)",
     ylab="Initial Risk Metric",
     axes = list(grid=T, x=list(rug = FALSE, cex=1), y=list(cex=1)),
     lwd = 4,      colors = c("black")
)
trellis.focus("panel", 1, 1)
panel.points(model2_data$Distance_to_edge_of_heath, model2_data$initial_risk_metric, pch=15, col = berryFunctions::addAlpha('dodgerblue',0.3))
panel.points(model2_data$Distance_to_edge_of_heath, model2_data$initial_risk_metric,pch = 0, col = berryFunctions::addAlpha('grey15',0.25))
trellis.unfocus()


plot(allEffects(model2, confidence.level=0.95), ylim=c(-0.1, 5.1), selection=6, main=NULL, xlab="Heath degradation Score", 
     ylab="Initial Risk Metric",
     axes = list(grid=T, x=list(rug = FALSE, cex=1), y=list(cex=1)),
     #lines=list(lty=0),
     colors = c("black"), lwd = 4
     #symbols=list(pch=16, cex=1.65)
)
trellis.focus("panel", 1, 1)
panel.points(model2_data$Heath_degradation_Score, model2_data$initial_risk_metric, pch=15, col = berryFunctions::addAlpha('dodgerblue',0.5))
panel.points(model2_data$Heath_degradation_Score, model2_data$initial_risk_metric,pch = 0, col = berryFunctions::addAlpha('black',0.5))
trellis.unfocus()




plot(allEffects(model2, confidence.level = 0.95), selection=4, ylim=c(-0.1, 5.1), 
     family="Arial Unicode MS",
     xlab="Distance to the nearest building (m)", 
     main = NULL,
     ylab="Initial Risk Metric",
     axes = list(grid=T, x=list(rug = FALSE, cex=1), y=list(cex=1)),
     lwd = 4,      colors = c("black")
)
trellis.focus("panel", 1, 1)
panel.points(model2_data$Distance_to_building, model2_data$initial_risk_metric, pch=15, col = berryFunctions::addAlpha('dodgerblue',0.5))
panel.points(model2_data$Distance_to_building, model2_data$initial_risk_metric,pch = 0, col = berryFunctions::addAlpha('grey15',0.5))
trellis.unfocus()


plot(allEffects(model2, confidence.level = 0.95), selection=3, ylim=c(-0.1, 5.1), 
     family="Arial Unicode MS",
     xlab="Distance to the nearest path (m)", 
     main = NULL,
     ylab="Initial Risk Metric",
     axes = list(grid=T, x=list(rug = FALSE, cex=1), y=list(cex=1)),
     lwd = 4,      colors = c("black")
)
trellis.focus("panel", 1, 1)
panel.points(model2_data$Distance_to_path, model2_data$initial_risk_metric, pch=15, col = berryFunctions::addAlpha('dodgerblue',0.5))
panel.points(model2_data$Distance_to_path, model2_data$initial_risk_metric,pch = 0, col = berryFunctions::addAlpha('grey15',0.5))
trellis.unfocus()



