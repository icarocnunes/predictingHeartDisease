# Projeto 1: Heart Disease UCI (Kaggle)
# Desenvolvido por Icaro Caze Nunes

# Windpows
setwd("~/dataProjects/Heart/Data") # Remember to change to a local folder that contains the data file 

# 1. Importing:


# 1.1. Packages:

library(dplyr)
library(tidyverse)
library(funModeling)
library(caret)
library(ggplot2)
library(viridis)
library(ggridges)
library (hrbrthemes)
library(ggpubr)
library(dplyr)

# 1.2. Functions:

# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })


geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# 1.3. Database:

heartData <- read.csv("processed.cleveland.data", header = F, stringsAsFactor = TRUE)

## Setting Column names:
colnames(heartData) <- c('age', 'sex', 'chestPain', 'bloodPressure', 'cholesterol', 'bloodSugar', 'ecg',
                         'maxHeartRate', 'indAngina', 'stDepression', 'stSlope', 'numVessels', 'scintigraphy', 'condition')

## Data checking:
summary(heartData)

# 2. Correcting data Mistakes

# 2.1. Deleting rows with wrong values:

heartData <- heartData[!(heartData$numVessels == '?' | heartData$scintigraphy == "?"),]

# 2.2. Sex (0- Female; 1- Male):

heartData <- heartData %>% 
  mutate(sex = ifelse(sex == 1,'Male','Female'))

heartData$sex <- as.factor(heartData$sex)

# 2.3. Chest Pain Type (1- Typical Angina; 2- Atypical Angina; 3- Non-Anginal Pain; 4- Asymptomatic):

heartData <- heartData %>% 
  mutate(chestPain = ifelse(chestPain == 1,'Typical Angina', ifelse(chestPain == 2, 'Atypical Angina',
                                                                    ifelse(chestPain == 3,'Non-Anginal Pain', 'Asymptomatic'))))

heartData$chestPain <- as.factor(heartData$chestPain)

# 2.4. Fasting Blood Sugar > 120 mg/dl (0- False; 1- True):

heartData <- heartData %>% 
  mutate(bloodSugar = ifelse(bloodSugar == 1,'> 120 mg/dl','< 120 mg/dl'))

heartData$bloodSugar <- as.factor(heartData$bloodSugar)

# 2.5. Resting Electrocardiographic Results Value (0- Normal; 1- having ST-T wave abnormality [T wave inversions and/or 
#        ST elevation or depression of > 0.05 mV]; 2- Showing probable or definite left ventricular hypertrophy by Estes' criteria):

heartData <- heartData %>% 
  mutate(ecg = ifelse(ecg == 0,'Normal', 
                      ifelse(ecg == 1, 'ST-T wave abnormality', 
                             'Left ventricular hypertrophy')))

heartData$ecg <- as.factor(heartData$ecg)

# 2.6. Exercise Induced Angina (0- No; 1- yes):

heartData <- heartData %>% 
  mutate(indAngina = ifelse(indAngina == 1,'Yes','No'))

heartData$indAngina <- as.factor(heartData$indAngina)

# 2.7. The Slope of the Peak Exercise ST Segment (1- Upsloping; 2- Flat; 3- Downsloping):

heartData <- heartData %>% 
  mutate(stSlope = ifelse(stSlope == 1,'Upsloping', 
                          ifelse(stSlope == 2, 'Flat', 'Downsloping')))

heartData$stSlope <- as.factor(heartData$stSlope)

# 2.8. Scintigraphy (3- Normal; 6- Fixed defect; 7- Reversable defect):

heartData <- heartData %>% 
  mutate(scintigraphy = ifelse(scintigraphy == '3.0', 'Normal', 
                               ifelse(scintigraphy == '6.0', 'Fixed Defect', 'Reversable Defect')))
heartData$scintigraphy <- as.factor(heartData$scintigraphy)

# 2.9. Number of Major Vessels Colored by Flourosopy:

heartData$numVessels <- factor(heartData$numVessels)


# 2.10. Diagnosis of Heart Disease (0- < 50% diameter narrowing; 1- > 50% diameter narrowing):

heartData <- heartData %>% 
  mutate(condition = ifelse(condition == 0,'CAD-','CAD+'))
heartData$condition <- as.factor(heartData$condition)

## Data checking:
glimpse(heartData)

summary(heartData)

# 3. Exploratory Data Analysis

# 3.1. Number of Patients with and without Heart Disease:

countCard <- count(heartData, condition)
percCard <- round((countCard$n*100)/sum(countCard$n),2) 
pieDF <- data.frame(stats = countCard$condition, n = countCard$n, perc = percCard)

pieChart <- ggplot(data = pieDF, aes(x = "", y = perc, fill = stats)) +
  geom_bar(stat = "identity") + coord_polar("y", start = 0) +
  scale_fill_manual(values = c("darkturquoise","salmon"), labels = c("Without Coronary Artery Disease",
                                                                     "With Coronary Artery Disease"))
blank_theme <- theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18))

pieChart <- pieChart + blank_theme + coord_polar("y",start=0) + 
  geom_text(aes(label = paste0(perc,"%")), position = position_stack(vjust = 0.5))

pieChart + theme(axis.text.x=element_blank()) + 
  labs(fill = "", title = 'Figure 1')

# 3.2. Age:

ggplot(heartData, aes(x=age, fill=condition)) +
  geom_histogram( alpha=0.6, position = 'identity', binwidth=1) +
  scale_fill_manual(values=c("darkturquoise","salmon")) +
  labs(fill="Diagnosis", x= 'Age', y = 'Number of patients', title = 'Figure 2') +
  scale_y_continuous(limits = c(0, 12), n.breaks = 8)

# 3.3. Sex:

# 3.3.1. Number of patients grouped by sex:

ggplot(heartData, aes(sex, fill = condition)) + 
  geom_bar() + labs(fill = "Diagnosis", x = "Sex", y = "Number of patients") +
  labs(title = 'Figure 3') +
  scale_fill_manual(values=c("darkturquoise","salmon"))

# 3.3.2. Probability of occurrence (Calculation):

pMale <- (count(heartData[heartData$condition == 'CAD+' & heartData$sex == 'Male', ])/
            count(heartData[heartData$sex == 'Male',]))
pFemale <- (count(heartData[heartData$condition == 'CAD+' & heartData$sex == 'Female', ])/
              count(heartData[heartData$sex == 'Female',]))

p <- cbind(rbind(pFemale,pMale), levels(heartData$sex))
colnames(p) <- c('prob', 'sex')

# 3.3.3. Probabilidade de ocorrência (Plot):

ggplot(data = p, aes(x = sex, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity") + 
  labs(fill="", title="Figure 4") +  
  theme(axis.text.x = element_text( vjust=0.6))  + ylab("CAD incidence rate") +
  xlab("Sex") +   scale_y_continuous(limits = c(0, 0.6), n.breaks = 4)

# 3.4. Cholesterol:


# 3.4.1. Cholesterol x Diagnosis x Age x Sex:

ggplot(heartData, aes(x = age, y = cholesterol, size = sex)) +
  geom_point(aes(color = condition)) + geom_point(alpha = 0.1)+ xlab("Age group") +
  ylab("Cholesterol (mg/dl)") + scale_color_manual(values=c("limegreen", "red2")) +
  labs(color = "Diagnosis", size = "Sex", title = 'Figure 5') +
  geom_hline(yintercept=200, linetype="dashed", color = "black", size=1)

# 3.5. Blood Pressure and Heart Rate:

ggplot(heartData, aes(x = maxHeartRate, y = bloodPressure, size = sex)) +
  geom_point(aes(color = condition),  alpha = 0.7) + xlab("Maximum Heart Rate (bpm)") +
  ylab("Blood Pressure (mmHg)") + scale_color_manual(values=c("limegreen", "red2")) +
  labs(color = "Diagnosis", size = "Sex", alpha = '', title = 'Figure 6')

# 3.6. Fasting Blood Sugar:

# 3.6.1. Number of patients grouped:

ggplot(heartData, aes(bloodSugar, age, fill = condition )) + geom_split_violin(trim = TRUE) +
  scale_fill_manual(values=c("darkturquoise","salmon")) +   xlab('Blood Sugar (mg/dL)') + ylab("Age") +
  theme(axis.text.x=element_text(angle=30, vjust=.8, hjust=0.8)) + labs(title = 'Figure 7')

# 3.6.1. Probability of occurrence (Calculation):

desire <- (heartData %>%  filter(bloodSugar == '< 120 mg/dl' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(bloodSugar == '< 120 mg/dl') %>% summarise(n()))
heartData %>%  filter(bloodSugar == '> 120 mg/dl'  & condition == 'CAD+') %>% summarise(n())

desire <- (heartData %>%  filter(bloodSugar == '< 120 mg/dl' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(bloodSugar == '< 120 mg/dl') %>% summarise(n()))

high <- (heartData %>%  filter(bloodSugar == '> 120 mg/dl' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(bloodSugar == '> 120 mg/dl') %>% summarise(n()))

sugarLevel  <- cbind(rbind(desire, high), levels(heartData$bloodSugar))
colnames(sugarLevel) <- c('prob', 'pSugar')

# 3.6.1. Probability of occurrence (Plot):

ggplot(data = sugarLevel, aes(x = pSugar, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity", width = 0.5) + labs(fill="", title = 'Figure 8') +  
  theme(axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8))  + ylab("Maximum Heart Rate") +
  xlab("Blood Sugar (mg/dl)") +   scale_y_continuous(limits = c(0, 0.6), n.breaks = 4)

# 3.7. Exercise Induced Angina:

# 3.7.1. Number of patients grouped by:

ggplot(heartData, aes(x = indAngina, fill = condition)) + geom_bar() + 
  labs(fill="Diagnosis", title = 'Figure 9') + 
  scale_fill_manual(values=c("darkturquoise","salmon")) + 
  xlab("Exercise Induced Angina") + ylab("Number of patients")+
  theme(axis.text.x=element_text( angle=30, vjust=.8, hjust=0.8))

# 3.7.2. Probability of occurrence (Calculation):

anginaYes <- (heartData %>%  filter(indAngina == 'Yes' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(indAngina == 'Yes') %>% summarise(n()))

anginaNo <- (heartData %>%  filter(indAngina == 'No' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(indAngina == 'No') %>% summarise(n()))

pAngina  <- cbind(rbind(anginaNo, anginaYes), levels(heartData$indAngina))
colnames(pAngina) <- c('prob', 'angina')

# 3.7.3. Probability of occurrence (Plot):

ggplot(data = pAngina, aes(x = angina, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity") + labs(fill="", title = 'Figure 10') +  
  theme(axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8))  + ylab("CAD incidence rate") +
  xlab("Exercise Induced Angina") 

# 3.8. ECG:

ggplot(heartData, aes(x = age, y = ecg, fill = ecg))+
  geom_density_ridges(aes(point_color = ecg, point_fill = ecg, point_shape = ecg),
                      alpha = .2, jittered_points = TRUE) + scale_point_color_hue(l = 50) +
  scale_discrete_manual(aesthetics = 'point_shape', values = c(21,22,23)) +
  theme(legend.position = "none") + scale_x_continuous(limits = c(20, 90), n.breaks = 8) + facet_grid(.~condition) +
  xlab('Age group') + ylab('ECG Results') + labs(title = 'Figure 11')

ecgNormal <- (heartData %>%  filter(ecg == 'Normal' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(ecg == 'Normal') %>% summarise(n()))

ecgHyper <- (heartData %>%  filter(ecg == 'Left ventricular hypertrophy' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(ecg == 'Left ventricular hypertrophy') %>% summarise(n()))

ecgAnorm <- (heartData %>%  filter(ecg == 'ST-T wave abnormality' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(ecg == 'ST-T wave abnormality') %>% summarise(n()))

ecgRate  <- cbind(rbind(ecgAnorm, ecgHyper, ecgNormal), levels(heartData$ecg))
colnames(ecgRate) <- c('prob', 'ecg')

ggplot(data = ecgRate, aes(x = ecg, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity", width = 0.6) + labs(fill="", title = 'Figure 12') +  
  theme(axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8))  + ylab("CAD incidence rate") +
  xlab("ECG resuts") 

# 3.9. Chest pain type and Depression induced by exercise:

ggline(heartData, x = "chestPain", y = "stDepression", add = c("mean_sd", "jitter"),
       color = "condition", palette = c("darkturquoise","salmon")) +
  labs(x = 'Chest Pain Type', y = 'Exercise-induced Abnormality of the ST wave', color = 'Diagnosis', title = 'Figure 13') + theme(legend.position = "top") 


# 3.10. Inclination and Abnormalities of the Exerice Induced ST Wave:

ggplot(heartData, aes(x = fct_reorder(stSlope,stDepression), y = stDepression, fill = stSlope)) +
  geom_boxplot(show.legend = F) + xlab("Inclination for the ST segment") + labs(fill = "Diagnosis", title = 'Figure 14')+
  stat_summary(geom = "text", fun.data=function(x)  {return(c(y=median(x)*1.7,label=length(x)))}) +
  facet_grid(.~condition) + ylab("Depression for the ST segment") + geom_jitter(width=0.1, alpha = 0.1) +
  stat_summary(fun=median, geom="line", aes(group=1))  + 
  stat_summary(fun=median, geom="point")

# 3.11. Number of Vessels Colored Through Fluoroscopy:

# 3.11.1. Number of patients grouped by:

ggplot(heartData, aes(x = numVessels, fill = condition)) + geom_bar(width = 0.5) + 
  labs(fill="Diagnosis", title = 'Figure 15') + 
  scale_fill_manual(values=c("darkturquoise","salmon")) + 
  xlab("Scintigraphy results") + ylab("Number of patients")+
  theme(axis.text.x=element_text( angle=30, vjust=.8, hjust=0.8))

# 3.11.2. Probability of occurrence (Calculation):

num0 <- (heartData %>%  filter(numVessels == '0.0' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(numVessels == '0.0') %>% summarise(n()))

num1 <- (heartData %>%  filter(numVessels == '1.0' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(numVessels == '1.0') %>% summarise(n()))

num2 <- (heartData %>%  filter(numVessels == '2.0' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(numVessels == '2.0') %>% summarise(n()))

num3 <- (heartData %>%  filter(numVessels == '3.0' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(numVessels == '3.0') %>% summarise(n()))

# summary(heartData$numVessels)

vesselRate  <- cbind(rbind(num0, num1, num2, num3), levels(heartData$numVessels))
colnames(vesselRate) <- c('prob', 'vessel')

# 3.11.3. Probability of occurrence (Plot):

ggplot(data = vesselRate, aes(x = vessel, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity", width = 0.6) + labs(fill="", title = 'Figure 16') +  
  theme(axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8))  + ylab("CAD incidence rate") +
  xlab("ECG resuts") 

# 3.12. SCINTIGRAPHY

# 3.12.1. Number of patients grouped by:

ggplot(heartData, aes(x = scintigraphy, fill = condition)) + geom_bar(width = 0.5) + 
  labs(fill="Diagnosis", title = 'Figure 17') + 
  scale_fill_manual(values=c("darkturquoise","salmon")) + 
  xlab("Scintigraphy results") + ylab("Number of patients")+
  theme(axis.text.x=element_text( angle=30, vjust=.8, hjust=0.8))

# 3.12.2. Probability of occurrence (Calculation):

scinFixed <- (heartData %>%  filter(scintigraphy == 'Fixed Defect' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(scintigraphy == 'Fixed Defect') %>% summarise(n()))

scinNormal<- (heartData %>%  filter(scintigraphy == 'Normal' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(scintigraphy == 'Normal') %>% summarise(n()))

scinRevers <- (heartData %>%  filter(scintigraphy == 'Reversable Defect' & condition == 'CAD+') %>% summarise(n()))/
  (heartData %>%  filter(scintigraphy == 'Reversable Defect') %>% summarise(n()))

scinRate  <- cbind(rbind(scinFixed, scinRevers, scinNormal), levels(heartData$scintigraphy))
colnames(scinRate) <- c('prob', 'scin')

# 3.12.3. Probability of occurrence (Plot):

ggplot(data = scinRate, aes(x = scin, y = prob, fill = ..y..)) + 
  geom_bar(stat="identity", width = 0.6) + labs(fill="", title = 'Figure 18') +  
  theme(axis.text.x = element_text(angle=30, vjust=.8, hjust=0.8))  + ylab("CAD incidence rate") +
  xlab("ECG resuts") 

# 4. Data Model Test
levels(heartData$condition)

levels(heartData$condition) <- c("negativo","positivo")

# Creating training data as 80% of the dataset 

set.seed(123) 
random_sample <- createDataPartition(heartData$condition,  
                                     p = 0.75, list = FALSE) 

# Generating training dataset from the random_sample 

training_dataset  <- heartData[random_sample, ] 

# Generating testing dataset from rows which are not included in random_sample

testing_dataset <- heartData[-random_sample, ] 
testing_dataset <- as.data.frame(testing_dataset)

# Model Test

train_control <- trainControl(method = 'cv', number = 10, savePredictions = 'all', classProbs = T)
set.seed(123)
model <- train(condition ~., data = training_dataset, trControl = train_control, method = 'glm',
               family = binomial())
print(model)
summary(model)

predictions <- predict(model, newdata = testing_dataset)
confusionMatrix(predictions, testing_dataset$condition, 
                positive = 'positivo')