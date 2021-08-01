#### Project: Latin America Fertility Project ####
# This script creates ASFR and TFR for LA countries using data from IPUMS (https://international.ipums.org/international/)
#and human life table database (https://www.lifetable.de/cgi-bin/index.php). I obtained unabridged lifetables using the package MORTPAK
# Author: Juliana Jaramillo Echeverri
# Date: 20 July 2021
# Updated: 

# Contents: 
  #1. Calculate mother age at birth for each child. 
  #2. Calculate the Age Specific Fertility Rates

# INSTALL PACKAGES
library(dplyr)
require(ggplot2)
require("ipumsr")

#set pathway
setwd("C:/Users/HOME/Dropbox/Webpage/LAFP/data")

#data from ipums: Contains all LA countries that have data for the 1970 period.
# NOTE: To load data, you must download both the extract's data and the DDI
# and also set the working directory to the folder with these files (or change the path).

ddi <- read_ipums_ddi("ipumsi_00028.xml")
fulldata <- read_ipums_micro(ddi)

fulldata$household_id <- paste(fulldata$SAMPLE, fulldata$SERIAL, sep="")

#### Calculate the ASFR using the OCM from Reid et. al. and using the MOMLOC in ipums developed by ####

#1. Calculate mother age at birth for each child using MOMLOC. 

#Retrieve the age of the mother and store in fullMotherAgeAtBirth the value motherAge-age

##profiling variables
timestart = proc.time()
peoplePerBlock <- 10000;
nextPerson <- 1;
fullMotherAgeAtBirth <- integer(nrow(fulldata))

profilingPeople = 100000;
maxLoads = profilingPeople/peoplePerBlock;
numLoads = 0;
profiling = 0;
while(nextPerson <= nrow(fulldata)) {
  firstPerson <- nextPerson;
  census <- fulldata[firstPerson:min(nrow(fulldata), nextPerson + peoplePerBlock - 1), ]
  
  lastHid <- census$household_id[nrow(census)]
  census[census$household_id != lastHid, ]
  
  lastPerson <- firstPerson + nrow(census) - 1;
  nextPerson <- lastPerson + 1;
  
  theseHouseholds <- unique(census$household_id);
  
  for(hid in theseHouseholds) {
    household_census_rows = (census$household_id == hid)
    indeces <- which(household_census_rows)
    household = census[household_census_rows, ]
    # for each potential child, retrieve the age and store in fullMotherAgeAtBirth the value motherAge-age
    for (p in 1:nrow(household)) {
      person2 <- household[p,];
      if(person2$MOMLOC != 0) {
        # this person is a child
        mother <- subset(household, PERNUM == person2$MOMLOC);
        if(nrow(mother) != 1){ 
          next;
        }
        motherAge = mother$AGE;
        childAge <- household[[p, "AGE"]]
        fullIdx = indeces[p] + firstPerson - 1
        fullMotherAgeAtBirth[fullIdx] <- motherAge - childAge
      }
    }
  }
}

code_to_country <- function(code) {
  labels_ <- ipums_val_labels(fulldata$COUNTRY)  
  for (n_ in 1:nrow(labels_)) {
    if(labels_$val[n_] == code)
      return (labels_$lbl[n_]);
  }
  return ("ERROR")
}
#2. Calculate ASFR for 14 years: adjusting by children and maternal and unmatched children 
#Adjusting by unmatched children

unmatched <- fulldata %>% 
  filter(AGE<=14 & fullMotherAgeAtBirth==0) %>% 
  group_by(AGE, COUNTRY) %>% 
  summarise(n_child_unmatch=n())

matched <- fulldata %>% 
  filter(AGE<=14 & fullMotherAgeAtBirth!=0) %>% 
  group_by(AGE, COUNTRY) %>% 
  summarise(n_child_match=n())


for (i in 1:15) {
  unmatched$prop[i] <- matched$n_child_match[i]/sum(fullCensus$AGE==i-1)
}

#for each country()

asfr <- matrix(nrow=41, ncol = 3)
asfr <- as.data.frame(asfr)
colnames(asfr) <- c("numerator", "denominator", "asfr")
asfr$age <- c(15:55)

probsurv_child <- read.csv("country_ProbSurv.csv")
probsurv_women <- read.csv("country_WomenProbSurv.csv")

allWomenAge <- fullCensus %>% 
  filter(SEX==2) %>% 
  filter() %>%  #by country
  count(AGE)

#Calculating the ASFR per each year
ASFR_TFR <- fullCensus %>% 
  filter(fullMotherAgeAtBirth!=0) %>% 
  select(AGE, fullMotherAgeAtBirth) %>% 
  mutate(fullMotherAgeAtcensus=AGE+fullMotherAgeAtBirth) %>% 
  group_by(fullMotherAgeAtcensus, AGE) %>% 
  summarise(n_child = n(), .groups = 'drop') %>% 
  #spread(AGE, n_child)
  ungroup()

ASFR_TFR$n_child_adjusted <- NA
ASFR_YEARS <- as.data.frame(matrix(0, nrow = 0, ncol = ncol(ASFR_TFR)+5))
colnames(ASFR_YEARS) <- c(colnames(ASFR_TFR), 'n_child_adjusted', 'women', 'asfr', 'year', 'prob')

motherAgeAtChildbirth <- ASFR_TFR$fullMotherAgeAtcensus - ASFR_TFR$AGE;
for(n in 1:nrow(ASFR_TFR)) {
  ASFR_TFR$prob[n] <- probsurv_women$l.x[probsurv_women$Age == ASFR_TFR$fullMotherAgeAtcensus[n]] /
    probsurv_women$l.x[probsurv_women$Age == motherAgeAtChildbirth[n]]
}

for (childAge in 0:14) {
  motherAgeMin = childAge + 15
  motherAgeMax = childAge + 50
  
  ASFR_TFR$n_child_adjusted[ASFR_TFR$AGE==childAge] <- ASFR_TFR$n_child[ASFR_TFR$AGE==childAge]/unmatched$prop[unmatched$AGE==childAge]
  ASFR_TFR$n_child_adjusted[ASFR_TFR$AGE==childAge] <- ASFR_TFR$n_child_adjusted[ASFR_TFR$AGE==childAge]/probsurv_child$ProbSurv[probsurv_child$Age==childAge]
  
  ASFR_TMP <- ASFR_TFR %>%
    filter(AGE==childAge) %>%
    filter(fullMotherAgeAtcensus>=motherAgeMin & fullMotherAgeAtcensus<=motherAgeMax)
  
  ASFR_TMP <- ASFR_TMP %>% 
    mutate(women=allWomenAge$n[allWomenAge$AGE %in% c(ASFR_TMP$fullMotherAgeAtcensus)],
           women_adjusted= women/prob,
           asfr=n_child_adjusted/women_adjusted)
  ASFR_TMP$year <- 1972 - childAge
  ASFR_YEARS <- rbind(ASFR_YEARS, ASFR_TMP)
  rm(ASFR_TMP)
}

ASFR_YEARS$Ageatbirth <- ASFR_YEARS$fullMotherAgeAtcensus - ASFR_YEARS$AGE

#ASFR_antlow <- ASFR_YEARS
#ASFR_boy <- ASFR_YEARS
#ASFR_cau <- ASFR_YEARS
#ASFR_cund <- ASFR_YEARS
ASFR_nari <- ASFR_YEARS


ASFR_YEARS_total <- rbind(ASFR_boy, ASFR_cau)
ASFR_YEARS_total <- rbind(ASFR_YEARS_total, ASFR_cund)
ASFR_YEARS_total <- rbind(ASFR_YEARS_total, ASFR_nari)
ASFR_YEARS_total <- rbind(ASFR_YEARS_total, ASFR_antlow)

ASFR_all <- ASFR_YEARS_total %>% 
  group_by(fullMotherAgeAtcensus, AGE, year) %>% 
  summarise(total_women=sum(women_adjusted),
            total_child=sum(n_child_adjusted), .groups = 'drop') %>% 
  mutate(asfr=total_child/total_women)

#3. Calculate Total Fertility Rate
TFR <- ASFR_all %>% 
  group_by(year) %>% 
  summarise(tfr=sum(asfr))

#4. Graph the outcome
tfr_graph <- ggplot(TFR, aes(x=year, y=tfr, group=as.factor(Census), color=as.factor(Census))) + 
  geom_smooth(method = "loess", alpha=0, size=1.5) + 
  #geom_line() +
  scale_x_continuous(breaks = seq(1958,1970, by=2)) + ylab("TFR") +  xlab("year") +
  scale_color_brewer(palette = "Dark2", type = "qual") +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())


tfr_graph

ggsave("TFR_groups.png")
