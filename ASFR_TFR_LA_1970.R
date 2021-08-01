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
printf <- function(...) cat(sprintf(...))

#set pathway
setwd("C:/Users/HOME/Dropbox/Webpage/LAFP/data")

#data from ipums: Contains all LA countries that have data for the 1970 period.
# NOTE: To load data, you must download both the extract's data and the DDI
# and also set the working directory to the folder with these files (or change the path).

if (!exists("fulldata")) {
  ddi <- read_ipums_ddi("ipumsi_00028.xml")
  fulldata <- read_ipums_micro(ddi)
  fulldata$household_id <- paste(fulldata$SAMPLE, fulldata$SERIAL, sep="")
  fulldata$SAMPLE <- NULL
  fulldata$SERIAL <- NULL
}
#### Calculate the ASFR using the OCM from Reid et. al. and using the MOMLOC in ipums developed by ####

#1. Calculate mother age at birth for each child using MOMLOC. 

#Retrieve the age of the mother and store in fullMotherAgeAtBirth the value motherAge-age

##profiling variables
peoplePerBlock <- 10000;
nextPerson <- 1;
fullMotherAgeAtBirth <- integer(nrow(fulldata))

profilingPeople = 100000;
maxLoads = profilingPeople/peoplePerBlock;
numLoads = 0;
profiling = 0;
households <- unique(fulldata$household_id)
householdsOriginal <- households
householdTotal <- length(households)
householdCount <- 0;
timestart = proc.time()

while(nextPerson <= nrow(fulldata)) {
  firstPerson <- nextPerson;
  census <- fulldata[firstPerson:min(nrow(fulldata), nextPerson + peoplePerBlock - 1), ]
  census <- apply(as.matrix.noquote(census), 2, as.numeric)
  lastHid <- census[nrow(census), "household_id"]
  
  lastPerson <- firstPerson + nrow(census) - 1;
  nextPerson <- lastPerson + 1;
  
  theseHouseholds <- unique(census[, "household_id"]);
  
  for(hid in theseHouseholds) {
    household_census_rows = (census[, "household_id"] == hid)
    indeces <- which(household_census_rows)
    household = census[household_census_rows,,drop=FALSE]
    # for each potential child, retrieve the age and store in fullMotherAgeAtBirth the value motherAge-age
    for (p in 1:nrow(household)) {
      person2 <- household[p,];
      if(person2["MOMLOC"] != 0) {
        # this person is a child
        mother <- household[household[,"PERNUM"] == person2["MOMLOC"],,drop=FALSE];
        if(nrow(mother) != 1){ 
          next;
        }
        motherAge = mother[,"AGE"];
        childAge <- household[p, "AGE"]
        fullIdx = indeces[p] + firstPerson - 1
        fullMotherAgeAtBirth[fullIdx] <- motherAge - childAge
      }
    }
    if(householdCount %% 1000 == 0) {
      timecur <- proc.time()
      elapsed <- as.numeric(timecur[1] - timestart[1] + timecur[2] - timestart[2])/3600
      printf("Done %d households out of %d (%.2f%%) [approx %d people out of %d] in %.3fhrs (projection: %.2fh)\n", 
             householdCount, householdTotal, householdCount/householdTotal * 100, firstPerson, nrow(fulldata), elapsed, elapsed*householdTotal/householdCount)
    }
    householdCount <- householdCount + 1
  } # for hid
}

fulldata$fullMotherAgeAtBirth <- fullMotherAgeAtBirth
fulldata$household_id <- NULL

code_to_country <- function(code) {
  labels_ <- ipums_val_labels(fulldata$COUNTRY)  
  for (n_ in 1:nrow(labels_)) {
    if(labels_$val[n_] == code)
      return (labels_$lbl[n_]);
  }
  return ("ERROR")
}
#2. Calculate ASFR for 14 years: adjusting by children and maternal and unmatched children 


#for each country()
ASFR_YEARS_total <- NULL
theseCountries <- unique(fulldata$COUNTRY)

for (c in 1:length(theseCountries)) {
asfr <- matrix(nrow=41, ncol = 3)
asfr <- as.data.frame(asfr)
colnames(asfr) <- c("numerator", "denominator", "asfr")
asfr$age <- c(15:55)

code <- theseCountries[c]
country <- code_to_country(code)

#Import lifetables data by country
probsurv_child <- read.csv(sprintf("%s_ProbSurv.csv", country))
probsurv_women <- read.csv(sprintf("%s_WomenProbSurv.csv", country))

fullCensus <- fulldata %>%       # fullCensus is now the full sample of one country
  filter(COUNTRY==code)

#Adjusting by unmatched children

unmatched <- fullCensus %>% 
  filter(AGE<=14 & fullMotherAgeAtBirth==0) %>% 
  group_by(AGE, COUNTRY) %>% 
  summarise(n_child_unmatch=n())

matched <- fullCensus %>% 
  filter(AGE<=14 & fullMotherAgeAtBirth!=0) %>% 
  group_by(AGE, COUNTRY) %>% 
  summarise(n_child_match=n(), .groups = "keep")


for (i in 1:15) {
  unmatched$prop[i] <- matched$n_child_match[i]/sum(fullCensus$AGE==i-1)
}

allWomenAge <- fullCensus %>% 
  filter(SEX==2) %>% 
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
  ASFR_TMP$year <- fullCensus$YEAR[1] - childAge
  ASFR_YEARS <- rbind(ASFR_YEARS, ASFR_TMP)
  rm(ASFR_TMP)
}

ASFR_YEARS$Ageatbirth <- ASFR_YEARS$fullMotherAgeAtcensus - ASFR_YEARS$AGE
ASFR_YEARS$COUNTRY <- code
if (exists("ASFR_YEARS_total"))
  ASFR_YEARS_total <- rbind(ASFR_YEARS_total, ASFR_YEARS)
else
  ASFR_YEARS_total <- ASFR_YEARS
} # for each country

ASFR_all <- ASFR_YEARS_total %>% 
  group_by(fullMotherAgeAtcensus, AGE, year, COUNTRY) %>% 
  summarise(total_women=sum(women_adjusted),
            total_child=sum(n_child_adjusted), .groups = 'drop') %>% 
  mutate(asfr=total_child/total_women)

#3. Calculate Total Fertility Rate
TFR <- ASFR_all %>% 
  group_by(year, COUNTRY) %>% 
  summarise(tfr=sum(asfr)) %>% 
  ungroup()

#4. Graph the outcome: I am separating the countries based on their geographical location or historical background

for (n in 1:nrow(TFR)) {
  TFR$name[n] <- code_to_country(TFR$COUNTRY[n])
}

#Argentina, Chile, Uruguay
TFR_conossur <- TFR %>% 
  filter(COUNTRY=="Argentina" | COUNTRY=="Chile" | COUNTRY=="Uruguay")

#Colombia, Venezuela, Ecuador, Panama

#Bolivia, Brazil, Paraguay

#Mexico, Costa Rica, Honduras, Guatemala,

#Trinidad y Tobago, Haiti


tfr_graph <- ggplot(TFR, aes(x=year, y=tfr, group=as.factor(COUNTRY), color=as.factor(COUNTRY))) + 
  geom_smooth(method = "loess", alpha=0, size=1.5) + 
  #geom_line() +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  #scale_color_brewer(palette = "Dark2", type = "qual") +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())


tfr_graph

ggsave("TFR_groups.png")

