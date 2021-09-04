#### Project: Latin America Fertility Project ####
# This script creates ASFR and TFR for LA countries using data from IPUMS (https://international.ipums.org/international/)
#and human life table database (https://www.lifetable.de/cgi-bin/index.php). I obtained unabridged lifetables using the package MORTPAK
# Author: Juliana Jaramillo Echeverri
# Date: 20 July 2021
# Updated: 4 September 2021: To make it public for the blog.

# Contents: 

  #1. Calculate mother age at birth for each child. 
  #2. Calculate the Age Specific Fertility Rates for 14 years
  #3. Calculates TFR for 14 years
  #4. Graph outcomes
  #5. Create some maps

# INSTALL PACKAGES
library(dplyr)
require(ggplot2)
require("ipumsr")
library("viridis")
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
#### Calculations following Reid et. al. (2020) [https://www.tandfonline.com/doi/full/10.1080/00324728.2019.1630563] and using the MOMLOC in IPUMS (https://international.ipums.org/international/resources/misc_docs/pointer_working_paper_2009.pdf) ####

####  1. Calculate mother age at birth for each child using MOMLOC. ####

#Retrieve the age of the mother and store in fullMotherAgeAtBirth the value motherAge-age

peoplePerBlock <- 10000;
nextPerson <- 1;
fullMotherAgeAtBirth <- integer(nrow(fulldata))
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
#### 2. Calculate ASFR for 14 years: adjusting by children and maternal and unmatched children  ####

#for each country
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
  filter(fullMotherAgeAtBirth>0 & fullMotherAgeAtBirth<100) %>% #This is the case for Brazil and others
  select(AGE, fullMotherAgeAtBirth) %>% 
  mutate(fullMotherAgeAtcensus=AGE+fullMotherAgeAtBirth) %>% 
  group_by(fullMotherAgeAtcensus, AGE) %>% 
  summarise(n_child = n(), .groups = 'drop') %>% 
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
  ASFR_TMP$year <- fullCensus$YEAR[1] - 1 - childAge
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

#### 3. Calculate Total Fertility Rate ####
TFR <- ASFR_all %>% 
  group_by(year, COUNTRY) %>% 
  summarise(tfr=sum(asfr)) %>% 
  ungroup()

#### 4. Graph the outcome: I am separating the countries based on their geographical location or historical background ####

for (n in 1:nrow(TFR)) {
  TFR$name[n] <- code_to_country(TFR$COUNTRY[n])
}

#Argentina, Chile, Uruguay
TFR_conossur <- TFR %>% 
  filter(name=="Argentina" | name=="Chile" | name=="Uruguay")

tfr_graph_cono <- ggplot(TFR_conossur, aes(x=year, y=tfr)) +
geom_point(aes(color = name)) +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  scale_y_continuous(limits=c(2,8.5))+
  geom_smooth(aes(color = name, fill = name), method = "loess", alpha=0.1) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())

tfr_graph_cono
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/cono_sur.png", width = 7.4, height = 4.25)

#Colombia, Venezuela, Ecuador, Panama
TFR_granco <- TFR %>% 
  filter(name=="Colombia" | name=="Venezuela" | name=="Ecuador" | name=="Panama")

tfr_graph_granco <- ggplot(TFR_granco, aes(x=year, y=tfr)) +
  geom_point(aes(color = name)) +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  scale_y_continuous(limits=c(2,8.5))+
  geom_smooth(aes(color = name, fill = name), method = "loess", alpha=0.1) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())

tfr_graph_granco
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/grancol.png", width =7.4 , height =4.25  )

#Bolivia, Brazil, Paraguay
TFR_bolbra <- TFR %>% 
  filter(name=="Bolivia" | name=="Brazil" | name=="Paraguay")

tfr_graph_bolbra <- ggplot(TFR_bolbra, aes(x=year, y=tfr)) +
  geom_point(aes(color = name)) +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  scale_y_continuous(limits=c(2,8.5))+
  geom_smooth(aes(color = name, fill = name), method = "loess", alpha=0.1) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())

tfr_graph_bolbra
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/bolbra.png", width =7.4 , height = 4.25 )

#Mexico, Costa Rica, Honduras, Guatemala,
TFR_central <- TFR %>% 
  filter(name=="Mexico" | name=="Costa Rica" | name=="Honduras" | name=="Guatemala")

tfr_graph_central <- ggplot(TFR_central, aes(x=year, y=tfr)) +
  geom_point(aes(color = name)) +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  scale_y_continuous(limits=c(2,8.5))+
  geom_smooth(aes(color = name, fill = name), method = "loess", alpha=0.1) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())

tfr_graph_central
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/central.png", width =7.4 , height =4.25  )


#Trinidad y Tobago, Haiti

TFR_carib <- TFR %>% 
  filter(name=="Trinidad and Tobago" | name=="Haiti")

tfr_graph_carib <- ggplot(TFR_carib, aes(x=year, y=tfr)) +
  geom_point(aes(color = name)) +
  scale_x_continuous(breaks = seq(1956,1976, by=2)) + ylab("TFR") +  xlab("year") +
  scale_y_continuous(limits=c(2,8.5))+
  geom_smooth(aes(color = name, fill = name), method = "loess", alpha=0.1) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()  + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      axis.title =element_text(size = 20) ,
                      axis.text =element_text(size = 20),
                      plot.title = element_text(size = 15),
                      legend.text=element_text(size=20),
                      legend.position = "bottom",
                      legend.title = element_blank())

tfr_graph_carib
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/carib.png", width =7.4 , height = 4.25 )


#### 5. Map TFR in 1962 ####

TFR_1962 <- TFR %>% 
  filter(year==1962)

map.la <- map_data("world") %>% 
filter(region=="Colombia" | region=="Venezuela" | region=="Ecuador" | region=="Panama" | region=="Argentina" | region=="Chile" | region=="Uruguay" |
         region=="Bolivia" | region=="Brazil" | region=="Paraguay" | region=="Trinidad and Tobago" | region=="Haiti" | region=="Belize" | region== "El Salvador" |
         region=="Nicaragua" | region == "French Guiana" | region == "Guyana" | region=="Peru" |region=="Suriname" |
         region=="Mexico" | region=="Costa Rica" | region=="Honduras" | region=="Guatemala") %>% 
fortify

map.la <- map.la %>% 
  left_join(TFR_1962, by=c("region"="name"))

Maptfr1962 <- ggplot() +
  geom_map(data = map.la, map = map.la,
           aes(x = long, y = lat, group = group, map_id=region),
           fill = NA, colour = "#7f7f7f", size=0.1) +
  geom_map(data = map.la, map=map.la,
           aes(fill=tfr, map_id=region),
           colour="#7f7f7f", size=0.1) +
  scale_fill_continuous(type = "viridis", na.value="#BDBDBD") +
  scale_y_continuous(breaks=c()) +
  scale_x_continuous(breaks=c()) +
  labs(fill="TFR", title="Total Fertility Rate, 1962", x="", y="") +
  theme_minimal() +  theme(legend.position = "bottom", legend.key.width = unit(2.5, 'cm'))

Maptfr1962
ggsave("C:/Users/HOME/Dropbox/Webpage/LAFP/graphs/map_1962.png")
