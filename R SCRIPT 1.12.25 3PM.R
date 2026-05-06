###############################################################################
# 0. Load Packages


# install.packages("devtools")
devtools::install_github("r-lib/conflicted")

library(tidyverse)    # For data manipulation and plotting
library(lubridate)    # For date parsing
library(janitor)      # For tabyl(), clean_names(), etc.
library(psych)        # For describe()
library(car)          # For some post-hoc tests, VIF, etc.
library(nnet)         # For multinomial logistic regression if needed
library(rstatix)      # For pairwise tests, normality checks, etc.

library(dplyr)

R.version.string

if (!require(car)) {
  install.packages("car")
  library(car)
} else {
  library(car)
}

if (!require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
} else {
  library(dplyr)
}
if (!require(forcats)) {
  install.packages("forcats")
  library(forcats)
} else {
  library(forcats)
}
###############################################################################
# 1. Data Import
 
data <- read_csv("VPGH.Callers.Total.Original.csv")

# Quick overview
glimpse(data)

###############################################################################
# 2. Data Cleaning


# 2.1. Convert DATE to YEAR
# DATE format Month/Day/Year format 
data <- data %>%
  mutate(YEAR = format(as.Date(DATE, format = "%m/%d/%Y"), "%Y")) %>%
  select(-DATE)  # year doesn't work correctly.  Shows format for 2018 as, "0018"

# 2.2. Rename MARITAL STATUS to RELATIONSHIP STATUS
data <- data %>%
  rename(`RELATIONSHIP STATUS` = `MARITAL STATUS`)

# 2.3. Replace UNKN with NA
cols_to_na <- c("GENDER", "AGE", "INCOME", "EDUCATION",
                "EMPLOYMENT", "RELATIONSHIP STATUS",
                "48HR FOLLOW-UP", "1 WEEK FOLLOWUP", "ONE MONTH",
                "48HR GAMBLING STATUS", "1 WEEK GAMBLING STATUS",
                "1 MONTH GAMBLING STATUS")

for (colname in cols_to_na) {
  data[[colname]] <- ifelse(data[[colname]] == "UNKN", NA, data[[colname]])
}

# 2.4. Convert relevant variables to factors (categorical)
data <- data %>%
  mutate(
    YEAR = as.factor(YEAR),
    GENDER = as.factor(GENDER),
    AGE = as.factor(AGE),
    INCOME = as.factor(INCOME),
    EDUCATION = as.factor(EDUCATION),
    EMPLOYMENT = as.factor(EMPLOYMENT),
    `RELATIONSHIP STATUS` = as.factor(`RELATIONSHIP STATUS`),
    `48HR FOLLOW-UP` = as.factor(`48HR FOLLOW-UP`),
    `1 WEEK FOLLOWUP` = as.factor(`1 WEEK FOLLOWUP`),
    `ONE MONTH` = as.factor(`ONE MONTH`)) # was missing )
    
  # no longer needed errored out code, checked for you.
    
    data_2018 <- data %>%
      filter(YEAR == "2018") 
    
    ###############################################################################
    # 3. Descriptive Statistics (RQ1)
    
    # A) Overall Descriptive Statistics for demographic variables
    #    Do for the entire dataset, for 2018, and for 2023.
    
    # 3.1. Subset data for 2018 and 2023
    data_2018 <- data %>% filter(YEAR == "2018")
    data_2023 <- data %>% filter(YEAR == "2023")
    
    # 3.2. Summarize Age 
    tabyl(data_2018$AGE) # does not work, see line 55
    tabyl(data_2023$AGE) # does not work, see line 55
    tabyl(data$AGE) # all years
    
# you can just run code above to see the tables in-line 
  
    
 # 3.4. Call duration (continuous):
    
describe(data_2018$`CONTACT LENGTH MINUTES`) 
describe(data_2023$`CONTACT LENGTH MINUTES`) 
describe(data$`CONTACT LENGTH MINUTES`) # this works with describe(data, check = TRUE), see ?describe for rationale RE check
describe(data, check = TRUE)

class(data$`CONTACT LENGTH MINUTES`)


# Address Character:
data_2018 <- data_2018 %>%
      mutate(`CONTACT LENGTH MINUTES` = as.numeric(`CONTACT LENGTH MINUTES`)) 
    
data_2018 <- data_2018 %>%
      mutate(`CONTACT LENGTH MINUTES` = as.numeric(as.character(`CONTACT LENGTH MINUTES`))) 
    
library(readr)
data_2018 <- data_2018 %>%
      mutate(`CONTACT LENGTH MINUTES` = parse_number(`CONTACT LENGTH MINUTES`)) 
    
library(psych)
describe(data_2018$`CONTACT LENGTH MINUTES`)


# Repeat 3.4. Call duration (continuous):
describe(data_2018$`CONTACT LENGTH MINUTES`) 
describe(data_2023$`CONTACT LENGTH MINUTES`) 
data_2023 <- data_2023 %>%
      mutate(`CONTACT LENGTH MINUTES` = parse_number(`CONTACT LENGTH MINUTES`)) 
describe(data_2023$`CONTACT LENGTH MINUTES`)


# Address character error again...
data <- data %>%
  mutate(`CONTACT LENGTH MINUTES` = as.numeric(as.character(`CONTACT LENGTH MINUTES`)))


    describe(data$`CONTACT LENGTH MINUTES`)
#    vars    n    mean    sd   median  trimmed   mad   min  max  range skew   kurtosis   se
#X1   1    2677  19.46  17.58    15    16.73    13.34   0   176   176  2.35     9.38   0.34   
    data <- data %>%
      
    NAs introduced by coercion
    sum(is.na(data$`CONTACT LENGTH MINUTES`))
    
    data %>%
      filter(is.na(`CONTACT LENGTH MINUTES`)) %>%
      select(`CONTACT LENGTH MINUTES`) %>%
      distinct()
    

# Re-create subsets
    data_2018 <- data %>% filter(YEAR == "2018")
    data_2023 <- data %>% filter(YEAR == "2023")
    
 3.5. Summarize employment 
    tabyl(data_2018$EMPLOYMENT)
    tabyl(data_2023$EMPLOYMENT)
    tabyl(data$EMPLOYMENT)  # all years 
    

###################################################################################     
# 3.6. Employment (categorical):
##############################   STUCK ###########
    
describe(data_2018$`EMPLOYMENT`)
                class(data$`EMPLOYMENT`)
describe(data_2023$`EMPLOYMENT`)
describe(data$`EMPLOYMENT`)
data <- data %>%
  mutate(`EMPLOYMENT` = as.numeric(as.character(`EMPLOYMENT`)))

#___________________________________________________________________________

# 3.7. TYPE OF GAMBLING (categorical):
# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize TYPE OF GAMBLING 
class(data$`TYPE OF GAMBLING`)
tabyl(data_2018$ 'TYPE OF GAMBLING')
tabyl(data_2023$ 'TYPE OF GAMBLING')
tabyl(data$ 'TYPE OF GAMBLING')  # all years 

2018    TYPE OF GAMBLING  n     percent
                   Bingo  1 0.006134969
        CASINO, INTERNET  1 0.006134969
          CASINO/LOTTERY  2 0.012269939
        CHARITABLE/LOTTO  1 0.006134969
 CHARITABLE/LOTTO/CASINO  1 0.006134969
                   Cards  2 0.012269939
                  Casino  9 0.055214724
           Casino/Sports  1 0.006134969
                Internet  2 0.012269939
                 Lottery 46 0.282208589
Lottery/CHARITABLE/SPORTS/CASINO  1 0.006134969
          Lottery/Casino  1 0.006134969
    Lottery/STOCK MARKET  1 0.006134969
 Lottery/online gambling  1 0.006134969
           Slot Machines  2 0.012269939
Slot Machines at Casino/Track 21 0.128834356
Slot Machines at Casino/Track, Lottery,  2 0.012269939
Slot Machines at Casino/Track, Lottery, Bingo  1 0.006134969
Slot Machines at Casino/Track, Lottery, Bingo, Internet  1 0.006134969
Slot Machines at Casino/Track, lotto, sport  1 0.006134969
Slot Machines at Casino/Track/Sports  1 0.006134969
Slot machine, casino, lottery  1 0.006134969
Slot machines at Track/Casino, lottery, online  1 0.006134969
Sports  3 0.018404908
Sports/Internet  1 0.006134969
Table Games at Casino  9 0.055214724
Table games Casino/ lottery  1 0.006134969
UNKN 41 0.251533742
VLTs/Gray Machines  1 0.006134969
Video Gaming  1 0.006134969
casino, internet  1 0.006134969
casino/cards  2 0.012269939
internet  1 0.006134969
slot machines at casino/lotto  1 0.006134969

#REDUCE NOISE
data <- data %>%
  mutate(
    `TYPE OF GAMBLING'= factor_collapse (`TYPE OF GAMBLING`),
CASINO =c("CASINO, INTERNET”, “CASINO/LOTTERY”, ”Casino”, ”Casino/Sports”, ”Table Games at Casino”, ”Table games Casino/ lottery”, ”casino, internet”, “casino/cards”),
Lottery =c(“Lottery”, “Lottery/CHARITABLE/SPORTS/CASINO”, “Lottery/Casino”, “Lottery/STOCK MARKET”, “Lottery/online gambling”)
'Slot Machines' =c(“Slot Machines at Casino/Track”, “Slot Machines at Casino/Track, Lottery”, “Slot Machines at Casino/Track, Lottery, Bingo”, “Slot Machines at Casino/Track, Lottery, Bingo, Internet”, “Slot Machines at Casino/Track, lotto, sport”, “Slot Machines at Casino/Track/Sports”, “Slot machine, casino, lottery”, “Slot machines at Track/Casino, lottery, online”, “slot machines at casino/lotto”)  
Internet =c(“internet”)
Sports =c(“Sports/Internet")
)


# Corrected R script to collapse 'TYPE OF GAMBLING' categories
data <- data %>%
  mutate(
    `TYPE OF GAMBLING` = recode(`TYPE OF GAMBLING`,
                                "CASINO, INTERNET" = "CASINO",
                                "CASINO/LOTTERY" = "CASINO",
                                "Casino" = "CASINO",
                                "Casino/Sports" = "CASINO",
                                "Table Games at Casino" = "CASINO",
                                "Table games Casino/ lottery" = "CASINO",
                                "casino, internet" = "CASINO",
                                "casino/cards" = "CASINO",
                                
                                "Lottery" = "Lottery",
                                "Lottery/CHARITABLE/SPORTS/CASINO" = "Lottery",
                                "Lottery/Casino" = "Lottery",
                                "Lottery/STOCK MARKET" = "Lottery",
                                "Lottery/online gambling" = "Lottery",
                                
                                "Slot Machines at Casino/Track" = "Slot Machines",
                                "Slot Machines at Casino/Track, Lottery" = "Slot Machines",
                                "Slot Machines at Casino/Track, Lottery, Bingo" = "Slot Machines",
                                "Slot Machines at Casino/Track, Lottery, Bingo, Internet" = "Slot Machines",
                                "Slot Machines at Casino/Track, lotto, sport" = "Slot Machines",
                                "Slot Machines at Casino/Track/Sports" = "Slot Machines",
                                "Slot machine, casino, lottery" = "Slot Machines",
                                "Slot machines at Track/Casino, lottery, online" = "Slot Machines",
                                "slot machines at casino/lotto" = "Slot Machines",
                                
                                "internet" = "Internet",
                                "Sports/Internet" = "Sports",
                                
                                .default = "Other"  # Assign all other cases to "Other"
    )
  )

#________________________________________________________________________________
# Well, that messed it up.
data <- read_csv("VPGH.Callers.Total.Original.csv")

# Quick overview
glimpse(data)

# clean it up again  Data Cleaning


# 2.1. Convert DATE to YEAR
#     Assuming DATE is in Month/Day/Year format (e.g., 01/15/2018)
data <- data %>%
  mutate(YEAR = format(as.Date(DATE, format = "%m/%d/%Y"), "%Y")) %>%
  select(-DATE)  

# 2.2. Rename MARITAL STATUS to RELATIONSHIP STATUS
data <- data %>%
  rename(`RELATIONSHIP STATUS` = `MARITAL STATUS`)

# 2.3. Replace UNKN with NA in relevant columns
#     (example for multiple columns; adjust as needed)
cols_to_na <- c("GENDER", "AGE", "INCOME", "EDUCATION",
                "EMPLOYMENT", "RELATIONSHIP STATUS",
                "48HR FOLLOW-UP", "1 WEEK FOLLOWUP", "ONE MONTH",
                "48HR GAMBLING STATUS", "1 WEEK GAMBLING STATUS",
                "1 MONTH GAMBLING STATUS")

for (colname in cols_to_na) {
  data[[colname]] <- ifelse(data[[colname]] == "UNKN", NA, data[[colname]])
}

# 2.4. Convert relevant variables to factors (categorical)
data <- data %>%
  mutate(
    YEAR = as.factor(YEAR),
    GENDER = as.factor(GENDER),
    AGE = as.factor(AGE),
    INCOME = as.factor(INCOME),
    EDUCATION = as.factor(EDUCATION),
    EMPLOYMENT = as.factor(EMPLOYMENT),
    `RELATIONSHIP STATUS` = as.factor(`RELATIONSHIP STATUS`),
    `48HR FOLLOW-UP` = as.factor(`48HR FOLLOW-UP`),
    `1 WEEK FOLLOWUP` = as.factor(`1 WEEK FOLLOWUP`),
    `ONE MONTH` = as.factor(`ONE MONTH`)
    
    data <- data %>%
      mutate(
        `ONE_MONTH` = as.factor(`ONE MONTH`)
        data_2018 <- data %>% filter(YEAR == "2018")
      )
    data <- data %>%
      mutate(
        `ONE MONTH` = as.factor(`ONE MONTH`)
        # ... other mutate transformations ...
      )
    
    data_2018 <- data %>%
      filter(YEAR == "2018")
    
    
    

#Summarize TYPE OF GAMBLING 

# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize TYPE OF GAMBLING 
class(data$`TYPE OF GAMBLING`)
tabyl(data_2018$ 'TYPE OF GAMBLING')

#  data_2018$         "TYPE OF GAMBLING"  n     percent
#                                   Bingo  1 0.006134969
#                        CASINO, INTERNET  1 0.006134969
#                          CASINO/LOTTERY  2 0.012269939
#                        CHARITABLE/LOTTO  1 0.006134969
#                 CHARITABLE/LOTTO/CASINO  1 0.006134969
#                                   Cards  2 0.012269939
#                                  Casino  9 0.055214724
#                           Casino/Sports  1 0.006134969
#                                Internet  2 0.012269939
#                                 Lottery 46 0.282208589
#        Lottery/CHARITABLE/SPORTS/CASINO  1 0.006134969
#                          Lottery/Casino  1 0.006134969
#                    Lottery/STOCK MARKET  1 0.006134969
#Lottery/online gambling  1 0.006134969
#Slot Machines  2 0.012269939
#Slot Machines at Casino/Track 21 0.128834356
#Slot Machines at Casino/Track, Lottery,  2 0.012269939
#Slot Machines at Casino/Track, Lottery, Bingo  1 0.006134969
#Slot Machines at Casino/Track, Lottery, Bingo, Internet  1 0.006134969
#Slot Machines at Casino/Track, lotto, sport  1 0.006134969
#Slot Machines at Casino/Track/Sports  1 0.006134969
#Slot machine, casino, lottery  1 0.006134969
#Slot machines at Track/Casino, lottery, online  1 0.006134969
#Sports  3 0.018404908
#Sports/Internet  1 0.006134969
#Table Games at Casino  9 0.055214724
#Table games Casino/ lottery  1 0.006134969
#UNKN 41 0.251533742
#VLTs/Gray Machines  1 0.006134969
#Video Gaming  1 0.006134969
#casino, internet  1 0.006134969
#casino/cards  2 0.012269939
#internet  1 0.006134969
#slot machines at casino/lotto  1 0.006134969

tabyl(data_2023$ 'TYPE OF GAMBLING')
#> tabyl(data_2023$ 'TYPE OF GAMBLING')
#data_2023$"TYPE OF GAMBLING"   n     percent
#Bingo   2 0.004376368
#Cards   1 0.002188184
#Cards at home or with friends   3 0.006564551
#Casino   1 0.002188184
#Casino/Lottery/Mobile Betting   1 0.002188184
#Internet  22 0.048140044
#Internet- Non sports  14 0.030634573
#Lottery  37 0.080962801
##Roulette   1 0.002188184
#Slot Machines 140 0.306345733
#Sports  66 0.144420131
#Stock Market   3 0.006564551
#Table Games at Casino  33 0.072210066
#Track   1 0.002188184
#Trading Cards   1 0.002188184
#UNKN  43 0.094091904
#VLTs/Gray Machines  87 0.190371991
#Video Gaming   1 0.002188184

tabyl(data$ 'TYPE OF GAMBLING')  # all years 
#> tabyl(data$ 'TYPE OF GAMBLING')  # all years 
#data$"TYPE OF GAMBLING"   n      percent
#Bingo   9 0.0031835868
#CASINO, INTERNET   1 0.0003537319
#CASINO/LOTTERY   2 0.0007074637
#CHARITABLE/LOTTO   1 0.0003537319
#CHARITABLE/LOTTO/CASINO   1 0.0003537319
#Cards   9 0.0031835868
#Cards at Home or with Friends   2 0.0007074637
#Cards at home or with friends   3 0.0010611956
#Casino  21 0.0074283693
#Casino/Lottery   2 0.0007074637
#Casino/Lottery/Mobile Betting   1 0.0003537319
#Casino/Sports   1 0.0003537319
#Dice   1 0.0003537319
#Internet 144 0.0509373895
#Internet- Non sports  29 0.0102582243
#Internet/App Real Sports Events   2 0.0007074637
#Lottery 344 0.1216837637
#Lottery & Sports   1 0.0003537319
#Lottery, Bingo, Sports   1 0.0003537319
#Lottery, Table games at Casino   1 0.0003537319
#Lottery, slot machines   1 0.0003537319
#Lottery/CHARITABLE/SPORTS/CASINO   1 0.0003537319
#Lottery/Casino   1 0.0003537319
#Lottery/STOCK MARKET   1 0.0003537319
#Lottery/online gambling   1 0.0003537319
#Racetrack   2 0.0007074637
#Racing / Lottery   1 0.0003537319
#Racing/sports   1 0.0003537319
#Roulette   2 0.0007074637
#Skill games  38 0.0134418111
#Skill gamess/Gray Machines  27 0.0095507605
#Slot Machines 210 0.0742836930
#Slot Machines Casino/Track & Video Terminals   2 0.0007074637
#Slot Machines at Casino/Track 542 0.1917226742
#Slot Machines at Casino/Track & Video Terminals   1 0.0003537319
#Slot Machines at Casino/Track, Lottery,   2 0.0007074637
#Slot Machines at Casino/Track, Lottery, Bingo   1 0.0003537319
#Slot Machines at Casino/Track, Lottery, Bingo, Internet   1 0.0003537319
#Slot Machines at Casino/Track, lotto, sport   1 0.0003537319
#Slot Machines at Casino/Track/Sports   1 0.0003537319
#Slot machine, casino, lottery   1 0.0003537319
#Slot machines at Casino/Track - Video Terminals   3 0.0010611956
#Slot machines at Casino/Track, Bingo   1 0.0003537319
#Slot machines at Casino/Track, Lottery   4 0.0014149275
#Slot machines at Casino/Track, Sports   1 0.0003537319
#Slot machines at Track/Casino, lottery, online   1 0.0003537319
#Slot machines, Table Games at Casino   1 0.0003537319
#Slot machines, Table Games at Casino/Track   2 0.0007074637
#Sports 336 0.1188539087
#Sports Betting at Casino   1 0.0003537319
#Sports/Internet   1 0.0003537319
#Stock Market  26 0.0091970287
#Table Games at Casino 230 0.0813583304
#Table Games at Casino, Sports   1 0.0003537319
#Table games Casino/ lottery   1 0.0003537319
#Table games at Casino  10 0.0035373187
#Table games at Casino / Day Trading   1 0.0003537319
#Tip Tickets   1 0.0003537319
#Track  23 0.0081358330
#Track, Sports   1 0.0003537319
#Trading Cards   1 0.0003537319
#UNKN 370 0.1308807924
#VLTs/Gray Machines 100 0.0353731871
#Video Gaming  16 0.0056597099
#Video Terminals 276 0.0976299965
#casino, internet   1 0.0003537319
#casino/cards   2 0.0007074637
#internet   1 0.0003537319
#slot machines at casino/lotto   1 0.0003537319


#Summarize GENDER 

# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize GENDER 
class(data$`GENDER`)
tabyl(data_2018$GENDER)   
tabyl(data_2023$GENDER)   
tabyl(data$GENDER)   

#> tabyl(data_2018$GENDER) 
#data_2018$GENDER  n    percent valid_percent
#F 60 0.36809816     0.3947368
#M 92 0.56441718     0.6052632
#<NA> 11 0.06748466            NA
#> tabyl(data_2023$GENDER) 
#data_2023$GENDER   n    percent valid_percent
#F 163 0.35667396     0.3808411
#M 265 0.57986871     0.6191589
#<NA>  29 0.06345733            NA
#> tabyl(data$GENDER)   
#data$GENDER    n    percent valid_percent
#F  913 0.32295720      0.350211
#M 1694 0.59922179      0.649789
#<NA>  220 0.07782101            NA

#Summarize GENDER 

# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize CALLER 
#class(data$`CALLER`)
#tabyl(data_2018$CALLER)   
#tabyl(data_2023$CALLER)   
#tabyl(data$CALLER)  
#> tabyl(data_2018$CALLER)   
#data_2018$CALLER   n    percent
#Family/Friend  21 0.12883436
#Other   2 0.01226994
#Self 123 0.75460123
#Spouse/SO  15 0.09202454
#UNKN   2 0.01226994
#> tabyl(data_2023$CALLER)   
#data_2023$CALLER   n     percent
#Family/Friend  46 0.100656455
#Other   3 0.006564551
#Self 347 0.759299781
#Spouse/SO  24 0.052516411
#UNKN  37 0.080962801
#> tabyl(data$CALLER) 
#data$CALLER    n     percent
#Family/Friend  423 0.149628582
#Other   13 0.004598514
#Self 2021 0.714892112
#Spouse/SO  161 0.056950831
#UNKN  209 0.073929961

# 3.1. Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

  # 3.2. Summarize RELATIONSHIP STATUS 
class(data$`RELATIONSHIP STATUS`)
tabyl(data_2018$`RELATIONSHIP STATUS`)
tabyl(data_2023$'RELATIONSHIP STATUS')
tabyl(data$'RELATIONSHIP STATUS')  # all years

#> tabyl(data_2018$`RELATIONSHIP STATUS`)
#data_2018$`RELATIONSHIP STATUS`  n     percent valid_percent
#Divorced  7 0.042944785    0.10606061
#Married 40 0.245398773    0.60606061
#Remarried  1 0.006134969    0.01515152
#Single 17 0.104294479    0.25757576
#Widowed  1 0.006134969    0.01515152
#<NA> 97 0.595092025            NA
#> tabyl(data_2023$'RELATIONSHIP STATUS')
#data_2023$"RELATIONSHIP STATUS"   n     percent valid_percent
#Divorced  18 0.039387309   0.060200669
#Living Together/Engaged  45 0.098468271   0.150501672
#Married 107 0.234135667   0.357859532
#Separated  10 0.021881838   0.033444816
#Seperated   1 0.002188184   0.003344482
#Single 105 0.229759300   0.351170569
#Widowed  13 0.028446389   0.043478261
#<NA> 158 0.345733042            NA
#> tabyl(data$'RELATIONSHIP STATUS')  # all years
#data$"RELATIONSHIP STATUS"    n      percent valid_percent
#Divorced  143 0.0505836576   0.085577499
#Living Together/Engaged   89 0.0314821365   0.053261520
#Married  671 0.2373540856   0.401555955
#Remarried  127 0.0449239476   0.076002394
#Separated   21 0.0074283693   0.012567325
#Seperated    2 0.0007074637   0.001196888
#Single  573 0.2026883622   0.342908438
#Widowed   45 0.0159179342   0.026929982
#<NA> 1156 0.4089140432            NA

data <- data %>%
  mutate(
    `RELATIONSHIP STATUS` = factor_collapse(`RELATIONSHIP STATUS`,
                                            Separated = c("Separated", "Seperated") 
    )
                                            
    class(data$`RELATIONSHIP STATUS`)
    tabyl(data_2018$`RELATIONSHIP STATUS`)
    tabyl(data_2023$'RELATIONSHIP STATUS')
    tabyl(data$'RELATIONSHIP STATUS')  # all years
    
#    >     tabyl(data_2018$`RELATIONSHIP STATUS`)
#    data_2018$`RELATIONSHIP STATUS`  n     percent valid_percent
#    Divorced  7 0.042944785    0.10606061
#    Married 40 0.245398773    0.60606061
#    Remarried  1 0.006134969    0.01515152
#    Single 17 0.104294479    0.25757576
#    Widowed  1 0.006134969    0.01515152
#    <NA> 97 0.595092025            NA
#    >     tabyl(data_2023$'RELATIONSHIP STATUS')
#    data_2023$"RELATIONSHIP STATUS"   n     percent valid_percent
#    Divorced  18 0.039387309   0.060200669
#    Living Together/Engaged  45 0.098468271   0.150501672
#    Married 107 0.234135667   0.357859532
#    Separated  10 0.021881838   0.033444816
#    Seperated   1 0.002188184   0.003344482
#    Single 105 0.229759300   0.351170569
#    Widowed  13 0.028446389   0.043478261
#    <NA> 158 0.345733042            NA
 #   >     tabyl(data$'RELATIONSHIP STATUS')
#    data$"RELATIONSHIP STATUS"    n      percent valid_percent
#    Divorced  143 0.0505836576   0.085577499
#    Living Together/Engaged   89 0.0314821365   0.053261520
#    Married  671 0.2373540856   0.401555955
#    Remarried  127 0.0449239476   0.076002394
#    Separated   21 0.0074283693   0.012567325
#    Seperated    2 0.0007074637   0.001196888
#    Single  573 0.2026883622   0.342908438
#    Widowed   45 0.0159179342   0.026929982
#    <NA> 1156 0.4089140432            NA
    
    # 3.1. Subset data for 2018 and 2023
    data_2018 <- data %>% filter(YEAR == "2018")
    data_2023 <- data %>% filter(YEAR == "2023")
    
    # 3.2. Summarize INCOME 
    class(data$INCOME)
    tabyl(data_2018$INCOME)
    tabyl(data_2023$INCOME)
    tabyl(data$INCOME)  # all years
    
#__________________________________________________
#Income would benefit from cleaning/noise reduction.   
 
 
#    >    class(data$INCOME)
#    [1] "character"
#    >     tabyl(data_2018$INCOME)
#    data_2018$INCOME   n     percent valid_percent
#    0-10K   2 0.012269939    0.11764706
#    0-10k   1 0.006134969    0.05882353
#    100+   1 0.006134969    0.05882353
#    101+   4 0.024539877    0.23529412
#    11-20K   1 0.006134969    0.05882353
#    11-20k   2 0.012269939    0.11764706
#    31-40   1 0.006134969    0.05882353
#    41-50k   1 0.006134969    0.05882353
#    51-100   1 0.006134969    0.05882353
#    51-100k   2 0.012269939    0.11764706
#    D   1 0.006134969    0.05882353
#    <NA> 146 0.895705521            NA
#    >     tabyl(data_2023$INCOME)
#    data_2023$INCOME   n     percent valid_percent
#    0-10k yr   2 0.004376368    0.02666667
#    101K + yr   1 0.002188184    0.01333333
#    101K+ yr   2 0.004376368    0.02666667
#    101k + yr   4 0.008752735    0.05333333
#    101k+ yr   6 0.013129103    0.08000000
#    21-30k yr   5 0.010940919    0.06666667
#    21k-30k yr   4 0.008752735    0.05333333
#    31-40k yr   5 0.010940919    0.06666667
#    31-40kyr   1 0.002188184    0.01333333
#    31k-40k yr   1 0.002188184    0.01333333
#    41-50k yr  13 0.028446389    0.17333333
 #   41k-50k yr   5 0.010940919    0.06666667
#    51-100k yr  10 0.021881838    0.13333333
#    51k-100k yr  15 0.032822757    0.20000000
#    Disability   1 0.002188184    0.01333333
#    <NA> 382 0.835886214            NA
#    >     tabyl(data$INCOME)  # all years
#    data$INCOME    n      percent valid_percent
#    0-10K    3 0.0010611956   0.022058824
#    0-10k    1 0.0003537319   0.007352941
#    0-10k yr    4 0.0014149275   0.029411765
#    100+    2 0.0007074637   0.014705882
#    101 +    1 0.0003537319   0.007352941
#    101+    9 0.0031835868   0.066176471
#    101K + yr    1 0.0003537319   0.007352941
#    101K+ yr    2 0.0007074637   0.014705882
#    101k + yr    5 0.0017686594   0.036764706
#    101k+ yr    8 0.0028298550   0.058823529
#    11-20K    1 0.0003537319   0.007352941
#    11-20k    4 0.0014149275   0.029411765
#    21-30    1 0.0003537319   0.007352941
#    21-30K    3 0.0010611956   0.022058824
#    21-30k yr    5 0.0017686594   0.036764706
#    21k-30k yr    4 0.0014149275   0.029411765
#    31-40    1 0.0003537319   0.007352941
#    31-40k yr    8 0.0028298550   0.058823529
#    31-40kyr    1 0.0003537319   0.007352941
#    31k-40k yr    1 0.0003537319   0.007352941
#    41-50    2 0.0007074637   0.014705882
#    41-50K    1 0.0003537319   0.007352941
#    41-50K/yr    2 0.0007074637   0.014705882
#    41-50k    1 0.0003537319   0.007352941
#    41-50k yr   18 0.0063671737   0.132352941
#    41k-50k yr    5 0.0017686594   0.036764706
#    51-100    1 0.0003537319   0.007352941
#    51-100K    1 0.0003537319   0.007352941
#    51-100k    3 0.0010611956   0.022058824
#    51-100k yr   20 0.0070746374   0.147058824
#    51k-100k yr   15 0.0053059781   0.110294118
#    D    1 0.0003537319   0.007352941
#    Disability    1 0.0003537319   0.007352941
#    <NA> 2691 0.9518924655            NA

    
    
    
 # Save all datasets into a single RData file
    save(data_2018, data_2023, data, file = "all_data_backup.RData")
    
 # To load the datasets back into R later, use:
  # load("all_data_backup.RData")
    # This will load data_2018, data_2023, and data into your environment
    
    #Just in case:
    # CSV file
    write.csv(data_2018, file = "data_2018_backup.csv", row.names = FALSE)
    write.csv(data_2023, file = "data_2023_backup.csv", row.names = FALSE)
    write.csv(data, file = "data_backup.csv", row.names = FALSE)
    
    # To load the datasets back into R later, use:
    # data_2018 <- read.csv("data_2018_backup.csv", stringsAsFactors = FALSE)
    # data_2023 <- read.csv("data_2023_backup.csv", stringsAsFactors = FALSE)
    # data <- read.csv("data_backup.csv", stringsAsFactors = FALSE)
    
    