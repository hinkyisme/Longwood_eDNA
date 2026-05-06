###############################################################################
#Packages

# install.packages("devtools")


library(tidyverse)    # For data manipulation and plotting
library(lubridate)    # For date parsing
library(janitor)      # For tabyl(), clean_names(), etc.
library(psych)        # For describe()
library(car)          # For some post-hoc tests, VIF, etc.
library(nnet)         # For multinomial logistic regression if needed
library(rstatix)      # For pairwise tests, normality checks, etc.

library(dplyr)
devtools::install_github("r-lib/conflicted")


R.version.string
library (forcats)


if (!require(car)) {
  install.packages("car")
  library(car)
} else {
  library(car)
}


#______________________________________________________________________
# 1. Data 

data <- read_csv("VPGH.Callers.Total.Original.csv")

# Quick overview
glimpse(data)

#__________________________________________________________________

# 2. Data Cleaning


# 2.1. Convert DATE to YEAR
# DATE format Month/Day/Year format 
data <- data %>%
  mutate(YEAR = format(as.Date(DATE, format = "%m/%d/%Y"), "%Y")) %>%
  select(-DATE)  

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
    `ONE MONTH` = as.factor(`ONE MONTH`)
    
  )
    
    data_2018 <- data %>%
      filter(YEAR == "2018")
    
#____________________________________________________________
# 3. Descriptive Statistics (RQ1)
    
# A) Overall Descriptive Statistics for demographic variables
#    Entire dataset, for 2018, and for 2023.
    
# 3.1. Subset data for 2018 and 2023
    data_2018 <- data %>% filter(YEAR == "2018")
    data_2023 <- data %>% filter(YEAR == "2023")
    
    # 3.2. Summarize Age 
    tabyl(data_2018$AGE)
    tabyl(data_2023$AGE)
    tabyl(data$AGE)  # all years
    
# 3.4. Call duration (continuous):
    
class(data$`CONTACT LENGTH MINUTES`)
    
    
# Address Character:
library(readr)
library(psych)

data <- data %>%
  mutate(`CONTACT LENGTH MINUTES` = as.numeric(as.character(`CONTACT LENGTH MINUTES`)))
data %>%
  filter(is.na(`CONTACT LENGTH MINUTES`)) %>%
  select(`CONTACT LENGTH MINUTES`) %>%
  distinct()

data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

describe(data$`CONTACT LENGTH MINUTES`)
describe(data_2018$`CONTACT LENGTH MINUTES`)
describe(data_2023$`CONTACT LENGTH MINUTES`)
 
   
#Employment
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

tabyl(data_2018$EMPLOYMENT)
tabyl(data_2023$EMPLOYMENT)
tabyl(data$EMPLOYMENT)  # all years

glimpse(data)
library(readr)
library(psych)

summary(data$EMPLOYMENT)

#____________________________________________________________


# 3.7. TYPE OF GAMBLING (categorical):
# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize TYPE OF GAMBLING 
class(data$`TYPE OF GAMBLING`)
tabyl(data_2018$ 'TYPE OF GAMBLING')
tabyl(data_2023$ 'TYPE OF GAMBLING')
tabyl(data$ 'TYPE OF GAMBLING')  # all years 


# Quick overview
glimpse(data)

#Create new colomn to keep raw and clean gambling type data
library(dplyr)
library(stringr)  # for str_detect()

glimpse(data)

janitor::tabyl(data$TYPE_GAMBLING_CLEAN)
class(data$TYPE_GAMBLING_CLEAN)
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

tabyl(data_2018$ TYPE_GAMBLING_CLEAN)
tabyl(data_2023$ TYPE_GAMBLING_CLEAN)
tabyl(data$ TYPE_GAMBLING_CLEAN)  # all years 

#Summarize GENDER#_____________________________________# 

# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#Summarize GENDER 
class(data$`GENDER`)
tabyl(data_2018$GENDER)   
tabyl(data_2023$GENDER)   
tabyl(data$GENDER)   


# Re-create subsets
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

#___________________________________________________________
# 3. Summarize RELATIONSHIP STATUS 
class(data$`RELATIONSHIP STATUS`)
tabyl(data_2018$`RELATIONSHIP STATUS`)
tabyl(data_2023$'RELATIONSHIP STATUS')
tabyl(data$'RELATIONSHIP STATUS')  # all years

library(dplyr)
library(stringr)  # for str_detect, if needed

data <- data %>%
  mutate(
    RELATIONSHIP_STATUS_CLEAN = case_when(
      # Fix the spelling of "Seperated" to "Separated"
      `RELATIONSHIP STATUS` == "Seperated" ~ "Separated",
      
          # Keep other categories as they are:
      TRUE ~ `RELATIONSHIP STATUS`
    )
  )

#Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

class(data$RELATIONSHIP_STATUS_CLEAN)
tabyl(data_2018$RELATIONSHIP_STATUS_CLEAN)
tabyl(data_2023$RELATIONSHIP_STATUS_CLEAN)
tabyl(data$RELATIONSHIP_STATUS_CLEAN)  # all years

#_______________________________________________________
#4. INCOME

#Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

class(data$INCOME)
tabyl(data_2018$INCOME)
tabyl(data_2023$INCOME)
tabyl(data$INCOME)  # all years

library(dplyr)
library(stringr)



data <- data %>%
  mutate(
    INCOME_CLEAN = case_when(
      # 0–10k forms
      str_detect(tolower(INCOME), "0-10") ~ "0-10k",
      str_detect(tolower(INCOME), "0-10k") ~ "0-10k",
      
      # 11–20k forms
      str_detect(tolower(INCOME), "11-20") ~ "11-20k",
      
      # 21–30k forms
      str_detect(tolower(INCOME), "21-30") ~ "21-30k",
      str_detect(tolower(INCOME), "21-30K") ~ "21-30k",
      str_detect(tolower(INCOME), "21-30k") ~ "21-30k",
      str_detect(tolower(INCOME), "21-30k yr") ~ "21-30k",
      
      # 31–40k forms
      str_detect(tolower(INCOME), "31-40") ~ "31-40k",
      str_detect(tolower(INCOME), "31-40k yr") ~ "31-40k",
      str_detect(tolower(INCOME), "31-40kyr") ~ "31-40k",
      str_detect(tolower(INCOME), "31k-40k yr") ~ "31-40k",
      
      # 41–50k forms
      str_detect(tolower(INCOME), "41-50") ~ "41-50k",
      str_detect(tolower(INCOME), "41-50K") ~ "41-50k",
      str_detect(tolower(INCOME), "41-50k/yr") ~ "41-50k",
      str_detect(tolower(INCOME), "41-50k") ~ "41-50k",
      str_detect(tolower(INCOME), "41k-50k yr") ~ "41-50k",
      str_detect(tolower(INCOME), "41k-50k yr") ~ "41-50k",
      
      # 51–100k forms
      str_detect(tolower(INCOME), "51-100") ~ "51-100k",
      str_detect(tolower(INCOME), "51-100K") ~ "51-100k",
      str_detect(tolower(INCOME), "51-100k") ~ "51-100k",
      str_detect(tolower(INCOME), "51-100k yr") ~ "51-100k",
      str_detect(tolower(INCOME), "51k-100k yr") ~ "51-100k",
      
      # 101k+ forms
      str_detect(tolower(INCOME), "100+") ~ "101k+",
      str_detect(tolower(INCOME), "101") ~ "101k+",
      str_detect(tolower(INCOME), "101+") ~ "101k+",
      str_detect(tolower(INCOME), "101 +") ~ "101k+",
      str_detect(tolower(INCOME), "101K + yr") ~ "101k+",
      str_detect(tolower(INCOME), "101K+yr") ~ "101k+",
      str_detect(tolower(INCOME), "101k+yr") ~ "101k+",
      str_detect(tolower(INCOME), "101k + yr") ~ "101k+",
      # "D" or "Disability" or any other single-letter codes
      INCOME %in% c("D", "Disability") ~ "Disability",
      
      # If it doesn't match any of the above, keep it as-is (or mark as "Other")
      TRUE ~ INCOME
    )
  )

#Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

tabyl(data_2018$INCOME_CLEAN)
tabyl(data_2023$INCOME_CLEAN)
tabyl(data$INCOME_CLEAN)  # all years

library(dplyr)
library(stringr)


#___________________________________________________________
#5. EDUCATION

#Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

tabyl(data_2018$EDUCATION)
tabyl(data_2023$EDUCATION)
tabyl(data$EDUCATION)  # all years

library(dplyr)

data <- data %>%
  mutate(
    EDUCATION_CLEAN = case_when(
      # 1) "<HS" stays as is
      EDUCATION == "<HS" ~ "<HS",
      
      # 2) HS / GED variants
      EDUCATION %in% c("HS", "HS Diploma", "GED") ~ "HS/GED",
      
      # 3) Some College variants
      EDUCATION %in% c("Some College", "Some college") ~ "Some College",
      
      # 4) Grad School variants
      EDUCATION %in% c("Grad School", "Completed Graduate School") ~ "Graduate School",
      
      # 5) "College" remains "College"
      EDUCATION == "College" ~ "College",
      
      # 6) Vocational
      EDUCATION == "Vocational" ~ "Vocational",
      
      # 7) "UnKN" can be recoded to "Unknown" (optional)
      EDUCATION == "UnKN" ~ "Unknown",
      
      # If none of the above apply, keep the original label
      TRUE ~ EDUCATION
    )
  )
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

tabyl(data_2018$EDUCATION_CLEAN)
tabyl(data_2023$EDUCATION_CLEAN)
tabyl(data$EDUCATION_CLEAN)  # all years

# Quick overview
glimpse(data)

#_________________________________________
#6. CALLER

#Subset data for 2018 and 2023
data_2018 <- data %>% filter(YEAR == "2018")
data_2023 <- data %>% filter(YEAR == "2023")

class(data$CALLER)
tabyl(data_2018$CALLER)
tabyl(data_2023$CALLER)
tabyl(data$CALLER)  # all years

library(dplyr)
library(stringr)
#______________________________________________

#Tables

library(janitor)
library(readr)  # or just use write.csv()

my_table1 <- tabyl(data_2018$AGE)
my_table2 <- tabyl(data_2023$AGE)
my_table3 <- tabyl(data$AGE) 

my_table4 <- tabyl(data_2018$GENDER)
my_table5 <- tabyl(data_2023$GENDER)
my_table6 <- tabyl(data$GENDER)

my_table7 <- tabyl(data_2018$'CONTACT LENGTH MINUTES')
my_table8 <- tabyl(data_2023$'CONTACT LENGTH MINUTES')
my_table9 <- tabyl(data$'CONTACT LENGTH MINUTES') 

my_table10 <- tabyl(data_2018$EMPLOYMENT)
my_table11 <- tabyl(data_2023$EMPLOYMENT)
my_table12 <- tabyl(data$EMPLOYMENT)

my_table13 <- tabyl(data_2018$TYPE_GAMBLING_CLEAN)
my_table14 <- tabyl(data_2023$TYPE_GAMBLING_CLEAN)
my_table15 <- tabyl(data$TYPE_GAMBLING_CLEAN)

my_table16 <- tabyl(data_2018$RELATIONSHIP_STATUS_CLEAN)
my_table17 <- tabyl(data_2023$RELATIONSHIP_STATUS_CLEAN)
my_table18 <- tabyl(data$RELATIONSHIP_STATUS_CLEAN)

my_table19 <- tabyl(data_2018$INCOME_CLEAN)
my_table20 <- tabyl(data_2023$INCOME_CLEAN)
my_table21 <- tabyl(data$INCOME_CLEAN)

my_table22 <- tabyl(data_2018$EDUCATION_CLEAN)
my_table23 <- tabyl(data_2023$EDUCATION_CLEAN)
my_table24 <- tabyl(data$EDUCATION_CLEAN)

my_table25 <- tabyl(data_2018$CALLER)
my_table26 <- tabyl(data_2023$CALLER)
my_table27 <- tabyl(data$CALLER)

write_csv(my_table1, "my_table1.csv")
write_csv(my_table2, "my_table2.csv")
write_csv(my_table3, "my_table3.csv")
write_csv(my_table4, "my_table4.csv")
write_csv(my_table5, "my_table5.csv")
write_csv(my_table6, "my_table6.csv")
write_csv(my_table7, "my_table7.csv")
write_csv(my_table8, "my_table8.csv")
write_csv(my_table9, "my_table9.csv")
write_csv(my_table10, "my_table10.csv")
write_csv(my_table11, "my_table11.csv")
write_csv(my_table12, "my_table12.csv")
write_csv(my_table13, "my_table13.csv")
write_csv(my_table14, "my_table14.csv")
write_csv(my_table15, "my_table15.csv")
write_csv(my_table16, "my_table16.csv")
write_csv(my_table17, "my_table17.csv")
write_csv(my_table18, "my_table18.csv")
write_csv(my_table19, "my_table19.csv")
write_csv(my_table20, "my_table20.csv")
write_csv(my_table21, "my_table21.csv")
write_csv(my_table22, "my_table22.csv")
write_csv(my_table23, "my_table23.csv")
write_csv(my_table24, "my_table24.csv")
write_csv(my_table25, "my_table25.csv")
write_csv(my_table26, "my_table26.csv")
write_csv(my_table27, "my_table27.csv")

###############################################################
# 0) LOAD PACKAGES 
###############################################################

library(janitor)  # for tabyl()
library(psych)    # for describe()

###############################################################
#CREATE A NEW DATASET WITH CALLER == "Self"
# keep only those rows where CALLER is exactly "Self."

data_self <- subset(data, CALLER == "Self")


nrow(data_self) #2021 rows

tabyl(data_self$CALLER) #shows only self

# CREATE SUBSETS FOR 2018 AND 2023 FROM data_self

data_2018_self <- subset(data_self, YEAR == "2018")
data_2023_self <- subset(data_self, YEAR == "2023")

# Quick checks
nrow(data_2018_self) #123 ROWS
nrow(data_2023_self) #347 ROWS


# 3) DESCRIPTIVE ANALYSES (AGE, CONTACT LENGTH, etc.)
#    *ALL* FOR CALLER == "Self" ONLY


# AGE
#---------------------------
tabyl(data_self$AGE)       # All years, Self only
tabyl(data_2018_self$AGE)  # 2018 only, Self
tabyl(data_2023_self$AGE)  # 2023 only, Self

my_table28 <- tabyl(data_self$AGE)
my_table29 <- tabyl(data_self$AGE)
my_table30 <- tabyl(data_2023_self$AGE) 
#---------------------------
# CONTACT LENGTH MINUTES
#---------------------------
# Check class
class(data_self$`CONTACT LENGTH MINUTES`)

# Descriptive stats
describe(data_self$`CONTACT LENGTH MINUTES`)       # All years, Self
describe(data_2018_self$`CONTACT LENGTH MINUTES`)  # 2018, Self
describe(data_2023_self$`CONTACT LENGTH MINUTES`)  # 2023, Self

my_table31 <- tabyl(data_self$`CONTACT LENGTH MINUTES`)
my_table32 <- tabyl(data_2018_self$`CONTACT LENGTH MINUTES`)
my_table33 <- tabyl(data_2023_self$`CONTACT LENGTH MINUTES`)

#---------------------------
# EMPLOYMENT
#---------------------------
tabyl(data_self$EMPLOYMENT)
tabyl(data_2018_self$EMPLOYMENT)
tabyl(data_2023_self$EMPLOYMENT)

my_table34 <- tabyl(data_self$EMPLOYMENT)
my_table35 <- tabyl(data_2018_self$EMPLOYMENT)
my_table36 <- tabyl(data_2023_self$EMPLOYMENT) 

#---------------------------
# TYPE OF GAMBLING (CLEAN)
#---------------------------
tabyl(data_self$TYPE_GAMBLING_CLEAN)
tabyl(data_2018_self$TYPE_GAMBLING_CLEAN)
tabyl(data_2023_self$TYPE_GAMBLING_CLEAN)

my_table37 <- tabyl(data_self$TYPE_GAMBLING_CLEAN)
my_table38 <- tabyl(data_2018_self$TYPE_GAMBLING_CLEAN)
my_table39 <- tabyl(data_2023_self$TYPE_GAMBLING_CLEAN)

#---------------------------
# GENDER
#---------------------------
tabyl(data_self$GENDER)
tabyl(data_2018_self$GENDER)
tabyl(data_2023_self$GENDER)

my_table40 <- tabyl(data_self$GENDER)
my_table41 <- tabyl(data_2018_self$GENDER)
my_table42 <- tabyl(data_2023_self$GENDER)

#---------------------------
# RELATIONSHIP STATUS (CLEAN)
#---------------------------
tabyl(data_self$RELATIONSHIP_STATUS_CLEAN)
tabyl(data_2018_self$RELATIONSHIP_STATUS_CLEAN)
tabyl(data_2023_self$RELATIONSHIP_STATUS_CLEAN)

my_table43 <- tabyl(data_self$RELATIONSHIP_STATUS_CLEAN)
my_table44 <- tabyl(data_2018_self$RELATIONSHIP_STATUS_CLEAN)
my_table45 <- tabyl(data_2023_self$RELATIONSHIP_STATUS_CLEAN)

#---------------------------
# INCOME (CLEAN)
#---------------------------
tabyl(data_self$INCOME_CLEAN)
tabyl(data_2018_self$INCOME_CLEAN)
tabyl(data_2023_self$INCOME_CLEAN)

my_table46 <- tabyl(data_self$INCOME_CLEAN)
my_table47 <- tabyl(data_2018_self$INCOME_CLEAN)
my_table48 <- tabyl(data_2023_self$INCOME_CLEAN)

library(dplyr)
library(stringr)

# Address 21-30k discrepancies:21-30k versus 21k-30k yr

data_self$INCOME_CLEAN <- ifelse(
  grepl("21k-30k yr|21-30k yr", tolower(data_self$INCOME_CLEAN)),
  "21-30k",
  data_self$INCOME_CLEAN
)

# 1) Identify rows that contain the string
i <- grepl("21k-30k yr|21-30k yr", tolower(data_self$INCOME_CLEAN))

# 2) For those rows, set INCOME_CLEAN to "21-30k"
data_self$INCOME_CLEAN[i] <- "21-30k"


#double check work and subsets
# CREATE SUBSETS FOR 2018 AND 2023 FROM data_self

data_2018_self <- subset(data_self, YEAR == "2018")
data_2023_self <- subset(data_self, YEAR == "2023")

library(janitor)  # for tabyl()
library(psych)    # for describe()

# Quick checks
nrow(data_2018_self) #123 ROWS
nrow(data_2023_self) #347 ROWS

tabyl(data_2018_self$INCOME_CLEAN)
tabyl(data_2023_self$INCOME_CLEAN)
tabyl(data_self$INCOME_CLEAN)  # all years

my_table46 <- tabyl(data_self$INCOME_CLEAN)
my_table47 <- tabyl(data_2018_self$INCOME_CLEAN)
my_table48 <- tabyl(data_2023_self$INCOME_CLEAN)

#---------------------------
# EDUCATION (CLEAN)
#---------------------------
tabyl(data_self$EDUCATION_CLEAN)
tabyl(data_2018_self$EDUCATION_CLEAN)
tabyl(data_2023_self$EDUCATION_CLEAN)

my_table49 <- tabyl(data_self$EDUCATION_CLEAN)
my_table50 <- tabyl(data_2018_self$EDUCATION_CLEAN)
my_table51 <- tabyl(data_2023_self$EDUCATION_CLEAN)

###############################################################
# Viewing and Inspecting a Data Frame in Base R
###############################################################

#See just the column names (variables)
names(data_self)

# 5) Summary of each column 
summary(data_self)

# 6) open a spreadsheet-like viewer in a new tab:
View(data_self)


