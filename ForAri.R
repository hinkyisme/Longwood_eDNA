#Jameson Hinkle for Ari
#12/18/24
#Testing for independence and/or relationship of 
#Calls to Gambling

# Dependencies
require(tidyverse)

test <- read_csv("VPGH.Callers.Total.Original.csv")

# Check for NAs in Date, data string.
contains_na <- any(is.na(test$DATE))
if (contains_na) {
  print("This list contains NAs")
 } else{
    print("The list does not contain NAs")
 }

# Convert character class to year class only.
test <- test %>% mutate(year = DATE %>%
as.Date(format = "%m/%d/%y") %>% format("%Y"))

# Make year a factor for downstream analysis

test$year <- as.factor(test$year)

# Convert age ranges to its own factor/variable

test$AGE <- as.factor(test$AGE)

# Create contigency table and perform Chi-square test for independence RE RQ1

contingency_table <- test %>%
  count(year, LOCATION, `TYPE OF GAMBLING`, GENDER, AGE, CALLER, `MARITAL STATUS`, INCOME, EMPLOYMENT, EDUCATION, CONTACT, `48HR FOLLOW-UP`, `48HR GAMBLING STATUS`, `ONE MONTH`, `1 MONTH GAMBLING STATUS`, REFERRAL, `HOW HEARD`, `TYPE OF CONTACT`) %>%   # Count occurrences for each combination
  pivot_wider(names_from = c(LOCATION, `TYPE OF GAMBLING`, GENDER, AGE, CALLER, `MARITAL STATUS`, INCOME, EMPLOYMENT, EDUCATION, CONTACT, `48HR FOLLOW-UP`, `48HR GAMBLING STATUS`, `ONE MONTH`, `1 MONTH GAMBLING STATUS`, REFERRAL, `HOW HEARD`, `TYPE OF CONTACT`), values_from = n, values_fill = 0) %>%
  column_to_rownames("year") %>%
  as.matrix()

chi_square_result <- chisq.test(contingency_table)
