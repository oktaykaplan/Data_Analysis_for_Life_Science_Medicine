# BENG 628 - Lecture 2: Exploratory Data Analysis for Life Sciences

# Lecture Goals

#Deepen understanding of data manipulation with dplyr and tidyverse
#Explore more complex visualizations with ggplot2
#Perform multivariate statistical analysis & survival analysis
#Students need to learn how to clean data with these powerful tools using concise and expressive code
#Practical Data Wrangling Tools
#dplyr for Data Manipulation:
#Filtering: filter() to select specific rows.
#Selecting: select() to choose specific columns.
#Mutating: mutate() to create new variables.
#Grouping: group_by() and summarise() to calculate group-level statistics.

#You can download a dataset from Kaggle or WHO's Open Data Portal:
#Kaggle: Global Health Statistics Dataset
#https://www.kaggle.com/datasets/malaiarasugraj/global-health-statistics
#Global Health Statistics.csv(134.4 MB)

# Practical Data Wrangling Tools

## Load Required Libraries

# Part 1: R Fundamentals and Advanced Data Wrangling

# Load essential packages
library(tidyverse)  # Includes ggplot2, dplyr, tidyr
library(ggpubr)     # For publication-ready plots
library(rstatix)    # For statistical analysis
library(pheatmap)   # For heatmaps
library(DESeq2)     # For RNA-seq analysis
library(biomaRt)    # For accessing biological databases
library(viridis)    # For colorblind-friendly palettes
library(ComplexHeatmap) # For advanced heatmaps
library(factoextra)  # For PCA visualization
library(readr)


#Once you've downloaded the dataset, let's clean it step by step.
# Load dataset (update path as needed)
health_data <- read_csv("/Users/sebihacevik/Downloads/Global Health Statistics.csv")
# View the first few rows
head(health_data)
colnames(health_data)
print(health_data)
View(health_data)

#1. Check for Missing Values
# Count missing values in each column
colSums(is.na(health_data))


# Visualize missing data
library(naniar)
gg_miss_var(health_data)


#2. Remove Duplicates
# Check for duplicate rows
sum(duplicated(health_data))

# Remove duplicate rows
health_data <- health_data %>% distinct()
health_data 


#3. Fix Inconsistent or Incorrect Data
#(a) Correct Disease Categories
#Some diseases are incorrectly classified

unique(health_data$Country)
unique(health_data$`Disease Name`)
# Check unique categories
unique(health_data$`Disease Category`)

# Correct specific values
#This part uses the mutate() function (from the dplyr package, part of the tidyverse)
#health_data <- health_data %>% : This assigns the modified data frame back to the health_data variable, effectively updating the data
#The %>% is the pipe operator (also from dplyr).
#mutate(Disease Category= ...): This creates or modifies the Disease Category column. 
#The recode() function is used to change the values within this column.
#Wherever you find the value "Genetic", replace it with "Infectious".
#Wherever you find the value "Bacterial", replace it with "Infectious"
health_data <- health_data %>%
  mutate(`Disease Category` = recode(`Disease Category`,
                                     "Genetic" = "Infectious",  # Example of misclassification
                                     "Bacterial" = "Infectious"
  ))
#(b) Standardize "Availability of Vaccines/Treatment
# Check inconsistent values
table(health_data$`Availability of Vaccines/Treatment`)
table(health_data$`Disease Name`)

# Convert "Yes", "Available", "TRUE" → "Yes"
#ifelse(...): This is a conditional statement. It checks a condition and returns one value if the condition is true, and another value if the condition is false.
#%in% is an operator that tests if a value is present within a vector.
health_data$`Availability of Vaccines/Treatment` <- ifelse(
  health_data$`Availability of Vaccines/Treatment` %in% c("Yes", "Available", "TRUE"), 
  "Yes", 
  "No"
)

head(health_data)
health_data$`Availability of Vaccines/Treatment`
health_data$`Availability of Vaccines/Treatment`[1:10]  # Shows the first 10 values
health_data$`Availability of Vaccines/Treatment`[c(1, 5, 10, 20)] # Shows values at specific row numbers

#4. Fix Column Data Types
#Sometimes, data is imported into R with incorrect data types. This needs to be corrected
#✅ Fixes issues where some numerical columns might be treated as text.
#`Incidence Rate (%)` = as.numeric(`Incidence Rate (%)`): 
#`Mortality Rate (%)` = as.numeric(`Mortality Rate (%)`)  Converts the Mortality Rate (%) column to numeric
colnames(health_data)
str(health_data)

health_data <- health_data %>%
  mutate(
    Year = as.integer(Year),
    `Prevalence Rate (%)` = as.numeric(`Prevalence Rate (%)`),
    `Incidence Rate (%)` = as.numeric(`Incidence Rate (%)`),
    `Mortality Rate (%)` = as.numeric(`Mortality Rate (%)`),
    `Population Affected` = as.numeric(`Population Affected`),
    `Per Capita Income (USD)` = as.numeric(`Per Capita Income (USD)`),
    `Urbanization Rate (%)` = as.numeric(`Urbanization Rate (%)`)
  )

# Select Country, Year, and Disease Name
health_data %>% dplyr::select(Country, Year, `Disease Name`)

# Select all columns EXCEPT Population Affected and DALYs
health_data %>% dplyr::select(-`Population Affected`, -DALYs)

# Select columns starting with "P" (using a regular expression)
health_data %>% dplyr::select(starts_with("P"))

#arrange(): Sorting rows

# Sort by Prevalence Rate (%) in descending order
health_data %>% dplyr::arrange(desc(`Prevalence Rate (%)`))

# Sort by Country, then by Year (ascending)
health_data %>% dplyr::arrange(Country, Year)
# Sort by Country, then by Year (ascending)
health_data %>% dplyr::arrange(`Treatment Type`)
colnames(health_data)

health_data %>% dplyr::arrange(Country, Year, `Disease Name`)

# Sort by Disease Category, then by Incidence Rate (%) (descending)
health_data %>% dplyr::arrange(`Disease Category`, desc(`Incidence Rate (%)`))

#group_by() and summarize(): Grouped calculations

# Calculate the mean prevalence rate for each disease category
health_data %>%
  dplyr::group_by(`Disease Category`) %>%
  summarize(mean_prevalence = mean(`Prevalence Rate (%)`, na.rm = TRUE),
            median_prevalence = median(`Prevalence Rate (%)`, na.rm = TRUE),
            n = n()) # n() counts the number of observations in each group

# Calculate the total population affected for each country
health_data %>%
  group_by(Country) %>%
  summarize(total_affected = sum(`Population Affected`, na.rm = TRUE))

# Calculate the average treatment cost and recovery rate for each disease category and year.
health_data %>%
  dplyr::group_by(`Disease Category`, Country) %>%
  summarize(avg_cost = mean(`Average Treatment Cost (USD)`, na.rm = TRUE),
            avg_recovery = mean(`Recovery Rate (%)`, na.rm = TRUE))


#Joining data frames
# Create a second data frame (example) - replace with your actual data

country_data <- tibble(
  Country = c("Italy", "France", "Turkey", "Indonesia", "Saudi Arabia"),
  Continent = c("Europe", "Europe", "Asia", "Asia", "Asia"),
  Population = c(60000000, 65000000, 85000000, 270000000, 35000000)
)

View(country_data)

# Inner join: Keep only rows where the country is present in BOTH data frames
joined_data <- health_data %>%
  inner_join(country_data, by = "Country")

View(joined_data)
# Left join: Keep all rows from health_data, and matching rows from country_data.
# If there's no match in country_data, you'll have NA values.
joined_data_left <- health_data %>%
  left_join(country_data, by = "Country")

# Right join: Keep all rows from country_data, and matching rows from health_data.
joined_data_right <- health_data %>%
  right_join(country_data, by = "Country")

# Full join: Keep all rows from BOTH data frames.
joined_data_full <- health_data %>%
  full_join(country_data, by = "Country")


print(joined_data) # Print the joined data frame
print(joined_data_left) # Print the joined data frame
print(joined_data_right) # Print the joined data frame
print(joined_data_full) # Print the joined data frame



# Handle Outliers
library(ggplot2)
colnames(health_data)
head(health_data)
ggplot(health_data, aes(x = `Per.Capita.Income..USD.`)) +
  geom_boxplot()

#✅ Identify extreme outliers (e.g., Per Capita Income (USD) values that don't make sense).
#The filter() function selects rows from a data frame based on a condition or set of conditions.  Only rows that meet the condition(s) are kept; the rest are discarded.
health_data <- health_data %>%
  filter(`Per.Capita.Income..USD.` > 100 & `Per.Capita.Income..USD.` < 100000)



#6. Feature Engineering
#(a) Create a Health Risk Index
health_data <- health_data %>%
  mutate(Health_Risk_Index = (`Prevalence Rate (%)` + `Incidence Rate (%)` + `Mortality Rate (%)`) / 3)
colnames(health_data)

#✅ This new column gives an overall risk score for each disease.

#(b) Categorize Countries Based on Healthcare Access
health_data <- health_data %>%
  mutate(Healthcare_Access_Level = case_when(
    `Healthcare Access (%)` > 80 ~ "High",
    `Healthcare Access (%)` > 50 ~ "Medium",
    TRUE ~ "Low"
  ))

#7. Save Cleaned Data

write_csv(health_data, "/Users/sebihacevik/Downloads/WHO_Health_Cleaned.csv")

#1️⃣ Exploratory Data Analysis (EDA)
#Before jumping into modeling, understand the data through visualization and summary statistics.

#Summary statistics:
#This command summarizes each column in your dataset.
#It provides minimum, maximum, mean, median (50% quantile), and quartiles (25%, 75%) for numerical columns.
#For categorical variables, it shows the frequency of each category.
summary(health_data)


#Check distributions
#To visualize how data is spread, we use histograms and boxplots
#Creates a histogram (bar chart) of Prevalence Rate (%).
hist(health_data$Year)
hist(health_data$`Incidence Rate (%)`)
hist(health_data$`Prevalence Rate (%)`)
boxplot(health_data$`Incidence Rate (%)`)

#Correlation analysis
#Calculates the correlation matrix between all numeric columns.
#Measures how strongly two variables are related.
#Outputs values between -1 and 1:
#1 → Perfect positive correlation (e.g., height & weight).
#0 → No correlation (random relationship).
#-1 → Perfect negative correlation (e.g., smoking & lung capacity).

cor(health_data[, sapply(health_data, is.numeric)])


#2️⃣ Data Visualization
#Let's see disease distribution by country
library(ggplot2)
ggplot(health_data, aes(x = Country)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("This is cool")

#Prevalence vs. Incidence Rate scatter plot
ggplot(health_data, aes(x = `Prevalence Rate (%)`, y = `Incidence Rate (%)`, color = `Disease Category`)) +
  geom_point() +
  theme_minimal()


#3️⃣ Feature Engineering & Transformation
#Create a new column for Disease Severity Score
#You learn how to create a column
health_data$Disease_Severity <- health_data$`Mortality Rate (%)` * health_data$DALYs

#Convert categorical variables into factors
health_data$`Disease Category` <- as.factor(health_data$`Disease Category`)

#4️⃣ Predictive Modeling (if you want)
model <- lm(`Mortality Rate (%)` ~ `Prevalence Rate (%)` + `Incidence Rate (%)` + `Healthcare Access (%)`, data = health_data)
summary(model)
#Random Forest (for classification or regression)
library(randomForest)
#Since "Disease Name" contains spaces, it might not be recognized properly.
colnames(health_data) <- make.names(colnames(health_data))  # Converts spaces to dots
#If "Disease Name" is missing, check if it was accidentally removed.
colnames(health_data)
rf_model <- randomForest(`Mortality.Rate....` ~ ., data = health_data, importance = TRUE)
#Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()


#🔧 Solutions to Reduce Memory Usage
#1️⃣ Reduce the Number of Features (Columns)
#Some columns (like "Country" or "Disease Name") may not be useful for predicting mortality rate. Try keeping only numeric features:

health_data_clean <- health_data[, sapply(health_data, is.numeric)]
rf_model <- randomForest(`Mortality.Rate....` ~ ., data = health_data_clean, importance = TRUE)
#Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

#✅ 1️⃣ Use the ranger Package (Best for Large Datasets)
library(ranger)

rf_model <- ranger(`Mortality.Rate....` ~ ., 
                   data = health_data_clean, 
                   num.trees = 100, 
                   importance = "impurity", 
                   write.forest = FALSE)  # Avoids storing large trees in memory

rf_model

###################Data Wrangling and Cleaning 2###################

library(dplyr)
cleaned_data <- health_data %>%
  dplyr::filter(Prevalence.Rate.... > 1) %>%
  dplyr::select(Country, Year, `Prevalence.Rate....`)

cleaned_data

#2️⃣  Data Reshaping and Pivoting
#Use pivot_longer() and pivot_wider() from the tidyr package to reshape data efficiently.
library(tidyr)
health_data_long <- pivot_longer(health_data, 
                                 cols = c(`Prevalence.Rate....`, `Incidence.Rate....`, `Mortality.Rate....`), 
                                 names_to = "Rate_Type", 
                                 values_to = "Rate")

#Outlier Detection and Handling
boxplot(health_data$`Prevalence.Rate....`)

# Handling Missing Values
colnames(health_data)
health_data %>% filter(is.na(`Prevalence.Rate....`))

# Outlier Handling (Example)
# For quantiles
Q1 <- quantile(health_data$`Per.Capita.Income..USD.`, 0.25)
Q3 <- quantile(health_data$`Per.Capita.Income..USD.`, 0.75)

IQR <- Q3 - Q1
upper_bound <- Q3 + 1.5 * IQR
health_data <- health_data %>% filter(`Per Capita Income (USD)` < upper_bound)

# For binning
health_data <- health_data %>%
  mutate(Income_Level = cut(`Per.Capita.Income..USD.`, breaks = c(0, 10000, 50000, Inf), labels = c("Low", "Medium", "High")))

# One-Hot Encoding (Example)
health_data$Gender <- factor(health_data$Gender)  # Convert to factor
dummy_vars <- model.matrix(~ Gender - 1, data = health_data) # Create dummy variables
health_data <- cbind(health_data, dummy_vars) # Add to data frame

# Correlation Heatmap
library(ComplexHeatmap)
cor_matrix <- cor(health_data[, sapply(health_data, is.numeric)], use = "pairwise.complete.obs") # Handle NAs in correlation
Heatmap(cor_matrix,  cluster_rows = FALSE, cluster_columns = FALSE)

# t-test (Example)
library(rstatix)
health_data %>% t_test(`Prevalence.Rate....` ~ Gender)


#Ratios: Population size to number of doctors.
#Binning: Group continuous variables into categories (e.g., low, medium, high income).
#Interaction terms: Combine features for richer models.

health_data_clean <- health_data %>%
  mutate(Doctor_to_Population_Ratio = Doctors.per.1000 / Population.Affected)
health_data_clean

#Handling Categorical Data for Modeling
#Encoding Categorical Variables:
#Learn how to factorize categorical data for regression and machine learning algorithms.
#Discuss One-Hot Encoding vs. Label Encoding
health_data$Gender <- factor(health_data$Gender)

#Data Exploration and Visualization
#Visualize data distribution with histograms or density plots (using ggplot2)
library(ggplot2)
ggplot(health_data, aes(x = `Prevalence.Rate....`, y = `Mortality.Rate....`)) +
  geom_point() +
  theme_minimal()

