#remove everything 
rm(list=ls()); gc()

library(tidyverse)

#read the data
data <- read.csv("Week-1/data/gapminder_clean.csv")

saveRDS(data,"Week-1/data/data.RDS") #save the data as RDS file
readr::write_csv(data,"Week-1/data/data.csv") #save the data as csv file

#ignore from pushing the data
usethis::use_git_ignore("Week-1/data/data.csv")
usethis::use_git_ignore("Week-1/data/data.RDS")

#read different formats of data
data <- read.csv("Week-1/data/data.csv")
data <- readRDS("Week-1/data/data.RDS")

data <- read.csv("Week-1/data/gapminder_clean.csv")
head(data)
dim(data)

#show the countries
table(data$Country.Name)

#filter for the year 1962
data1962 <- data %>%
  dplyr::filter(Year == 1962) %>%
  dplyr::rename("CO2_emissions" = "CO2.emissions..metric.tons.per.capita.")

data1962 <- data %>%
  filter(Year == 1962) %>%
  rename("CO2_emissions" = "CO2.emissions..metric.tons.per.capita.")


#plot
p <- ggplot(data1962, aes(
  x = CO2_emissions,
  y = gdpPercap
)) +
  geom_point() +
  labs(x = "log (CO2 emissions (metric tons/capita))",
       y = "log (GDP per cap)",
       title = "CO2 emissions and GDP per capita (1962)") +
  #scale_x_log10() +
  #scale_y_log10() +
  geom_smooth(method = "lm", se = FALSE)

p

anyNA(data) #check missing values


cont_energy <- data %>%
  rename("energy_use" = "Energy.use..kg.of.oil.equivalent.per.capita.")%>%
  drop_na(continent, energy_use, Year) %>%
  select(continent, energy_use, Year)

cont_energy %>%
ggplot(aes(x=continent, y=energy_use)) +
           geom_boxplot() +
           scale_y_log10() +
           ylab("Energy use (kg of oil equivalent/capita)") +
           ggtitle("Global energy use (1962-2007)") 

#Country with highest population density over the years
gapminderd <- data %>%
  rename("popdensity"= "pop","country" = "Country.Name")%>%
  group_by(country)%>%
  filter(Year >= 1962, Year <= 2007)%>%
  summarise(medpopden = median(popdensity)) %>%
  arrange(desc(medpopden))%>%
  top_n(5, medpopden) %>%
  ungroup()

gapminderd


ggplot(roundedGMD, aes(x=country, y = medpopden, fill = country)) +
           labs(title = "Countries with highest population density (people/sq.km) between 1962-2007",
                y = "median population density (1962-2007)") +
           geom_col() + geom_text(aes(label=medpopden)) +
           theme (axis.text.x = element_text(angle=30, size=9,vjust=.8, hjust=0.8))


#################################################################################
########################### Model Training #####################################

library(gapminder)

ggplot(gapminder_unfiltered, aes(gdpPercap, lifeExp)) +
  geom_point()

ggplot(gapminder_unfiltered, aes(gdpPercap, lifeExp)) +
  geom_point() + 
  geom_smooth(method = "lm")

lm(lifeExp ~ gdpPercap, data = gapminder_unfiltered)

#log transform the data
ggplot(gapminder_unfiltered, aes(log(gdpPercap), lifeExp)) +
  geom_point() + 
  geom_smooth(method = "lm")

#take the summary model
summary(lm(lifeExp ~ log(gdpPercap), data = gapminder_unfiltered))


#let's make some predictions
gapMod <- lm(lifeExp ~ log(gdpPercap) + continent + year + (log(gdpPercap)):continent, data = gapminder_unfiltered)

gapPred <- gapminder_unfiltered
gapPred <- gapPred %>%
  mutate(predict = predict(gapMod, newdata = gapPred))

head(gapPred)

ggplot(gapPred, aes(gdpPercap)) + 
  geom_point(aes(y = lifeExp)) + 
  geom_point(aes(y = predict), color = "blue", alpha = 0.25)

######################## CLASSIFICATION ######################################

head(mtcars)

#summary plot
ggplot(mtcars, aes(mpg, am)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F)

#train a logistic regression model
logmod <- glm(am ~ mpg, data = mtcars, family = binomial)
summary(logmod)
