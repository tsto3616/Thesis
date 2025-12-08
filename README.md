# **Uste-prevalence branch**

## **This branch details the inclusion criteria and the coding required to generate the map of Australia with overlapping prevalence rates of *Uncinaria stenocephala.***

### **Inclusion of *U. stenocephala* prevalence studies:**

A search of the prevalence rates of *U. stenocephala* was conducted in PubMed on the 17th of November 2025 using the search terms **"Uncinaria stenocephala Australia", "Uncinaria stenocephala New Zealand", "Hookworms Australia" and "Hookworms New Zealand"**. To be included in the map the Authors had to detail the location of the sampling with an identifiable GPS coordinate and had to detect at least one species specific identification of *U. stenocephala*. In the absence of either inclusion criteria the study would be excluded. By species specific identification the criteria was either molecular methods or morphology upon post mortems of the hosts.

### **Coding pipeline in R:** 

Note the supplementary data collected from the search of PubMed is available under this branch for the name "Uste-prevalence-spreadsheet.csv".In brief the code below details the pipeline for the generation of the map:

Load the necessary libraries:

```         
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
```

Load the country at a medium level and import the data:

```         
# Load Australia
australia <- ne_countries(scale = "medium", country = "Australia", returnclass = "sf")

# Uste prevalence data
data <- read.csv("Uste-prevalence-spreadsheet.csv")
```

Now for the graphing of the prevalence rates onto the map of Australia:

```         
Uste<- ggplot(data = australia) +
  geom_sf(fill = "white", color = "black") +
  geom_point(data = data, aes(x = Longitude, y = Latitude, color = Prevalence, size = Samples), alpha = 0.5) + scale_color_gradientn(colors = c("yellow", "red"), limits= c(0, 100)) + ggtitle("Australia with Uncinaria stenocephala prevalence rates") +
  scale_size_continuous(limits = c(0, 800), range = c(2, 6)) + coord_sf(ylim = c(-5, -45), xlim = c(110, 155)) + theme_minimal()

Uste
```

And **done**.
