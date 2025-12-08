# **Text-2025 branch:**

## **This branch includes the discussion avenues of *Uncinaria stenocephala* over the 2023-2025 period and comparisons to the historic-2022 period. All analyses were conducted in R.**

### **Inclusion criteria:**

The previous data from *U. stenocephala* text mining of the historic to 2022 period (branch "Text-2022") is mined with the same criteria, only with an additional key word ("resistance"), as such it is not featured in this branch of the repository. The novel data (from 2023-2025) for *U. stenocephala* was mined using a similar pipeline as the Text-2022 pipeline, regardless the code is featured below. The PubMed search was conducted on the 8th of December 2025.

### **PubMed data mining code:**

The following libraries are used in this branch:

```         
library(rentrez)
library(XML)
library(tm)
library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
```

The data mining code is provided below: for greater detail see the Text-2022 branch:

```         
2025_search <- entrez_search(db="pubmed", term="Uncinaria stenocephala[Title/Abstract]", mindate = 2023, maxdate = 2025, retmax = 1000)

2025_search

abstracts <- entrez_fetch(db="pubmed", id=2025_search$ids, rettype="abstract", retmode="text")

head(abstracts)

docs <- Corpus(VectorSource(abstracts))

inspect(docs)

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
docs <- tm_map(docs, removeWords, stopwords("en")) 
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
# Text stemming
#docs <- tm_map(docs, stemDocument)

dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 100)

target_words <- c("stenocephala", "gene", "protein", "proteins",
                  "molecular", "prevalence", "resistance")

freq_Uste <- data.frame(
  word = target_words,
  freq = sapply(target_words, function(w) {
    if(w %in% rownames(m)) rowSums(m)[w] else 0
  })
)

write.csv(freq_Uste, "Uste_word_counts-2025.csv", row.names = FALSE)
```

The csv files were merged in excel and the protein/proteins rows merged.

**And now for the graphing:**

```         
merged_csv<- read.csv("text-mining-2022-2025.csv")

head(merged_csv)

# Divisors for the numeric columns (skip the first categorical column)
divisors <- c(192, 28)

# Apply sweep only to numeric columns
df_divided <- merged_csv
df_divided[, 2:3] <- sweep(merged_csv[, 2:3], 2, divisors, FUN="/")

df_divided

pivoted<- df_divided %>% pivot_longer(cols=c(year_2022, year_2025), names_to="Hookworms", values_to="Word_freq")

pivoted

custom_order <- c("caninum", "stenocephala", "prevalence", "resistance", "molecular", "protein", "gene")

pivoted$word <- factor(pivoted$word, levels = custom_order)

plot<- ggplot(pivoted, aes(x=word, y=Word_freq, fill=Hookworms))+
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Word", y = "Frequency", title = "Average word Frequencies in Abstracts") +
  theme_minimal() + theme_prism(base_size = 12) +        # Prism-style theme
  scale_fill_prism(palette = "colors") +  # Prism color palette
  theme(axis.text.x = element_text(angle =90, hjust = 1))

plot
```

**The accompanying csv file can be found in this branch under "text-mining-2022-2025.csv"**
