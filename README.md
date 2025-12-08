# Text-2022 branch

## **This branch represents the inclusion and coding for the data mining of text for *Uncinaria stenocephala* and its search results in PubMed from 2022 and earlier. All analyses in this branch used R.**

### **Inclusion criteria:**

To be included in the data mining the results from the search of abstracts had to include the search terms "Uncinaria stenocephala" and "Ancylostoma caninum" (two separate searches).

### **PubMed data mining code**:

The data was mined and analysed using the following libraries:

```         
library(rentrez)
library(XML)
library(tm)
```

Then the data was mined via rentrez for both *U. stenocephala* and *A. caninum.* The number of retrieved articles differed greatly for the two species, with *U. stenocephala* retrieving 192 articles from the past to the end of 2022 and *A. caninum* retrieving 799 articles in the same timeline:

```         
Uste_search <- entrez_search(db="pubmed", term="Uncinaria stenocephala[Title/Abstract]", mindate = 0001, maxdate = 2022, retmax = 1000)

Uste_search

abstracts <- entrez_fetch(db="pubmed", id=Uste_search$ids, rettype="abstract", retmode="text")

head(abstracts)

Acan_search <- entrez_search(db="pubmed", term="Ancylostoma caninum[Title/Abstract]", mindate = 0001, maxdate = 2022, retmax = 1000, use_history=TRUE)

Acan_search #799 results

Acan_abstracts <- entrez_fetch(db="pubmed", web_history=Acan_search$web_history, rettype="abstract", retmode="text")

head(Acan_abstracts)
```

Next they were converted into a format that could be computationally manipulated in R:

```         
# first Uste
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

# now Acan
docs_Acan <- Corpus(VectorSource(Acan_abstracts))

inspect(docs_Acan)

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs_Acan <- tm_map(docs_Acan, toSpace, "/")
docs_Acan <- tm_map(docs_Acan, toSpace, "@")
docs_Acan <- tm_map(docs_Acan, toSpace, "\\|")

# Convert the text to lower case
docs_Acan <- tm_map(docs_Acan, content_transformer(tolower))
# Remove numbers
docs_Acan <- tm_map(docs_Acan, removeNumbers)
# Remove english common stopwords
docs_Acan <- tm_map(docs_Acan, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
docs_Acan <- tm_map(docs_Acan, removeWords, stopwords("en")) 
# Remove punctuations
docs_Acan <- tm_map(docs_Acan, removePunctuation)
# Eliminate extra white spaces
docs_Acan <- tm_map(docs_Acan, stripWhitespace)
# Text stemming
#docs <- tm_map(docs, stemDocument)

dtm_Acan <- TermDocumentMatrix(docs_Acan)
m_Acan <- as.matrix(dtm_Acan)
v_Acan <- sort(rowSums(m_Acan),decreasing=TRUE)
d_Acan <- data.frame(word = names(v_Acan),freq=v_Acan)
head(d_Acan, 100)
```

Now that the files are manipulated we can proceed to export key words as a csv file for graphing and analyses in R:

```         
target_words <- c("stenocephala", "caninum", "gene", "protein", "proteins",
                  "molecular", "prevalence", "australia", "dogs", "dog")

# For Uncinaria stenocephala abstracts
freq_Uste <- data.frame(
  word = target_words,
  freq = sapply(target_words, function(w) {
    if(w %in% rownames(m)) rowSums(m)[w] else 0
  })
)

# For Ancylostoma caninum abstracts
freq_Acan <- data.frame(
  word = target_words,
  freq = sapply(target_words, function(w) {
    if(w %in% rownames(m_Acan)) rowSums(m_Acan)[w] else 0
  })
)

write.csv(freq_Uste, "C:/Users/tsto3616/thesis-drafts/Uste_word_counts.csv", row.names = FALSE)
write.csv(freq_Acan, "C:/Users/tsto3616/thesis-drafts/Acan_word_counts.csv", row.names = FALSE)
```

**The resultant csv files were merged and had plural and singular key words merged - it can be found under this branch as "2022-text-mining.csv"**
