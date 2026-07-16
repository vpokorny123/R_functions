df<-data.table::fread('/Users/victorpokorny/Downloads/SearchResult_Export_07May2025_124958.csv')
#let's look at just k99 for now
k99_df<-df[df$Activity=='K99']

#remove repeats (like year 1 vs year 2 etc.)
k99_df<-k99_df[!duplicated(k99_df$`Contact PI Person ID`),]

str_spk99_df$`Project Abstract`

# Load necessary libraries
library(dplyr)
library(stringr)
library(tidytext)
library(tibble)

# Example: Assume your data frame is called df and the column of abstracts is named "abstract"
# Step 1: Convert the abstracts into a tidy format
word_counts <- df %>%
  # Make sure text is character
  mutate(abstract = as.character(`Project Abstract`)) %>%
  # Tokenize the abstracts into words
  unnest_tokens(word, abstract) %>%
  # Remove common stop words (optional)
  anti_join(stop_words, by = "word") %>%
  # Count the number of abstracts each word appears in
  distinct(word, row_number()) %>%
  count(word, sort = TRUE)

# Step 2: Find words common to *all* abstracts
# Count how many unique abstracts there are
num_abstracts <- nrow(df)

# Step 3: Filter to keep only words that appear in all abstracts
common_words <- word_counts %>%
  filter(n == num_abstracts)

# Output
print(common_words)