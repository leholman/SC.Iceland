
setwd("~/GitHubRepos/SC.Iceland/taxonomy/functionalAnno")

data <- readLines("alnout.txt")
df <- data.frame(data, stringsAsFactors = FALSE)

# Extract ASV and category

df$ASV <- regmatches(df$data, regexpr("ASV_\\d+", df$data))
df$Category <- sub("^ASV_\\d+\\s+([^;]+);.*", "\\1", df$data)

# Resolve conflicts
resolved <- aggregate(Category ~ ASV, df, function(x) if(length(unique(x)) > 1) "Unknown" else unique(x))

PR2assignments <- read.csv("EUK.cleaned.PR2.csv")

PR2assignments$cat <- resolved$Category[match(PR2assignments$X.1,resolved$ASV)] 

PR2assignments$cat[is.na(PR2assignments$cat)] <- "Unknown"

write.csv(PR2assignments,"EUK.cleaned.PR2.cat.csv")

