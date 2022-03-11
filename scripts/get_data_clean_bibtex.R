## we first create the order using the information contained in the bib file
## we then extract the citation data

rm(list = ls())
library(bib2df)
library(stringr)
library(dplyr)
filestr <- "data/SCC2016/4Journals_cleaned.bib"
df <- bib2df(filestr)
head(df)
dim(df)
colnames(df)

df <- as.data.frame(df)
# df[1,]

# df$DOI
# df$JOURNAL
# df$NUMBER
# df$PAGES
# df$VOLUME ## volume is equivalent to YEAR
# table(df$VOLUME[df$JOURNAL == journals[5]],df$YEAR[df$JOURNAL == journals[5]])

# sapply(df, function(x) any(!is.na(x)))
# sapply(df, function(x) any(is.na(x)))

path = "data/SCC2016/"
paperList = read.table(paste0(path, "paperList.txt"), header = T, sep = ",")
paperDOI = as.character(paperList$DOI)

paperCitAdj = read.table(paste0(path, "paperCitAdj.txt"), header = F, sep = " ")
paperCitAdj = as.matrix(paperCitAdj)
colnames(paperCitAdj) = NULL
rownames(paperCitAdj) = NULL

xx  = rowSums(paperCitAdj)      # number of citation received by each paper
yy  = colSums(paperCitAdj)      # number of citation given by each paper

# same order as paperList
all.equal(paperDOI, df$DOI)
# same order as paperCitAdj
all.equal(paperList$citCounts, xx)

journals <- unique(df$JOURNAL)
df$JOURNAL[df$JOURNAL == journals[4]] <- journals[5]
journals <- unique(df$JOURNAL)

unique(df$VOLUME[df$JOURNAL == journals[1]])
unique(df$NUMBER[df$JOURNAL == journals[1]])

table(df$VOLUME[df$JOURNAL == journals[1]],df$NUMBER[df$JOURNAL == journals[1]])

unique(df$PAGES)
df$PAGES[149] <- "2823--2856" ## originally "2823-2856"
tmp <- str_split(df$PAGES, "--")
df$PAGES1 <- lapply(tmp, function(x)x[1])
df$PAGES2 <- lapply(tmp, function(x)x[2])
df$PAGES1 <- as.numeric(df$PAGES1)
df$PAGES2 <- as.numeric(df$PAGES2)

table(df$YEAR[df$JOURNAL == journals[1]],df$NUMBER[df$JOURNAL == journals[1]])
df$NUMBER[1] <- "2" ## originally NA and X. == 2
df$NUMBER_mod <- df$NUMBER
df$NUMBER_mod[df$NUMBER_mod == "6A"] <- "6"
df$NUMBER_mod[df$NUMBER_mod == "6B"] <- "6"
df$NUMBER_mod[df$NUMBER_mod == "5A"] <- "5"
df$NUMBER_mod[df$NUMBER_mod == "5B"] <- "5"
subindex <- as.numeric(df$NUMBER_mod) > 400
df$NUMBER_mod <- as.numeric(df$NUMBER_mod)
df$NUMBER_mod[subindex] <- (df$NUMBER_mod[subindex] - 1) %% 4 + 1

## we don't need number because within a year the pages are increasing (and unique)
# df$PAGES1[df$JOURNAL == journals[1]]
# plot(df$YEAR[df$JOURNAL == journals[1]], df$PAGES1[df$JOURNAL == journals[1]])

# table(df$YEAR[df$JOURNAL == journals[3]],df$NUMBER[df$JOURNAL == journals[3]])

order1 <- order(df$YEAR, df$NUMBER_mod, df$PAGES1)
df[order1, c("JOURNAL", "YEAR", "NUMBER_mod", "PAGES1")]

## let's order using the months
# annals of stats publishes in: Feb, Apr, June, Aug, Oct, Dec (every 2 months)
# biometrika publishes in: March, June, Sept, Dec (every 3 months)
# JASA publishes in: March, June, Sept, Dec (every 3 months)
# JRSS-B publishes in: 
# Feb, May, Aug, Nov  (in 2003-2004)
# Feb, Apr, June, Sept, Nov  (2005-2007)
# Feb, Apr, July, Sept, Nov  (2008)
# Jan, Apr, June, Sept, Nov  (2009)
# Jan, March, June, Sept, Nov  (2010-2012)

df$MONTH <- NA

# AOS
subindex <- which(df$JOURNAL == journals[1])
months <- c(2,4,6,8,10,12)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

# BIO
subindex <- which(df$JOURNAL == journals[2])
months <- c(3,6,9,12)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

# JASA
subindex <- which(df$JOURNAL == journals[3])
months <- c(3,6,9,12)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

# JRSSB
subindex <- which((df$JOURNAL == journals[4]) &(df$YEAR %in% c(2003,2004)))
months <- c(2,5,8,11)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

subindex <- which((df$JOURNAL == journals[4]) &(df$YEAR %in% c(2005,2006,2007)))
months <- c(2,4,6,9,11)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

subindex <- which((df$JOURNAL == journals[4]) &(df$YEAR %in% c(2008)))
months <- c(2,4,7,9,11)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

subindex <- which((df$JOURNAL == journals[4]) &(df$YEAR %in% c(2009)))
months <- c(1,4,6,9,11)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

subindex <- which((df$JOURNAL == journals[4]) &(df$YEAR %in% c(2010,2011,2012)))
months <- c(1,3,6,9,11)
df$MONTH[subindex] <- months[df$NUMBER_mod[subindex]]

sum(is.na(df$MONTH))

order2 <- order(df$YEAR, df$MONTH, df$PAGES1, decreasing = TRUE)
head(df[order2, c("JOURNAL", "YEAR", "NUMBER", "PAGES1")]) #ordered dataset

### remember: rank 1 means more recent, rank k means oldest, so we used DECREASING
df$ranking <- match(1:nrow(df), order2) # this gives you the ordering.

bibtex_data <- df[, c(5,13,15,17,19,23,25,26,27,33,36,37,38,39)]

bibtex_data$N_CITATION_RECEIVED <- paperList$citCounts

bibtex_data[which.max(bibtex_data$ranking),] # 10.1214/aos/1046294456
save(bibtex_data, file = "data/bibtex_clean.rdata")



##########
## now that we have the order, we extract the citation data
##########

rm(list = ls())
load("data/bibtex_clean.rdata")
path = "data/SCC2016/"
paperList = read.table(paste0(path, "paperList.txt"), header = T, sep = ",")
paperCitAdj = read.table(paste0(path, "paperCitAdj.txt"), header = F, sep = " ")

tmp_order <- match(1:nrow(bibtex_data),bibtex_data$ranking)

Adj_ord <- paperCitAdj[tmp_order,tmp_order]
List_ord <- paperList[tmp_order,]

n_citing <- nrow(Adj_ord)

data <- data.frame(DOI_citing = NA, DOI_cited = NA, 
                   year_citing = NA, year_cited = NA,
                   tmp_order_citing = NA, tmp_order_cited = NA,
                   rank_citing = NA, rank_cited = NA)
count <- 1
for(i in 1:n_citing){
  # which papers is paper i citing? the columns give all the citations contained in that paper
  index_tmp <- which(Adj_ord[,i] == 1)
  for(j in index_tmp){
    data[count, ] <- c(List_ord[i,"DOI"],List_ord[j,"DOI"], 
                       List_ord[i,"year"],List_ord[j,"year"],
                       i, j,
                       bibtex_data$ranking[bibtex_data$DOI == List_ord[i,"DOI"]],
                       bibtex_data$ranking[bibtex_data$DOI == List_ord[j,"DOI"]])
    count <- count + 1
  }
}

# great, they are equal and correct.
# all.equal(data$tmp_order_citing, data$rank_citing)
# all.equal(data$tmp_order_cited, data$rank_cited)

data$year_citing <- as.numeric(data$year_citing)
data$year_cited <- as.numeric(data$year_cited)
data$tmp_order_citing <- as.numeric(data$tmp_order_citing)
data$tmp_order_cited <- as.numeric(data$tmp_order_cited)
data$rank_citing <- as.numeric(data$rank_citing)
data$rank_cited <- as.numeric(data$rank_cited)

ordered_data <- data

save(ordered_data, file = "data/ordered_cit_data.rdata")


DOI_oldest <- bibtex_data$DOI[which.max(bibtex_data$ranking)]
which(data$DOI_cited == DOI_oldest)
# the oldest paper is not cited:
sum(paperCitAdj[paperList$DOI == DOI_oldest,])

DOI_oldests <- bibtex_data$DOI[order(bibtex_data$ranking, decreasing = T)[1:5]]
rowSums(paperCitAdj[paperList$DOI %in% DOI_oldests,])
data[which(data$DOI_cited %in% DOI_oldests),]
