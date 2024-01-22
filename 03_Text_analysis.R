### Project for Quantitative Text Analysis

# (1) Preparations ------
# Clean work space
rm(list=ls())

# Set working directory
setwd("/Users/hannahfrank/Desktop/PhD/01_year/02/text_analysis/project")
getwd()

# Load packages
# https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
if(!require(quanteda)){
  install.packages("quanteda")
  library(quanteda)
}

if(!require(quanteda.textstats)){
  install.packages("quanteda.textstats")
  library(quanteda.textstats)
}

if(!require(quanteda.textplots)){
  install.packages("quanteda.textplots")
  library(quanteda.textplots)
}

if(!require(lubridate)){
  install.packages("lubridate")
  library(lubridate)
}

if(!require(stm)){
  install.packages("stm")
  library(stm)
}

if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}

if(!require(stargazer)){
  install.packages("stargazer")
  library(stargazer)
}

if(!require(AER)){
  install.packages("AER")
  library(AER)
}

# Load data
df <- read.csv("/Users/hannahfrank/Desktop/qca_peace_agreements_out/data/ucdp_full.csv",
                    stringsAsFactors=FALSE,
                    encoding = "utf-8")

# Fix date format
df$pa_date <- as.Date(df$pa_date, format = "%Y-%m-%d")

# Only include documents which are 'ready'
df_s <- df[df$ready==1, ]

# (2) Get corpus ----------

# Tokenize corpus
corpus_df <- corpus(df_s,
                    text_field = "fulltext")

# Summary statistics 
nrow(df_s)
summary(df_s$y_duration_months)

# Convert to lower case
toks <- quanteda::tokens(corpus_df,
                         include_docvars = TRUE) %>% 
                          tokens_tolower() 

# Remove tokens
toks <- quanteda::tokens(toks,
                         remove_numbers = TRUE,
                         remove_punct = TRUE,
                         remove_symbols = TRUE,
                         split_hyphens = FALSE,
                         remove_separators = TRUE,
                         remove_url = TRUE)
# Remove stop words
toks <- quanteda::tokens_remove(toks,
                                stopwords("english"),
                                padding = TRUE)
# Stem tokens
toks <- tokens_wordstem(toks)

toks <- tokens_select(toks, min_nchar=3L)

# Remove extra tokens manually
toks <- tokens_remove(toks, c("shall", 
                              "nshall",
                              "narticl",
                              "nbi",
                              "nand", 
                              "nof", 
                              "nin",
                              "nto", 
                              "also", 
                              "nfor", 
                              "nbe", 
                              "nthe",
                              "nwill", 
                              "nwith", 
                              "nthat", 
                              "nas", 
                              "nit"),  
                              valuetype = "fixed") 
# Find collocations
col <- textstat_collocations(toks,
                             size = 2, 
                             min_count = 10) 
# Sort collocations and check
col <- col[order(-col$count),]
head(col, 10)

# Choose cut-off at z>10 and add collocations
toks <- tokens_compound(toks,
                        pattern = col[col$z >10])
# Remove spaces
toks <- tokens_remove(quanteda::tokens(toks),
                      "")
# Create dfm
dfm <- dfm(toks)
dfm <- dfm_trim(dfm,
                min_docfreq = 20)

# Keep original full texts
dfm@docvars$original_text <- df_s$fulltext

# Keep durability, days and month
dfm@docvars$durability_days <- df_s$y_durability_days
dfm@docvars$durability_months <- df_s$y_durability_months

# Keep year
dfm@docvars$year <- df_s$year

# Print top 50 tokens
topfeatures(dfm, 50)

# Make word cloud for 500 top words
#pdf("corpus.pdf",
#    family="Times")
textplot_wordcloud(dfm, 
                   max_words = 500,
                   color="black")
#dev.off()

# Two empty documents which should not be the case
# Delete empty documents, this might be an incompatibility between R and Python
# <--- or an error TODO CHECK
row_sums<- apply(dfm, 1, sum)
dfm <- dfm[row_sums> 0, ] 

# Convert quanteda dfm to stm format
stm <- convert(dfm,
              to = "stm")

# (3) Implement topic model -------
# Find number of topics 
# https://juliasilge.com/blog/evaluating-stm/
kResult <- searchK(documents = stm$documents,
                   vocab = stm$vocab,
                   K=c(5,10,15,20,25,30,35,40),
                   init.type = "Spectral",
                   data = stm$meta,
                   seed = 2222,
                   prevalence = ~ s(as.numeric(year)),
                   cores = 5)
# Print results
kResult

# Plot results and store locally
#pdf("k.pdf",
#    family="Times")
par(mar=c(1,1,1,1))
plot(kResult)
#dev.off()

# Plot exclusivity versus semantic coherence
#pdf("k2.pdf",
#    height = 8,
#    width = 10,
#    family="Times")
par(mar=c(5,5,5,5))
plot(kResult$results$semcoh, 
     kResult$results$exclus, 
     cex=0.1,
     cex.axis = 1.8,
     cex.lab = 1.8, 
     xlab="Semantic Coherence",
     ylab="Exclusivity")
text(kResult$results$semcoh, kResult$results$exclus, 
     labels=kResult$results$K, cex=1.8)
#dev.off()
# This analysis suggests 15 or 20 topics, I will use 15

# Fit and save topic model with 15 topics
stm_fit <- stm(documents = stm$documents,
               vocab = stm$vocab,
               K = 15,
               prevalence = ~ s(as.numeric(year)),
               data = stm$meta,
               max.em.its = 500,
               init.type = "Spectral",
               seed = 1234,
               verbose = TRUE)

# Save
save(stm_fit,
     file = "stm_fit.Rdata")

# (4) Analyse topic model -------
# Print tokens with highest probability 
labelTopics(stm_fit)$prob[,1:7]

# Plot tokens with highest probability
par(mar=c(4,4,4,4))
#pdf("summary.pdf",
#    height = 8,
#    width = 10,
#    family="Times")
plot.STM(stm_fit,
         type = "summary",
         labeltype = "prob",
         text.cex = 1.8,
         cex.axis = 1.8,
         cex.lab = 1.8, 
         cex.main=1.8,
         xlab="Topic prevalence")
#dev.off()

# Make word cloud
par(mar=c(0.1,0.1,0.1,0.1))
#pdf("cloud_1.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 1, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Political power-sharing

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 1,
#             n = 1)

# Make word cloud
#pdf("cloud_2.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 2, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Demobilization

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 2,
#             n = 1)

# Make word cloud
#pdf("cloud_3.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 3, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Territory

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 3,
#             n = 1)

# Make word cloud
#pdf("cloud_4.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 4, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Government

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 4,
#             n = 1)

# Make word cloud
#pdf("cloud_5.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 5, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Parliament

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 5,
#             n = 1)

# Make word cloud
#pdf("cloud_6.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 6, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Nation-building

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 6,
#             n = 1)

# Make word cloud
#pdf("cloud_7.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 7, 
      scale = c(2.5, 0.3),
      max.words = 50)
#dev.off()
# Label: Dividing territory

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 7,
#             n = 1)

# Make word cloud
#pdf("cloud_8.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 8, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Transitional Government

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 8,
#             n = 1)

# Make word cloud
#pdf("cloud_9.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 9, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Governing cultural differences

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 9,
#             n = 1)

# Make word cloud
#pdf("cloud_10.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 10, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Procedural provisions

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 10,
#             n = 1)

# Make word cloud
#pdf("cloud_11.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 11, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Security sector reforms

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 11,
#             n = 1)

# Make word cloud
#pdf("cloud_12.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 12, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Involvement of United Nations

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 12,
#             n = 1)

# Search unit_nation
unit_nation <- kwic(corpus_df,
                    pattern = phrase("unit*"),
                    window = 15, 
                    case_insensitive = TRUE)
unit_nation

# Make word cloud
#pdf("cloud_13.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 13, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Security sector in more general terms

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 13,
#             n = 1)

# Make word cloud
#pdf("cloud_14.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 14, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Governing a state in more general terms

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 14,
#             n = 1)

# Make word cloud
#pdf("cloud_15.pdf",
#    family="Times")
cloud(stm_fit,
      topic = 15, 
      scale = c(2.5, 0.3),
      max.words = 50) 
#dev.off()
# Label: Reconciliation

# Look at the original documents 
#findThoughts(stm_fit,
#             texts = dfm@docvars$original_text,
#             topics = 15,
#             n = 1)

# Search victim
victim <- kwic(corpus_df,
                pattern = phrase("victim"),
                window = 15,
                case_insensitive = TRUE)
victim

### (5) Validation -------
# Make plot for semantic validity
topic_correlations <- topicCorr(stm_fit) 
par(mar=c(4,4,4,4))
#pdf("corr.pdf",
#    family="Times")
plot.topicCorr(topic_correlations,
               vlabels = seq(1:ncol(stm_fit$theta)),
               vertex.color = "white",
               main = "Topic correlations")
#dev.off()

# Make plot for coherence versus exclusivity
#pdf("versus.pdf",
#    family="Times")
par(mar=c(4,4,4,4))
topicQuality(model = stm_fit,
             documents = stm$documents,
             xlab = "Semantic Coherence",
             ylab = "Exclusivity",
             labels = 1:ncol(stm_fit$theta),
             M = 15) 
#dev.off()

### (6) Predict the durability of peace -------
# Get data
df_model <- data.frame(cbind(stm$meta$y_duration_months, 
                             stm$meta$pol.prov,
                             stm$meta$mil.prov,
                             stm$meta$incompatibility,
                             stm$meta$active_conflict,
                             stm$meta$out_iss,
                             stm_fit$theta[,1],
                             stm_fit$theta[,2],
                             stm_fit$theta[,12],
                             stm_fit$theta[,15],
                             stm_fit$theta[,3],
                             stm_fit$theta[,8],
                             stm_fit$theta[,10],
                             stm_fit$theta[,7],
                             stm_fit$theta[,11],
                             stm_fit$theta[,4],
                             stm_fit$theta[,5],
                             stm_fit$theta[,13],
                             stm_fit$theta[,9],
                             stm_fit$theta[,6],
                             stm_fit$theta[,14]
                             ))

# Rename columns
colnames(df_model) <- c("y_duration_months", 
                        "pol_prov",
                        "mil_prov",
                        "incompatibility",
                        "active_conflict",
                        "out_iss",
                        "topic_1",
                        "topic_2", 
                        "topic_12",
                        "topic_15",
                        "topic_3", 
                        "topic_8",
                        "topic_10",
                        "topic_7",
                        "topic_11",
                        "topic_4",
                        "topic_5",
                        "topic_13",
                        "topic_9",
                        "topic_6",
                        "topic_14")  

# Convert topic prevalence as numeric
df_model$y_duration_months <- as.numeric(df_model$y_duration_months)    
df_model$topic_1 <- as.numeric(df_model$topic_1)    
df_model$topic_2 <- as.numeric(df_model$topic_2)                                 
df_model$topic_12 <- as.numeric(df_model$topic_12)                                 
df_model$topic_15 <- as.numeric(df_model$topic_15)                                 
df_model$topic_3 <- as.numeric(df_model$topic_3)
df_model$topic_8 <- as.numeric(df_model$topic_8)  
df_model$topic_10 <- as.numeric(df_model$topic_10)   
df_model$topic_7 <- as.numeric(df_model$topic_7)
df_model$topic_11 <- as.numeric(df_model$topic_11)
df_model$topic_4 <- as.numeric(df_model$topic_4)  
df_model$topic_5 <- as.numeric(df_model$topic_5)   
df_model$topic_13 <- as.numeric(df_model$topic_13)   
df_model$topic_9 <- as.numeric(df_model$topic_9)   
df_model$topic_6 <- as.numeric(df_model$topic_6)   
df_model$topic_14 <- as.numeric(df_model$topic_14)   

# Make additive index 
df_model$strength <- df_model$topic_1+
  df_model$topic_2+
  df_model$topic_12+
  df_model$topic_15+
  df_model$topic_3

# Predict topic prevalence with UCDP coding
val_1 <- estimateEffect(c(1) ~ pol.prov,
                          stm_fit,
                          metadata = stm$meta,
                          uncertainty = "Global",
                          nsims = 25)
# Make effect plot
summary(val_1)
par(mar=c(4,2,2,2))
#pdf("val_1_effect.pdf",
#    height = 5,
#    width = 9,
#    family="Times")
plot(val_1, 
     "pol.prov",
     labeltype = "custom",
     custom.labels=c("Political Provisions", 
                     "No Political Provisions"),
     xlim = c(-0.02, .3),
     cex=1.8,
     cex.axis = 1.8,
     cex.lab = 1.8)
#dev.off()

# Predict topic prevalence with UCDP coding
val_2 <- estimateEffect(c(2) ~ mil.prov,
                        stm_fit,
                        metadata = stm$meta,
                        uncertainty = "Global",
                        nsims = 25)
# Make effect plot
summary(val_2)
par(mar=c(1,1,1,1))
#pdf("val_2_effect.pdf",
#    height = 5,
#    width = 9,
#    family="Times")
plot(val_2, 
     "mil.prov",
     labeltype = "custom",
     custom.labels=c("No Military Provisions", 
                     "Military Provisions"),
     xlim = c(-0.02, .18),
     cex=1.8,
     cex.axis = 1.8,
     cex.lab = 1.8)
#dev.off()

# Poisson regression
# http://biometry.github.io/APES//LectureNotes/2016-JAGS/Overdispersion/OverdispersionJAGS.html
m1_poission <- glm(y_duration_months ~ 
        topic_1+
        topic_2+  
        topic_12+
        topic_15+
        topic_3+
        topic_8+
        topic_10+
        topic_7+
        topic_11+
        topic_4+
        topic_5+
        topic_13+
        topic_9+
        topic_6,
       data = df_model,
       family = poisson)
summary(m1_poission)
dispersiontest(m1_poission) # overdispersion

# Run negative binomial count model 
# All topics, topic 14 is NA
summary(m1 <- glm.nb(y_duration_months ~ 
                       topic_1+
                       topic_2+
                       topic_12+
                       topic_15+
                       topic_3+
                       topic_8+
                       topic_10+
                       topic_7+
                       topic_11+
                       topic_4+
                       topic_5+
                       topic_13+
                       topic_9+
                       topic_6,
                      # topic_14,
                     data = df_model))

# Run negative binomial count model
# Most prevalent topics
summary(m2 <- glm.nb(y_duration_months ~ 
                       topic_1+
                       topic_2+
                       topic_12+
                       topic_15+
                       topic_3+
                       factor(incompatibility)+
                       factor(active_conflict),
                       data = df_model))

# Run negative binomial count model
# with significant topics
summary(m3 <- glm.nb(y_duration_months ~ 
                       topic_2+
                       topic_12+
                       topic_3+
                       topic_10+
                       topic_7+
                       factor(incompatibility)+
                       factor(active_conflict),
                     data = df_model))

# Run negative binomial count model
# with index
summary(m4 <- glm.nb(y_duration_months ~ 
                       strength+
                       factor(incompatibility)+
                       factor(active_conflict),
                       data = df_model))

# Get regression table
stargazer(m1,
          m4, 
          m3)



