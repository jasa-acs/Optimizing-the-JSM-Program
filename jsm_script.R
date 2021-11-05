## R 3.6.1
## packages dependencies:
## text2vec 0.5.1, data.table 1.13.4, stringr 1.4.0, stopwords 1.0, corpus 0.10.0,
## slam 0.1-45, topicmodels 0.2-8, ggplot2 3.3.2

library(text2vec)
library(data.table)

rm(list = ls())

# loading of trayvon_martin, abstract_info, seeded_words
load("inputdata_JSM.RData")

#--------------------------------------------------
# Settings
#--------------------------------------------------

# Pre-processing settings
pvalue_thr <- 0.005
min_cooccurences_words <- 5L
min_cooccurences_ngrams <- 20L
max_ngram <- 5L
greedy_n <- 100L

# Topic Modeling settings
SeedWeight <- 100L
burn_par <- 500L
num.iterations_par <- 400L
alpha_par <- 1
eta_par <- 1L
seed <- 1L
extra_topics <- 10L

# Optimization settings
gamma_thr <- 0.2
max_iterations_greedy <- 1e7 # max number of iterations of while loop for stopping
limit_no_change_greedy <- 1e4 # number of no change iterations for stopping



#--------------------------------------------------
# Functions
#--------------------------------------------------

# total variation distance
tv_distance <- function(X,Y){max(abs(X-Y))}

# rho
rho_index <- function(allocation_matrix,distance_matrix){
  rho <- c()
  for(i in 1:nrow(allocation_matrix)){
    ind <- allocation_matrix[i,][allocation_matrix[i,] != "SNA"]
    rho <- c(rho,sum(distance_matrix[ind,ind]))
  }
  return(sum(rho)/2)
}

check_incomp <- function(allocation_matrix, incomp_matrix){
  inc <- c()
  for(i in 1:nrow(allocation_matrix)){
    ind <- allocation_matrix[i,][allocation_matrix[i,] != "SNA"]
    mt <- incomp_matrix[ind,ind]
    diag(mt) <- FALSE
    inc <- c(inc,any(mt))
  }
  return(any(inc))
}


denoise <- function(X, gamma_thr){
  if(max(X) > gamma_thr){
    X[X<gamma_thr] <- 0
  } else {
    X[-which.max(X)] <- 0
  }
  X <- X/sum(X)
  return(X)
}

combination_matrix <- function(matrix, pair_coo, pair_label, fun, new_value){
  new_matrix <- matrix[-pair_coo,-pair_coo]
  newrow <- switch(fun,
                   "any" = apply(matrix[pair_coo,-pair_coo],2,any),
                   "mean" = apply(matrix[pair_coo,-pair_coo],2,mean))
  new_matrix <- rbind(new_matrix,newrow)
  new_matrix <- cbind(new_matrix,c(newrow,new_value))
  rownames(new_matrix)[nrow(new_matrix)] <- pair_label
  colnames(new_matrix)[ncol(new_matrix)] <- pair_label
  return(new_matrix)
}

remotion_stopwords <- function(text,stopwords){
  words <- unlist(strsplit(text, " "))
  cleaned_words <- words[sapply(words, function(X) !X %in% stopwords)]
  cleaned_words <- paste(cleaned_words, collapse = " ")
  return(cleaned_words)
}

distinctive_words <- function(lda_results){
  # P[ j | word ] = P[ word | j ] * P[ j ]/ sum_{k=1}^K P[ word | k ] * P[ k ].
  # with P[ j ] = sum_{k=1}^K ( P[ j | session_k ] * P[ session_k ] ), where
  # - P[ session_k ] = 1/(# of sessions)
  # - P[ j | session_k ] is the output of LDA
  
  gamma <- lda_results@gamma
  beta <- lda_results@beta
  terms <- lda_results@terms
  p_s <- 1/nrow(gamma)
  p_j <- t(rep(1,nrow(gamma)) %*% (gamma*p_s))
  num <- exp(beta) * (p_j %*% t(rep(1,ncol(beta))))
  p_jw <-  num / ( rep(1,nrow(beta)) %*% (rep(1,nrow(beta)) %*% num ))
  colnames(p_jw) <- terms
  return(p_jw)
}


#--------------------------------------------------
# Pre-processing
#--------------------------------------------------


# the JSM program
prog <- abstract_info[which(!to_exclude),list(start,end,session_id)]
prog <- prog[!duplicated(prog)]
prog <- prog[order(prog[,start])]
slots_prog <- prog[,.N,by=list(start,end)]

time_bands <- nrow(slots_prog)
max_parallels <- max(slots_prog$N)
allocation_matrix_empty <- matrix(NA_character_, time_bands, max_parallels)

for(i in 1:time_bands){
  if((slots_prog[i,N]+1) <= max_parallels){
    allocation_matrix_empty[i,(slots_prog[i,N]+1):max_parallels] <- "SNA"
  }
}

# escluding seeded words from trayvon_martin
trayvon_martin <- trayvon_martin[!trayvon_martin %in% seeded_words[,word]]

# creation of keywords2
abstract_info[,keywords2:=sapply(keywords, function(X)
  paste(gsub(" ","-",unlist(strsplit(X, "; "))), collapse = " "))]

# sessions text creation: abstracts of the same sessions are merged
sessions <- abstract_info[which(!to_exclude),
                          list(text=paste(abstract_title,abstract_title,abstract_title,
                                          keywords2,keywords2,keywords2,
                                          abstract,collapse = " "),
                               speaker = paste(speaker,collapse = " and ")), by=session_id]


list_speakers <- strsplit(sessions$speaker, " and ")
names(list_speakers) <- sessions$session_id

incomp_matrix <- matrix(FALSE,length(list_speakers),length(list_speakers))
rownames(incomp_matrix) <- colnames(incomp_matrix) <- sessions$session_id
for(i in sessions$session_id){
  speakers_ses <- list_speakers[[i]]
  incomp <- c()
  for(j in speakers_ses){
    incomp <- c(incomp,sessions$session_id[which(sapply(list_speakers, function(X) j %in% X))])
  }
  incomp <- unique(incomp)
  incomp_matrix[i,incomp] <- TRUE
}

#########################
# text cleaning
#########################

ses_clean <- gsub("[[:digit:]]","",sessions$text) # remotion of: numbers
ses_clean <- gsub('[“”‘’]'," ", ses_clean) # remotion of: special characters
ses_clean <- gsub('\"'," ", ses_clean) # '\"' into space
ses_clean <- gsub('–+',"1", ses_clean) # dashes into '1'
ses_clean <- gsub('-+',"1", ses_clean) # dashes into '1'
ses_clean <- gsub('[[:punct:]]'," ", ses_clean) # remotion of: punctuations
ses_clean <- gsub('1',"_", ses_clean) # '1' into '_'
ses_clean <- gsub("\\s+", " ", stringr::str_trim(ses_clean)) # remotion of: double space
ses_clean <- gsub(" _ ", " ", stringr::str_trim(ses_clean)) # remotion of: " _ "
ses_clean <- gsub(" _", "_", ses_clean) # remotion of space
ses_clean <- gsub("_ ", "_", ses_clean) # remotion of space
ses_clean <- tolower(ses_clean) # text into lowercase

# definition of the stop words from two sources ('snowball' and 'stopwords-iso')
stop_words <- c(stopwords::stopwords("en", source = "snowball"),
                stopwords::stopwords("en", source = "stopwords-iso"))

# remotion of stop words
ses_clean <- lapply(ses_clean, function(X) remotion_stopwords(X, stop_words))

step1 <- paste(ses_clean, collapse = " ")
(step1_bis <- length(unique(unlist(strsplit(step1, " ")))))

#-------------------------------------------------------
# At the end of step 1 (remotion of stop words):
# 21813 words
#-------------------------------------------------------

# lemmatization of the words
ses_clean <- unlist(lapply(ses_clean, function(X) paste(corpus::stem_snowball(unlist(strsplit(X, " "))), collapse = " ")))

step2 <- paste(ses_clean, collapse = " ")
step2_words <- unique(unlist(strsplit(step2, " ")))
(step2_length <- length(step2_words))
step2_ngrams <- sum(sapply(strsplit(step2_words,"_"), function(X) length(X)>1))

# tokenization of the text
tokens <- itoken(ses_clean,
                 tokenizer = word_tokenizer,
                 ids = sessions$session_id,
                 progressbar = FALSE)

vocab_step2 <- create_vocabulary(tokens)

pruned_vocab_step2 <- prune_vocabulary(vocab_step2
                                       ,doc_count_min = 2
                                       ,term_count_min = min_cooccurences_words
)

nrow(pruned_vocab_step2)

#-------------------------------------------------------
# At the end of step 2 (lemmatization of the words):
# 16331 words (of which 5264 n-grams)
#-------------------------------------------------------


# creation of a single vocabulary for each session
global_vocab1 <- data.table()

for(ss in 1:length(ses_clean)){
  tokens1 <- itoken(ses_clean[ss],
                    tokenizer = word_tokenizer,
                    ids = sessions$session_id[ss],
                    progressbar = FALSE)

  vocab1 <- create_vocabulary(tokens1#, stopwords = stop_words
  )

  vocab1 <- as.data.frame(vocab1[,c("term", "term_count")])
  vocab1$session <- ss

  global_vocab1 <- rbind(global_vocab1,vocab1)
}


# vocabulary creation for single words
vocab <- create_vocabulary(
  tokens, ngram = c(1L, max_ngram)
  #,stopwords = stop_words
)

pruned_vocab <- prune_vocabulary(vocab
                                 ,doc_count_min = 2
                                 ,term_count_min = min_cooccurences_words
)

ngrams_position <- sapply(strsplit(pruned_vocab$term,"_"), function(X) length(X)>1)


unigram_vocab <- data.table(pruned_vocab[!ngrams_position,])
ngram25_vocab <- data.table(pruned_vocab[ngrams_position,])
ngram25_vocab[,pvalue := NA_real_]

for(ngram in ngram25_vocab[,term]){
  terms <- unlist(strsplit(ngram,"_"))
  n <- length(unlist(strsplit(terms,"_")))
  count25 <- ngram25_vocab[ngram == term,term_count]

  terms_count <- global_vocab1[term %in% terms
                               # & session %in% common_sessions
                               ,
                               list(count=sum(term_count)),by=term]

  sessions_count <- global_vocab1[#session %in% common_sessions
    ,list(len_ses=length(term)), by=session]
  n_sessions <- nrow(sessions_count) #length(common_sessions)
  len1 <- sum(sessions_count$len_ses)

  len25 <- len1 - n*n_sessions + n_sessions

  ngram25_vocab[term == ngram, pvalue :=
                  prop.test(c(count25,prod(terms_count$count)),c(len25,len1^n),
                            alternative = "greater")$p.value]
}


ngrams_list <- pruned_vocab[ngrams_position,"term"]
ngrams2remove <- ngrams_list[!ngrams_list %in% ngram25_vocab[pvalue <= pvalue_thr, term]]
pruned_vocab <- pruned_vocab[!pruned_vocab$term %in% ngrams2remove,]

nrow(pruned_vocab)

(step3_ngrams <- sum(sapply(strsplit(pruned_vocab$term,"_"), function(X) length(X)>1)))
#-------------------------------------------------------
# At the end of step 3 (n-gramming):
# 17471 words (of which 14017 n-grams)
#-------------------------------------------------------

# remotion from vocab_1 trayvon martin corpus and 'statistician','statist'
pruned_vocab <- pruned_vocab[!pruned_vocab$term %in% c(trayvon_martin,"statistician","statist","_statist"),]

nrow(pruned_vocab)

#-------------------------------------------------------
# At the end of step 4 (Trayvon Martin corpus):
# 15755 words
#-------------------------------------------------------



#--------------------------------------------------
# Topic Modeling
#--------------------------------------------------


vectorizer <- vocab_vectorizer(pruned_vocab)

# document-term matrix construction and its coertion into data.frame class
dtm <- create_dtm(tokens, vectorizer)

ii <- jj <- c()

for(w in 1:nrow(seeded_words)){
  word <- seeded_words[w,]
  sw <- which(dtm@Dimnames[[2]] == word[,word])
  ii <- c(ii,rep(as.integer(word[,topic]),length(sw)))
  jj <- c(jj,sw)
}

ordering <- order(ii)

jj <- jj[ordering]
ii <- ii[ordering]

K <- length(unique(ii)) + extra_topics

deltaS <- slam::simple_triplet_matrix(ii, jj, v = rep(SeedWeight, length(ii)),
                                      nrow = K, ncol = ncol(dtm))

lda_results <- topicmodels::LDA(dtm, k = K, method = "Gibbs", seedwords = deltaS,
                                control = list(alpha = alpha_par, best = TRUE, seed = seed,
                                               verbose = 500, burnin = burn_par, iter = 100,
                                               thin = 100, prefix = character()))



#--------------------------------------------------
# Optimization
#--------------------------------------------------



# creation of the actual allocation matrix of JSM
actual_allocation_matrix <- t(allocation_matrix_empty)
actual_allocation_matrix[is.na(actual_allocation_matrix)] <- prog[,session_id]
actual_allocation_matrix <- t(actual_allocation_matrix)

# creation of a random allocation matrix of JSM without overlapping
random_allocation_matrix_list <- list()
counter <- 1
condition_random_overlap <- TRUE
while(condition_random_overlap){
  set.seed(counter)
  random_allocation_matrix <- allocation_matrix_empty
  random_allocation_matrix[is.na(random_allocation_matrix)] <- sample(prog[,session_id])
  no_overlap <- !check_incomp(random_allocation_matrix, incomp_matrix)
  if(no_overlap){
    random_allocation_matrix_list <- c(random_allocation_matrix_list,list(random_allocation_matrix))
  }
  condition_random_overlap <- length(random_allocation_matrix_list) < 100
  counter <- counter + 1
}



# creation of our allocation 

# gamma matrix
gamma_matrix <- t(apply(lda_results@gamma,1, function(X) denoise(X, gamma_thr)))

# id of sessions
id_sessions <- lda_results@documents

# number of sessions
n_sessions <- length(id_sessions)

distance_matrix <- matrix(NA_real_,n_sessions,n_sessions)
colnames(distance_matrix) <- rownames(distance_matrix) <- id_sessions
diag(distance_matrix) <- 0

# definition of distance matrix through total variation (hence TV) metric
for(i in id_sessions){
  for(j in id_sessions){
    if(which(i == id_sessions) < which(j == id_sessions)){
      distance_matrix[j,i] <- distance_matrix[i,j] <-
        tv_distance(gamma_matrix[id_sessions == i,], gamma_matrix[id_sessions == j,])
    }
  }
}

# The result of our allocation is in "allocation_matrix_list.RData"
# to save time it is possible to skip the lines 407-572 and load
# directly that. To run the allocation by your own, please change
# FALSE to TRUE in "if" statement below.

# 
if(FALSE) {
  # sessions indexes with speaker overlapping
  overlap_speakers <- which(apply(incomp_matrix, 1, sum) > 1)

  # incomp_matrix with at least one overlap
  incomp_matrix_overlap <- incomp_matrix[overlap_speakers,overlap_speakers]

  # distance_matrix with at least one overlap
  distance_matrix_overlap <- distance_matrix[overlap_speakers,overlap_speakers]

  # creation of an empty data.frame to collect sessions overlapping aggregation info
  aggregation <- data.frame(step = numeric(0), overlaps = numeric(0), pairs = character(0))

  iteration_aggregation <- 1
  condition_aggregation <- TRUE

  while(condition_aggregation){

    # setting lower triangle of incomp_matrix_overlap as TRUE
    incomp_matrix_overlap_lower <- incomp_matrix_overlap
    incomp_matrix_overlap_lower[upper.tri(incomp_matrix_overlap_lower)] <- TRUE

    # looking for not overlapping sessions to put in the same time slot
    pairs_ind <- pairs_id <- as.data.frame(which(!incomp_matrix_overlap_lower, arr.ind=TRUE), row.names = FALSE)

    # replacing session index with session id
    pairs_id[,"row"] <- rownames(incomp_matrix_overlap)[pairs_ind[,"row"]]
    pairs_id[,"col"] <- colnames(incomp_matrix_overlap)[pairs_ind[,"col"]]

    # adding of the distance between sessions in each pair
    for(i in 1:nrow(pairs_id)){pairs_id[i,"dist"] <- distance_matrix_overlap[pairs_id[i,"row"],pairs_id[i,"col"]]}

    pair_max <- pairs_id[which.max(pairs_id$dist),c("row","col")] # pair with maximum distance to combine
    pair_label <- paste(pair_max, collapse = "_") # label of the pair to combine
    pair_coo <- which(rownames(incomp_matrix_overlap) %in% pair_max) # indexes of the pair to combine

    # remotion, for incomp_matrix_overlap and distance_matrix_overlap,
    # of the info on the sessions of the pair to combine
    # and addition of the info of the combined pair
    incomp_matrix_overlap <- combination_matrix(
      incomp_matrix_overlap, pair_coo, pair_label, fun = "any", new_value = TRUE)

    distance_matrix_overlap <- combination_matrix(
      distance_matrix_overlap, pair_coo, pair_label, fun = "mean", new_value = 0)

    # inclusion of the sessions overlapping aggregation info
    aggregation <- rbind(aggregation,
                       data.frame(step=iteration_aggregation,overlaps = nrow(incomp_matrix_overlap),
                                  pairs=paste(rownames(incomp_matrix_overlap), collapse = "__"),
                                  stringsAsFactors=FALSE))

    # checking if no overlapping sessions still exist
    condition_aggregation <- any(!incomp_matrix_overlap[upper.tri(incomp_matrix_overlap)])

    iteration_aggregation <- iteration_aggregation + 1
  }

  # selection of sessions combination suitable for the number of slot time
  selected_aggregation <- aggregation[aggregation$overlaps <= nrow(allocation_matrix_empty),]

  # definition of the best aggregation
  rho_overlap <- c()
  for(p in selected_aggregation$pairs){
    am <- allocation_matrix_empty
    pairs <- unlist(strsplit(p, "__"))
    pairs2 <- lapply(pairs, function(X) unlist(strsplit(X, "_")))
    n <- length(pairs2)
    for(i in 1:n){am[i,1:length(pairs2[[i]])] <- pairs2[[i]]}
    am[is.na(am)] <- "SNA"
    rho_overlap <- c(rho_overlap,rho_index(am,distance_matrix)/n)
  }

  # selection of the best aggregation
  best_aggregation <- unlist(strsplit(selected_aggregation$pairs[which.max(rho_overlap)],"__"))
  best_aggregation_list <- lapply(best_aggregation, function(X) unlist(strsplit(X, "_")))

  # allocation of the best aggregation to allocation_matrix
  allocation_matrix_starting <- allocation_matrix_empty

  for(i in 1:length(best_aggregation_list)){
    allocation_matrix_starting[i,1:length(best_aggregation_list[[i]])] <- best_aggregation_list[[i]]}

  # definition of the sessions with no overlapping speakers (hence 'free sessions')
  # that can be used in the greedy optimization approach
  free_session_id <- as.vector(actual_allocation_matrix)
  free_session_id <- free_session_id[!free_session_id %in% c(unlist(pairs2),"SNA")]

  allocation_matrix_list <- iteration_greedy_list <- as.list(rep(NA,greedy_n))

  for(u in 1:greedy_n){

    am <- allocation_matrix_starting

    # random allocation of the free sessions to the residual empty cells of am
    set.seed(u)
    am[is.na(am)] <- sample(free_session_id)

    ## Greedy optimization approach

    new_allocation_matrix <- am # creation of a copy of the original allocation
    changing <- c() # creation of vector to store step changing

    # rho index at starting point
    starting_rho <- rho_index(new_allocation_matrix,distance_matrix)

    # counters
    iteration_greedy <- 1
    no_change <- 0

    # starting condition of the loop
    condition_greedy <- TRUE

    while(condition_greedy){

      set.seed(iteration_greedy-1)
      sampling <- sample(free_session_id, 2) # sampling of two sessions
  
      # the random second session
      set.seed(iteration_greedy)
      sampling <- sample(free_session_id, 2) # sampling of two sessions

      ind1 <- which(new_allocation_matrix == sampling[1])
      ind2 <- which(new_allocation_matrix == sampling[2])

      new_allocation_matrix_copy <- new_allocation_matrix

      # to change the position of the sessions
      new_allocation_matrix_copy[ind1] <- sampling[2]
      new_allocation_matrix_copy[ind2] <- sampling[1]

      if(check_incomp(new_allocation_matrix_copy, incomp_matrix))
        stop("New allocation has speakers overlapping")

      # to compute rho index of the new_allocation_matrix_copy
      rho_new <- rho_index(new_allocation_matrix_copy,distance_matrix)
      rho_previous <- rho_index(new_allocation_matrix,distance_matrix)

      if(rho_new > rho_previous){
        changing <- c(changing,iteration_greedy)
        new_allocation_matrix <- new_allocation_matrix_copy
        no_change <- 0
        gain_from_start <- rho_new - starting_rho
      } else {
        no_change <- no_change + 1
      }

      if(no_change == limit_no_change_greedy | iteration_greedy == max_iterations_greedy){
        condition_greedy <- FALSE
      } else {
        condition_greedy <- TRUE
      }
      iteration_greedy <- iteration_greedy + 1
    }

    allocation_matrix_list[[u]] <- new_allocation_matrix
    iteration_greedy_list[[u]] <- iteration_greedy
  }

  save(allocation_matrix_list, iteration_greedy_list,
      file = "allocation_matrix_list.RData")
}



load("allocation_matrix_list.RData")



rho_new <- sapply(allocation_matrix_list, function(X) rho_index(X,distance_matrix))

# the rho index is computed by summing the row distances: the goal is to maximized it
rho_random <- mean(sapply(random_allocation_matrix_list, function(X) rho_index(X, distance_matrix)))
rho_actual <- rho_index(actual_allocation_matrix,distance_matrix) # 10124.76

# definition of the theoretical maximum
parallels <- apply(allocation_matrix_empty,1,function(X) sum(is.na(X)))
rho_max <- sum(parallels*(parallels-1)/2)

# efficiency index (with range [0,1])
(rho_actual-rho_random)/(rho_max-rho_random) # JSM allocation
(max(rho_new)-rho_random)/(rho_max-rho_random) # our approach allocation

# incompatibility check: FALSE -> no overlapping!
any(sapply(random_allocation_matrix_list, function(X) check_incomp(X, incomp_matrix)))
any(sapply(allocation_matrix_list, function(X) check_incomp(X, incomp_matrix)))
check_incomp(actual_allocation_matrix, incomp_matrix)




#--------------------------------------------------
# Paper outputs
#--------------------------------------------------

dw <- distinctive_words(lda_results)

rt <- lapply(1:nrow(dw), function(X) sort(dw[X,],decreasing = TRUE)[1:5])

tabs_2_3 <- data.frame(a=1:5)
for(i in 1:length(rt)){
  df <- data.frame(prob=rt[[i]])
  df$token <- rownames(df)
  rownames(df) <- NULL
  df <- df[c("token","prob")]
  names(df)[1] <- paste0("token ","topic",i)
  
  tabs_2_3 <- cbind(tabs_2_3,df)
}
tabs_2_3$a <- NULL

tab_2 <- tabs_2_3[,1:82]
tab_3 <- tabs_2_3[,83:102]




allocation_new <- allocation_matrix_list[[which.max(rho_new)]]
#actual_allocation_matrix

tv_allocation_newS <- tv_allocation_newA <- n_sessions_new <- c()
for(i in 1:nrow(allocation_new)){
  ind <- allocation_new[i,][allocation_new[i,] != "SNA"]
  mm <- distance_matrix[ind,ind]
  n_sessions_new <- c(n_sessions_new,length(ind))
  tv_allocation_newA <- c(tv_allocation_newA,mean(mm[upper.tri(mm)]))
  tv_allocation_newS <- c(tv_allocation_newS,sum(mm[upper.tri(mm)]))
}

tv_allocation_actualS <- tv_allocation_actualA <- n_sessions_actual <- c()
for(i in 1:nrow(actual_allocation_matrix)){
  ind <- actual_allocation_matrix[i,][actual_allocation_matrix[i,] != "SNA"]
  mm <- distance_matrix[ind,ind]
  n_sessions_actual <- c(n_sessions_actual,length(ind))
  tv_allocation_actualA <- c(tv_allocation_actualA,mean(mm[upper.tri(mm)]))
  tv_allocation_actualS <- c(tv_allocation_actualS,sum(mm[upper.tri(mm)]))
}
identical(n_sessions_new,n_sessions_actual)


time_band_rho <- data.frame(Time_band = 1:13, n_parallel = n_sessions_new, 
                            Sum_new = format(round(tv_allocation_newS,2), nsmall = 2),
                            Avg_new = format(round(tv_allocation_newA,3), nsmall = 3),
                            Sum_actual =  format(round(tv_allocation_actualS,2), nsmall = 2),
                            Avg_actual = format(round(tv_allocation_actualA,3), nsmall = 3))

names(time_band_rho) <- c("Time band","# of parallel sessions", "Best assignment Sum",
                          "Best assignment Average", "JSM 2020Sum", "JSM 2020 Average")





# Table 1 is: seeded_words
# Table 2 is: tab_2
# Table 3 is: tab_3
# Table 5 is: time_band_rho

# Figure 1
library(ggplot2)

filename = paste('histogram.eps')
setEPS()
postscript(filename,height=8,width=10)

ggplot(data.frame(rho=rho_new), aes(x=rho)) + geom_histogram(color="black", fill="white", bins=15) + 
  theme_classic()+ theme(plot.title = element_text(size=30),
                         axis.title.x = element_text(size = 20),
                         axis.title.y = element_text(size = 20),
                         axis.text.x = element_text(size = 16),
                         axis.text.y = element_text(size = 16)) +  
  labs( x = expression(rho))

dev.off()







