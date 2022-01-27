rm(list = ls())
setwd("~/Dropbox (Gladstone)/GB-SM-1238/CVS time motif info")
require(dplyr)
require(magrittr)
require(markovchain)
data <- read.csv("210416 motif_time 25Hz v5.csv", header = TRUE)

mouse_ids <- data %>%
  .$Mouse %>%
  unique() %>%
  sort()

motifs <- data %>%
  .$Motif %>%
  unique() %>%
  sort()

pheno_data <- data %>%
  dplyr::select(Mouse,Genotype, Sex) %>%
  unique() %>%
  arrange(Mouse)


##load community mapping
motif_communities <- "MotifCommunities.csv" %>%
  read.csv(., header = TRUE)
data %<>% merge(., motif_communities)

mid <- 13
state_transition <- c(0,11)
resolution <- 0.04*16
time_step <- resolution/0.04

mode <- function(comm){
  names(which.max(table(comm)))
}

get_mouse_data <- function(mid, data, time_step) {
  mouse_data <- data %>%
    filter(Mouse == mid) %>%
    arrange(Seconds)
  
  rows_to_use <- seq(1, nrow(mouse_data), by=time_step)
  group_time <- rep(1:length(rows_to_use), each=time_step) %>%
    as.factor() 
  
  group_time <- group_time[1:nrow(mouse_data)]
  
  sub_community <- mouse_data %>% 
    slice(1:length(group_time)) %>% 
    mutate(group_time=group_time) %>%
    group_by(group_time) %>%
    summarise(mode = mode(SubCommunity))
  
  mouse_data %<>% slice(rows_to_use) %>%
    mutate(Motif=sub_community$mode)
  
  return(mouse_data)
}

get_transition_matrix <- function(mid, data, pheno_data, time_step) {
  mouse_data <- get_mouse_data(mid,data,time_step)

  mcFitMLE <- mouse_data %>%
    dplyr::slice(1:nrow(.)) %>% 
    .$Motif %>%
    as.character() %>%
    markovchainFit(data = .)
  
  # mouse_mc <- new("markovchain", states=as.character(motifs), transitionMatrix=mcFitMLE$estimate@transitionMatrix)
  transitionMatrix=mcFitMLE$estimate@transitionMatrix
  return(transitionMatrix)
}

get_motif_counts <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  counts <- mouse_data %>%
    .$Motif %>%
    table() 
  return(counts)
}

get_motif_transition_counts <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  sequence <- mouse_data %>%
    .$Motif %>%
    as.character()
  
  transMatrix <- markovchainFit(data = sequence)$estimate@transitionMatrix
  subSample <- sequence[1:(length(sequence) - (length(sequence)%%3))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 2)], 
                      subSample[2:(length(subSample) - 1)], subSample[3:(length(subSample))]), 
                    ncol = 3)
  temp <- as.data.frame(seqSet1)
  Nijk <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  seqSet2 <- seqSet1[, -3]
  temp2 <- as.data.frame(seqSet2)
  Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)]
  
  Nij %<>% as.data.frame()
  colnames(Nij) <- c("start_motif", "end_motif", "counts")
  return(Nij)
}

get_second_order_motif_transition_counts <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  sequence <- mouse_data %>%
    .$Motif %>%
    as.character()
  
  # transMatrix <- markovchainFit(data = sequence)$estimate@transitionMatrix
  # subSample <- sequence[1:(length(sequence) - (length(sequence)%%4))]
  # seqSet1 <- matrix(c(subSample[1:(length(subSample) - 3)], 
  #                     subSample[2:(length(subSample) - 2)], 
  #                     subSample[3:(length(subSample) - 1)], subSample[4:(length(subSample))]), 
  #                   ncol = 4)
  # temp <- as.data.frame(seqSet1)
  # Nijkl <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  
  subSample <- sequence[1:(length(sequence) - (length(sequence)%%3))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 2)], 
                      subSample[2:(length(subSample) - 1)], subSample[3:(length(subSample))]), 
                    ncol = 3)
  temp <- as.data.frame(seqSet1)
  Nijk <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  # seqSet2 <- seqSet1[, -3]
  # temp2 <- as.data.frame(seqSet2)
  # Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)]
  # 
  # Nij %<>% as.data.frame()
  Nijk %<>% as.data.frame()
  colnames(Nijk) <- c("start_motif", "middle_motif","end_motif", "counts")
  return(Nijk)
}

get_second_order_motif_transition_probability <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  sequence <- mouse_data %>%
    .$Motif %>%
    as.character()
  
  # transMatrix <- markovchainFit(data = sequence)$estimate@transitionMatrix
  # subSample <- sequence[1:(length(sequence) - (length(sequence)%%4))]
  # seqSet1 <- matrix(c(subSample[1:(length(subSample) - 3)], 
  #                     subSample[2:(length(subSample) - 2)], 
  #                     subSample[3:(length(subSample) - 1)], subSample[4:(length(subSample))]), 
  #                   ncol = 4)
  # temp <- as.data.frame(seqSet1)
  # Nijkl <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  
  subSample <- sequence[1:(length(sequence) - (length(sequence)%%3))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 2)], 
                      subSample[2:(length(subSample) - 1)], subSample[3:(length(subSample))]), 
                    ncol = 3)
  temp <- as.data.frame(seqSet1)
  Nijk <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  seqSet2 <- seqSet1[, -3]
  temp2 <- as.data.frame(seqSet2)
  Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)]

  Nij %<>% as.data.frame()
  Nijk %<>% as.data.frame()
  colnames(Nij) <- c("start_motif","end_motif", "counts")
  colnames(Nijk) <- c("first_motif", "second_motif","end_motif", "counts")
  total_counts <- vector(mode = "numeric")
  for(i in 1:nrow(Nijk)) {
    total_counts[i] <- Nij %>%
      filter(start_motif == as.character(Nijk$first_motif[i]) & end_motif == as.character(Nijk$second_motif[i])) %>%
      .$counts
  }
  Nijk %<>% mutate(!!paste0("probability_", mid):=counts/total_counts, transition=paste(first_motif, second_motif, end_motif, sep = "_")) %>% dplyr::select(-c(first_motif,second_motif,end_motif,counts))
  return(Nijk)
}

get_third_order_motif_transition_counts <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  sequence <- mouse_data %>%
    .$Motif %>%
    as.character()
  

  subSample <- sequence[1:(length(sequence) - (length(sequence)%%4))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 3)],
                      subSample[2:(length(subSample) - 2)], 
                      subSample[3:(length(subSample) - 1)], 
                      subSample[4:(length(subSample))]), 
                    ncol = 4)
  temp <- as.data.frame(seqSet1)
  Nijkl <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  # seqSet2 <- seqSet1[, -3]
  # temp2 <- as.data.frame(seqSet2)
  # Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)]
  # 
  # Nij %<>% as.data.frame()
  Nijkl %<>% as.data.frame()
  colnames(Nijkl) <- c("first_motif", "second_motif", "third_motif", "end_motif", "counts")
  return(Nijkl)
}

get_third_order_motif_transition_probability <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  sequence <- mouse_data %>%
    .$Motif %>%
    as.character()
  

  subSample <- sequence[1:(length(sequence) - (length(sequence)%%4))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 3)],
                      subSample[2:(length(subSample) - 2)], 
                      subSample[3:(length(subSample) - 1)], 
                      subSample[4:(length(subSample))]), 
                    ncol = 4)
  temp <- as.data.frame(seqSet1)
  Nijkl <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  seqSet2 <- seqSet1[, -4]
  temp2 <- as.data.frame(seqSet2)
  Nijk <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)]
  
  Nijk %<>% as.data.frame()
  Nijkl %<>% as.data.frame()
  colnames(Nijk) <- c("first_motif", "second_motif","end_motif", "counts")
  colnames(Nijkl) <- c("first_motif", "second_motif", "third_motif", "end_motif", "counts")
  total_counts <- vector(mode = "numeric")
  for(i in 1:nrow(Nijkl)) {
    total_counts[i] <- Nijk %>%
      filter(first_motif == as.character(Nijkl$first_motif[i]) & second_motif == as.character(Nijkl$second_motif[i]) & end_motif == as.character(Nijkl$third_motif[i])) %>%
      .$counts
  }
  Nijkl %<>% mutate(!!paste0("probability_", mid):=counts/total_counts, transition=paste(first_motif, second_motif, third_motif, end_motif, sep = "_")) %>% dplyr::select(-c(first_motif,second_motif,third_motif, end_motif,counts))
  return(Nijkl)
}

get_motif_probability <- function(mid, data, pheno_data) {
  mouse_data <- get_mouse_data(mid,data,time_step)
  
  probability <- mouse_data %>%
    .$Motif %>%
    table() %>%
    divide_by(nrow(mouse_data))
  return(probability)
}



mouse_ids <- pheno_data %>% filter(Genotype == 1) %>% .$Mouse
log_like_3 <- vector(mode = "numeric")
log_like_2 <- vector(mode = "numeric")
log_like_1 <- vector(mode = "numeric")
log_like_0 <- vector(mode = "numeric")
tempCount <- 1

for(mid in mouse_ids) {
  print(tempCount)
  
  ##third_order
  ##second_order
  third_order_trans_prob_results <- lapply(setdiff(mouse_ids,mid), get_third_order_motif_transition_probability, data, pheno_data)
  all_third_order_trans_prob_results <- merge(third_order_trans_prob_results[[1]], third_order_trans_prob_results[[2]], all=TRUE, by="transition")
  for(i in 3:length(third_order_trans_prob_results)) {
    all_third_order_trans_prob_results %<>% merge(.,third_order_trans_prob_results[[i]], all=TRUE, by="transition")
  }
  all_third_order_trans_prob_results[is.na(all_third_order_trans_prob_results)] <- 0
  
  probability <- all_third_order_trans_prob_results %>%
    dplyr::select(starts_with("probability")) %>%
    rowMeans()
  third_order_transition_probability <- data.frame(transition=all_third_order_trans_prob_results$transition, probability)
  Nijkl <- get_third_order_motif_transition_counts(mid, data, pheno_data) %>%
    mutate(transition=paste(first_motif, second_motif, third_motif, end_motif, sep = "_"))
  
  ##NOTE: this step (in the merge function) involves ignoring transitions observed only in the mid mouse and not in any of the other mice
  ##MAY NEED TO CHECK THIS LATER
  temp_3_order_like <- Nijkl %>%
    merge(., third_order_transition_probability) %>%
    mutate(temp_like = counts*log(probability))
  log_like_3[tempCount] <- temp_3_order_like %>%
    .$temp_like %>%
    sum()
  
  ##second_order
  second_order_trans_prob_results <- lapply(setdiff(mouse_ids,mid), get_second_order_motif_transition_probability, data, pheno_data)
  all_second_order_trans_prob_results <- merge(second_order_trans_prob_results[[1]], second_order_trans_prob_results[[2]], all=TRUE, by="transition")
  for(i in 3:length(second_order_trans_prob_results)) {
    all_second_order_trans_prob_results %<>% merge(.,second_order_trans_prob_results[[i]], all=TRUE, by="transition")
  }
  all_second_order_trans_prob_results[is.na(all_second_order_trans_prob_results)] <- 0
  
  probability <- all_second_order_trans_prob_results %>%
    dplyr::select(starts_with("probability")) %>%
    rowMeans()
  second_order_transition_probability <- data.frame(transition=all_second_order_trans_prob_results$transition, probability)
  Nijk <- get_second_order_motif_transition_counts(mid, data, pheno_data) %>%
    mutate(transition=paste(start_motif, middle_motif, end_motif, sep = "_"))
  
  ##NOTE: this step (in the merge function) involves ignoring transitions observed only in the mid mouse and not in any of the other mice
  ##MAY NEED TO CHECK THIS LATER
  temp_2_order_like <- Nijk %>%
    merge(., second_order_transition_probability) %>%
    mutate(temp_like = counts*log(probability))
  log_like_2[tempCount] <- temp_2_order_like %>%
    .$temp_like %>%
    sum()
  
  ## first_order
  trans_prob_results <- lapply(setdiff(mouse_ids,mid), get_transition_matrix , data, pheno_data)
  
  mean_trans_prob_results <- trans_prob_results[[1]]
  for(i in 2:length(trans_prob_results)) {
    mean_trans_prob_results <- mean_trans_prob_results + trans_prob_results[[i]]
  }
  mean_trans_prob_results <- mean_trans_prob_results/length(trans_prob_results)
  
  transition_counts <- get_motif_transition_counts(mid, data, pheno_data)
  
  log_like_1[tempCount] <- 0
  for(i in 1:nrow(transition_counts)) {
    temp_prob <- mean_trans_prob_results[which(row.names(mean_trans_prob_results)==as.character(transition_counts$start_motif[i])), which(row.names(mean_trans_prob_results)==as.character(transition_counts$end_motif[i]))]
    if(temp_prob > 1e-5) {
    log_like_1[tempCount] <- log_like_1[tempCount] + transition_counts$counts[i]*log(temp_prob)
    }
  }
  
  ##zero_order
  prob_results <- lapply(setdiff(mouse_ids,mid), get_motif_probability , data, pheno_data)
  mean_prob_results <- prob_results[[1]]
  for(i in 2:length(prob_results)) {
    mean_prob_results <- mean_prob_results + prob_results[[i]]
  }
  mean_prob_results <- mean_prob_results/length(prob_results)
  
  motif_counts <- get_motif_counts(mid, data, pheno_data)
  log_like_0[tempCount] <- sum(motif_counts*log(mean_prob_results))
  tempCount <- tempCount + 1
}

plot_data <- data.frame(log_likelihood=append(append(log_like_0, append(log_like_1, log_like_2)), log_like_3), markov_order=rep(c("0","1", "2", "3"), each=length(mouse_ids)), mouse_id=as.character(mouse_ids))

sub_pheno_data <- pheno_data %>% 
  filter(Genotype == 1) %>% 
  arrange(Mouse)

plot_data %<>% cbind(., Sex=as.factor(sub_pheno_data$Sex), Genotype=as.factor(sub_pheno_data$Genotype))

plot_data %>%
  write.csv(.,paste0("disease_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_", resolutionLabel,"second.csv"), row.names = FALSE)
resolutionLabel <- sub("[.]","_",as.character(resolution))

p <- (ggplot(plot_data, aes(x=markov_order, y=log_likelihood, group=mouse_id, color=Sex)) + 
  geom_point() +
  geom_line() +
  ylab("Leave-One-Out Log Likelihood") +
  xlab("markov order")) 
p %>% ggsave(file = paste0("disease_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_", resolutionLabel,"second.pdf"),
         plot = .,
         width = 7, 
         height = 7)

mouse_data %>%
  dplyr::slice(15000:30000) %>% 
  .$Motif %>%
  as.character() %>%
  verifyMarkovProperty(., verbose = TRUE)

mouse_data %>%
  .$Motif %>%
  as.character() %>%
  assessOrder(.)

mouse_data %>%
  .$Motif %>%
  as.character() %>%
  assessStationarity(.,1)

View(mcFitMLE$estimate)

motifs <- data %>%
  .$SubCommunity %>%
  unique() %>%
  sort()

mouse_data %>%
  .$Community %>%
  as.character() %>%
  verifyMarkovProperty(.)


test_sequence <- mouse_data %>%
  .$Motif %>%
  as.character()
sequence <- test_sequence


.findNijPjk<-function(Nijk = Nijk, Nij = Nij, trans, row = 1){
  i <- Nijk[row,1]
  j <- Nijk[row,2]
  k <- Nijk[row,3]
  
  fromCh <- as.character(j)
  toCh <- as.character(k)
  Pjk <- trans[fromCh,toCh]
  
  m1 <- which(Nij[, 1] == i)
  m2 <- which(Nij[, 2] == j)
  m <- c(m1, m2)
  return(Nij[m[anyDuplicated(m)], 3] * Pjk)
}

verifyMarkovProperty
function (sequence, verbose = TRUE) 
{
  transMatrix <- markovchainFit(data = sequence)$estimate@transitionMatrix
  subSample <- sequence[1:(length(sequence) - (length(sequence)%%3))]
  seqSet1 <- matrix(c(subSample[1:(length(subSample) - 2)], 
                      subSample[2:(length(subSample) - 1)], subSample[3:(length(subSample))]), 
                    ncol = 3)
  temp <- as.data.frame(seqSet1)
  Nijk <- aggregate(temp, by = temp, length)[1:(ncol(temp) + 
                                                  1)]
  seqSet2 <- seqSet1[, -3]
  temp2 <- as.data.frame(seqSet2)
  Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 
                                                   1)]
  test <- c(length = dim(Nijk)[1])
  invisible(lapply(seq_len(dim(Nijk)[1]), function(i) {
    foundNijPjk <- .findNijPjk(Nijk = Nijk, Nij = Nij, trans = transMatrix, 
                               row = i)
    
    test[i] <<- foundNijPjk
  }))
  
  Nijk <- cbind(Nijk, round(test))
  
  colnames(Nijk) <- c("start_state", "middle_state", "end_state", "observed", "expected")
  Nijk %<>% as.data.frame() %>% mutate(difference=observed-expected)
  
  statistic <- sum(test)
  doubles = numeric(length(sequence) - 1)
  for (i in 1:(length(doubles))) {
    doubles[i] = paste(sequence[i], sequence[i + 1], sep = "", 
                       collapse = NULL)
  }
  triples = numeric(length(sequence) - 2)
  for (i in 1:(length(triples))) {
    triples[i] = paste(sequence[i], sequence[i + 1], sequence[i + 
                                                                2], sep = "", collapse = NULL)
  }
  dof = length(unique(triples)) - length(unique(doubles)) + 
    length(unique(sequence)) - 1
  pvalue <- 1 - pchisq(q = statistic, df = dof)
  out <- list(statistic = statistic, dof = dof, p.value = pvalue)
  if (verbose == TRUE) {
    cat("Testing markovianity property on given data sequence\n")
    cat("Chi - square statistic is:", statistic, "\n")
    cat("Degrees of freedom are:", dof, "\n")
    cat("And corresponding p-value is:", pvalue, "\n")
  }
  invisible(out)
}
