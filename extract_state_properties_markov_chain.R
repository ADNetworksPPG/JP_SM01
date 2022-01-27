rm(list = ls())
setwd("~/Dropbox (Gladstone)/GB-SM-1238/CVS time motif info")
require(dplyr)
require(magrittr)
require(markovchain)
require(ggplot2)
data <- read.csv("210416 motif_time 25Hz v5.csv", header = TRUE)

time_step <- 16
print(time_step)
resolution <- 0.04*time_step

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
motif_communities <- paste0("MotifCommunities.csv") %>%
  read.csv(., header = TRUE)
data %<>% merge(., motif_communities)

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
    dplyr::summarise(mode = mode(SubCommunity))
  
  mouse_data %<>% slice(rows_to_use) %>%
    mutate(Motif=sub_community$mode)
  
  return(mouse_data)
}

description_motifs <- read.csv("Description_Motifs.csv", header = TRUE)
mid <- 13
state_transition <- c(0,11)
mouse_data <- get_mouse_data(mid, data, time_step)

##simulate states
## male-disease mouse
mid <- 105
mouse_data <- get_mouse_data(mid, data, time_step)
mouse_data %<>% merge(., description_motifs, all.x=TRUE) %>% 
  arrange(Seconds) %>%
  dplyr::select(Seconds, Description, Motif_order) %>%
  mutate(Minutes = Seconds/60, Description=factor(Description, 
                               levels = c("Groom", "Sit-down", "Sit-up", "Rear", "Slow-exlore", "Fast-explore", "Walk", "Walk-raise", "Run")))

#mouse_data %<>% mutate(Motif_order = as.factor(Motif_order))
(ggplot(mouse_data, aes(x=Minutes, y=Motif_order)) + 
  geom_step() +
  xlim(-2,5) +
  ylim(-1,10) +
  ylab("") +
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  theme(axis.line.y = element_blank()) +
  geom_text(x=-1, y=0, label=c("Groom")) +
  geom_text(x=-1, y=1, label=c("Sit-down")) +
  geom_text(x=-1, y=2, label=c("Sit-up")) +
  geom_text(x=-1, y=3, label=c("Rear")) +
  geom_text(x=-1, y=4, label=c("Slow-explore")) +
  geom_text(x=-1, y=5, label=c("Fast-explore")) +
  geom_text(x=-1, y=6, label=c("Walk")) +
  geom_text(x=-1, y=7, label=c("Walk-raise")) +
  geom_text(x=-1, y=8, label=c("Run")) +
  theme(text = element_text(size = 20))) %>%
  ggsave(plot=.,
         file = "state_transitions_mouse_105_male_control.pdf",
         width=10,
         height=7)

#set.seed(234)
rmc <- mouse_data %>% 
  nrow() %>%
  rmarkovchain(., mouse_mc[[20]])

random_mouse_data <- data.frame(Motif=rmc)
random_mouse_data %<>% mutate(Minutes=mouse_data$Minutes)
random_mouse_data %<>% merge(., description_motifs, all.x=TRUE)

(ggplot(random_mouse_data, aes(x=Minutes, y=Motif_order)) + 
    geom_step() +
    xlim(-2,5) +
    ylim(-1,10) +
    ylab("") +
    theme_classic() +
    theme(axis.text.y = element_blank()) +
    theme(axis.line.y = element_blank()) +
    geom_text(x=-1, y=0, label=c("Groom")) +
    geom_text(x=-1, y=1, label=c("Sit-down")) +
    geom_text(x=-1, y=2, label=c("Sit-up")) +
    geom_text(x=-1, y=3, label=c("Rear")) +
    geom_text(x=-1, y=4, label=c("Slow-explore")) +
    geom_text(x=-1, y=5, label=c("Fast-explore")) +
    geom_text(x=-1, y=6, label=c("Walk")) +
    geom_text(x=-1, y=7, label=c("Walk-raise")) +
    geom_text(x=-1, y=8, label=c("Run")) +
    theme(text = element_text(size = 20))) %>%
  ggsave(plot=.,
         file = "state_transitions_random_mouse_105_male_control.pdf",
         width=10,
         height=7)

get_transition_matrix <- function(mid, data, time_step) {
  mouse_data <- get_mouse_data(mid, data, time_step)
  
  mcFitMLE <- mouse_data %>%
    dplyr::slice(1:nrow(.)) %>% 
    .$Motif %>%
    as.character() %>%
    markovchainFit(data = .)
  
  # mouse_mc <- new("markovchain", states=as.character(motifs), transitionMatrix=mcFitMLE$estimate@transitionMatrix)
  transitionMatrix=mcFitMLE$estimate@transitionMatrix
  return(transitionMatrix)
}

get_mouse_markovchain <- function(mid, data, time_step) {
  mouse_data <- get_mouse_data(mid, data, time_step)
  
  mcFitMLE <- mouse_data %>%
    dplyr::slice(1:nrow(.)) %>% 
    .$Motif %>%
    as.character() %>%
    markovchainFit(data = .)
  
  mouse_mc <- new("markovchain", transitionMatrix=mcFitMLE$estimate@transitionMatrix)
  return(mouse_mc)
}

mouse_mc <- lapply(mouse_ids, get_mouse_markovchain, data,  time_step)

mean_recurrence_times <- lapply(mouse_mc, meanRecurrenceTime) %>% 
  do.call("rbind", .) %>%
  as.data.frame() 

states <- colnames(mean_recurrence_times)
mean_recurrence_times %<>%
  mutate(mouse_id = mouse_ids)

pheno_data %<>% rename(mouse_id = Mouse) 
pheno_data %<>% mutate(Genotype=as.factor(Genotype), Sex=as.factor(Sex))
mean_recurrence_times %<>% merge(., pheno_data)

mc_feature_data <- mean_recurrence_times
state <- "A1"

test_genotype_association <-  function(state, mc_feature_data) {
  #print(state)
  fit_data <- mc_feature_data %>% 
    dplyr::select(all_of(state), Genotype, Sex) %>% 
    rename(Y=all_of(state))
  
  lmFit <- lm(Y ~ Genotype + Sex + Genotype:Sex, data=fit_data)
  slmFit <-  summary(lmFit)
  res <- slmFit$coefficients[-1,] %>% as.data.frame()
  row.names(res) <-  paste(state, row.names(res), sep="_")
  return(res)
}

recur_res <- lapply(states, test_genotype_association, mean_recurrence_times) %>%
  do.call("rbind", .) %>%
  rename(pvalue="Pr(>|t|)") %>%
  mutate(padjust=p.adjust(pvalue, "BH")) %>%
  write.csv(., "difference_in_recurrence_times.csv")

plot_genotype_association <-  function(state, mc_feature_data, resolution) {
  #print(state)
  fit_data <- mc_feature_data %>% 
    dplyr::select(all_of(state), Genotype, Sex) %>% 
    rename(Y=all_of(state)) %>%
    mutate(Y=resolution*Y)
  
  dGenotype <- rep("control", nrow(fit_data))
  dGenotype[fit_data$Genotype == "1"] <- "disease"

  dSex <- rep("female", nrow(fit_data))
  dSex[fit_data$Sex == "1"] <- "male"
  
  fit_data %<>% mutate(dGenotype=dGenotype, Sex=dSex)
  

  print(ggplot(fit_data, aes(x=dGenotype, y=Y, color=Sex)) + 
          geom_boxplot() + 
          ggtitle(paste0("Return to state_", state)) +
          ylab("Time (in seconds)") +
          xlab("Genotype") +
          theme(text = element_text(size = 20)))
}

(ggplot(fit_data, aes(x=dGenotype, y=Y, color=Sex)) + 
        geom_boxplot() + 
        ggtitle(paste0("Return to Sit-down")) +
        ylab("Time (in seconds)") +
        xlab("Genotype") +
        theme(text = element_text(size = 20))) %>%
  ggsave(plot=.,
         file="average_time_to_sitting_down.pdf",
         width = 10,
         height=7)
plot_genotype_association("G", mean_recurrence_times)
mean_first_passage_times <- lapply(mouse_mc, meanFirstPassageTime)

test_mean_first_passage_genotype_association <-  function(state1, state2, mean_first_passage_times, pheno_data) {
  #print(state)
  times <- vector(mode = "numeric")
  for(i in 1:length(mean_first_passage_times)) {
    temp_data <- mean_first_passage_times[[i]] %>%
      as.data.frame() %>%
      filter(row.names(.)==state1) %>%
      dplyr::select(matches(state2))
    times[i] <- temp_data[1,1]
  }
  fit_data <- data.frame(times, pheno_data)
  lmFit <- lm(times ~ Genotype + Sex + Genotype:Sex, data=fit_data)
  slmFit <-  summary(lmFit)
  res <- slmFit$coefficients[-1,] %>% as.data.frame()
  row.names(res) <-  paste(state1, state2,row.names(res), sep="_")
  return(res)
}

plot_mean_first_passage_genotype_association <-  function(state1, state2, mean_first_passage_times, pheno_data, resolution) {
  #print(state)
  times <- vector(mode = "numeric")
  for(i in 1:length(mean_first_passage_times)) {
    temp_data <- mean_first_passage_times[[i]] %>%
      as.data.frame() %>%
      filter(row.names(.)==state1) %>%
      dplyr::select(matches(state2))
    times[i] <- temp_data[1,1]
  }
  fit_data <- data.frame(times, pheno_data)
  
  dGenotype <- rep("control", nrow(fit_data))
  dGenotype[fit_data$Genotype == "1"] <- "disease"
  
  dSex <- rep("female", nrow(fit_data))
  dSex[fit_data$Sex == "1"] <- "male"
  
  fit_data %<>% mutate(dGenotype=dGenotype, Sex=dSex, times = resolution*times)
  
  print(ggplot(fit_data, aes(x=dGenotype, y=times, color=Sex)) + 
          geom_boxplot() + 
          ggtitle(paste0("state_", state1, "_to_state_", state2)) +
          xlab("Genotype") +
          ylab("Time (in seconds)") +
          theme(text = element_text(size = 20)))
  
}

(ggplot(fit_data, aes(x=dGenotype, y=times, color=Sex)) + 
        geom_boxplot() + 
        ggtitle(paste0("Average time to change from Groom to Rear")) +
        xlab("Genotype") +
        ylab("Time (in seconds)") +
        theme(text = element_text(size = 20))) %>%
  ggsave(plot=.,
         file="average_time_to_change_from_Groom_to_Rear.pdf",
         width = 10,
         height=7)

first_passage_time_fits <- list()
temp_count <- 1
for(state1 in states) {
  for(state2 in states) {
    if(state1 != state2) {
      first_passage_time_fits[[temp_count]] <- test_mean_first_passage_genotype_association(state1, state2, mean_first_passage_times, pheno_data)
      temp_count <- temp_count + 1
    }
  }
}

first_passage_time_fits_res <- do.call("rbind", first_passage_time_fits) %>%
  rename(pvalue="Pr(>|t|)") %>%
  mutate(padjust=p.adjust(pvalue, "BH")) %>%
  write.csv(., "differences_in_first_passage_times.csv")
