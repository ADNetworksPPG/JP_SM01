rm(list = ls())
setwd("~/Gladstone Dropbox/Reuben Thomas/JP_SM02/data/")
require(dplyr)
require(magrittr)
data <- read.csv("250325 AppSAA Cohort 1 Bioinformatics Network Analysis Data Truncated.csv", header = TRUE)
#motif_communities <- read.csv("MotifCommunities.csv", header = TRUE)
motif_communities <- table(data$Community_Label) %>% names(.)

data %<>% merge(., motif_communities, all.x=TRUE, by="Motif")

mouse_ids <- data %>%
  .$Mouse %>%
  unique() %>%
  sort()

motifs <- data %>%
  .$Community_Label %>%
  unique() %>%
  sort()

pheno_data <- data %>%
  select(MouseID,GroupCode, Sex) %>%
  unique() %>%
  arrange(MouseID)


mid <- 2
state_transition <- c("A1","C")
get_transition_info <- function(mid, state_transition, data) {
  mouse_data <- data %>%
    filter(MouseID == mid) %>%
    arrange(Seconds)
  
  start_indices <- which(mouse_data$Community_Label == state_transition[1])
  start_indices <- start_indices[start_indices < nrow(mouse_data)]
  no_state_transitions <- sum(mouse_data$Community_Label[start_indices+1] == state_transition[2])
  total_transitions <- length(start_indices)
  return(c(no_state_transitions, total_transitions))
}

motifs %<>% as.character()
#data %<>% mutate(SubCommunity = as.character(SubCommunity))

state_transitions <- list()
temp_count <- 1
for(i in 1:length(motifs)) {
  for(j in 1:length(motifs)) {
    state_transitions[[temp_count]] <- c(motifs[i],motifs[j])
    temp_count <- temp_count + 1
  }
}

estimate_difference_transition_probability <- function(state_transition, mouse_ids, data, pheno_data) {
  require(lme4)
  print(state_transition)
  freq_data <- sapply(mouse_ids, get_transition_info, state_transition, data) %>%
    t()
  
  fit_data <- data.frame(mouse_ids, freq_data)
  colnames(fit_data) <- c("MouseID","Count", "Total")
  fit_data %<>% merge(.,pheno_data)
  fit_data %<>% mutate(GroupCode=as.factor(GroupCode)) %>% arrange(GroupCode)
  fit_data %<>% mutate(GroupCode = relevel(GroupCode, ref = "WT"))
  glmer_fit <- glmer(cbind(Count, (Total - Count)) ~ (1|MouseID) + GroupCode, family = "binomial", data = fit_data)
  sglmer_fit <-  summary(glmer_fit)
  res <- sglmer_fit$coefficients[-1,]
  
  
  summary_counts <- fit_data %>%
    group_by(GroupCode) %>%
    summarise(mean_counts=mean(Count), mean_total_counts=mean(Total))
  # ggplot(fit_data, aes(x=Genotype, y=Count)) + geom_boxplot() + ggtitle(paste0("Transition from motif ", state_transition[1], " to motif ", state_transition[2]))
  # pois_fit <- glm(Count ~ Genotype, data=fit_data, family = "poisson")
  # spoisfit <- summary(pois_fit)
  # require(lme4)
  # lmer_fit_data <- NULL
  # for(i in 1:nrow(fit_data)) {
  #   Y <- rep(0, fit_data$Total[i])
  #   if(fit_data$Count[i] > 0) {
  #     Y[1:fit_data$Count[i]] <- 1
  #   }
  #   MouseID <- rep(as.character(fit_data$Mouse[i]), fit_data$Total[i])
  #   GroupCode <- rep(as.character(fit_data$GroupCode[i]), fit_data$Total[i])
  #   Sex <- rep(as.character(fit_data$Sex[i]), fit_data$Total[i])
  #   temp_data <- data.frame(Y, MouseID, GroupCode, Sex)
  #   lmer_fit_data %<>% rbind(.,temp_data)
  # }
  # res <- rep(NA, 4)
  # res[4] <- 1
  # if(sum(lmer_fit_data$Y) > 0) {
  #   glmer_fit <- glmer(Y ~ (1|MouseID) + GroupCode, family = "binomial", data = lmer_fit_data)
  #   sglmer_fit <-  summary(glmer_fit)
  #   res <- sglmer_fit$coefficients[-1,]
  # }
  res %<>% append(., summary_counts[1,2])
  res %<>% append(., summary_counts[2,2])
  res %<>% append(., summary_counts[1,3])
  res %<>% append(., summary_counts[2,3])
  
  res %<>% as.numeric()
  names(res)[1:8] <- c("estimate","std.error","z.value","p.value", paste0("mean_counts_genotype_", 0:1), paste0("mean_total_counts_genotype_", 0:1))
  
  return(res)
}

trans_prob_diff <- sapply(state_transitions, estimate_difference_transition_probability, mouse_ids, data, pheno_data)
trans_prob_diff %<>% t()

trans_prob_diff %<>% as.data.frame()
colnames(trans_prob_diff)[4] <- "p.value"
mat_state_transitions <- do.call(rbind, state_transitions) %>% as.data.frame()
colnames(mat_state_transitions) <- c("start_community", "end_community")

results <- cbind(mat_state_transitions, trans_prob_diff)

results %<>% mutate(p.adjust=p.adjust(p.value, method = "BH"))
write.csv(results, "glmer_diff_trans_probabilities_community_250325_AppSAA_Cohort_1.csv", row.names = FALSE)
