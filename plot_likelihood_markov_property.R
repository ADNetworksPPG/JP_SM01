rm(list = ls())
setwd("~/Dropbox (Gladstone)/GB-SM-1238/CVS time motif info")
require(dplyr)
require(tidyr)
require(ggplot2)

##get gender information
gender_data <- read.csv("210416 motif_time 25Hz v5.csv", header = TRUE) %>%
  dplyr::select(Mouse, Genotype, Sex) %>% 
  unique() %>%
  mutate(mouse_id = paste0("mid_", Mouse))

time_step <- c(1, 6, 9, 12, 16, 25, 50, 125, 250) 
#time_step <- 16
for(i in 1:length(time_step)) {
  resolution <- 0.04*time_step[i]
  print(resolution)
  resolutionLabel <- sub("[.]","_",as.character(resolution))
  
  control_plot_data <- read.csv(paste0("updated_analyses_control_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel,"_seconds.csv"), header = TRUE)
  n_control <- ncol(control_plot_data)
  control_plot_data %<>% mutate(markov_order=as.character(0:(nrow(.)-1)),condition=rep("control", nrow(.)))
  control_plot_data %<>%
    gather(., mouse_id, log_like, 1:all_of(n_control))
  
  
  disease_plot_data <- read.csv(paste0("updated_analyses_disease_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel,"_seconds.csv"), header = TRUE)
  n_disease <- ncol(disease_plot_data)
  disease_plot_data %<>% mutate(markov_order=as.character(0:(nrow(.)-1)), condition=rep("disease", nrow(.)))
  disease_plot_data %<>%
    gather(., mouse_id, log_like, 1:all_of(n_disease))
  
  plot_data <- rbind(control_plot_data, disease_plot_data)
  
  plot_data %<>% merge(., gender_data, all.x=TRUE)
  plot_data %<>% mutate(Sex = as.factor(Sex))
  p <- (ggplot(plot_data, aes(x=markov_order, y=log_like, group=mouse_id, color=Sex)) + 
          facet_grid(~condition) +
          geom_point() +
         #ylim(-5000,-2500) +
          geom_line() +
          ylab("Leave-One-Out Log Likelihood") +
          xlab("markov order")) 
  p %>% ggsave(file = paste0("updated_analyses_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_", resolutionLabel,"seconds.pdf"),
               plot = .,
               width = 10, 
               height = 7)
  
  sub_plot_data <- plot_data %>%
    filter(markov_order == "0" | markov_order == "1" | markov_order == "2")
  p1 <- (ggplot(sub_plot_data, aes(x=markov_order, y=log_like, group=mouse_id, color=Sex)) + 
      facet_grid(~condition) +
      geom_point() +
      #ylim(-5000,-2500) +
      geom_line() +
      ylab("Leave-One-Out Log Likelihood") +
      xlab("markov order"))
  
  p1 %>% ggsave(file = paste0("updated_analyses_mouse_behavior_markov_order_0_to_2_likelihood_sub_communities_time_resolution_", resolutionLabel,"seconds.pdf"),
               plot = .,
               width = 10, 
               height = 7)
  
  sub_data <- plot_data %>%
    filter(markov_order == "1")
  
  print(t.test(log_like ~ condition, data = sub_data, var.equal=FALSE))
  
  p3 <- ggplot(sub_data, aes(x=condition, y=log_like, color=Sex)) +
    #geom_point() +
    #geom_jitter() +
    geom_boxplot()
  p3 %>% ggsave(file = paste0("updated_analyses_mouse_behavior_markov_order_1_likelihood_sub_communities_time_resolution_", resolutionLabel,"seconds.pdf"),
                plot = .,
                width = 7, 
                height = 7)
  summ_data <- sub_plot_data %>%
    group_by(condition, markov_order) %>%
    summarise(mean_log_like = mean(log_like), se_log_like = sd(log_like)/sqrt(nrow(.)))
  p2 <- (ggplot(summ_data, aes(x=markov_order, y=mean_log_like, color=condition, group=condition)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean_log_like - 1.96*se_log_like), ymax = (mean_log_like + 1.96*se_log_like)), width = 0.2))
  
  p2 %>% ggsave(file = paste0("updated_analyses_mouse_behavior_markov_order_0_to_2_mean_likelihood_sub_communities_time_resolution_", resolutionLabel,"seconds.pdf"),
                plot = .,
                width = 10, 
                height = 7)
}

##compute pseudo bayes factors for 1 vs 0, 2 vs 1 and 2 vs 0

control_pseudo_bayes_factors <- matrix(NA, length(time_step), 5)
control_pseudo_bayes_factors_sex0 <- matrix(NA, length(time_step), 5)
control_pseudo_bayes_factors_sex1 <- matrix(NA, length(time_step), 5)

disease_pseudo_bayes_factors <- matrix(NA, length(time_step), 5)
disease_pseudo_bayes_factors_sex0 <- matrix(NA, length(time_step), 5)
disease_pseudo_bayes_factors_sex1 <- matrix(NA, length(time_step), 5)

for(i in 1:length(time_step)) {
  resolution <- 0.04*time_step[i]
  resolutionLabel[i] <- sub("[.]","_",as.character(resolution))
  
  control_plot_data <- read.csv(paste0("control_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel[i],"second.csv"), header = TRUE)
  control_plot_data %<>% mutate(condition=rep("control", nrow(.)))
  
  disease_plot_data <- read.csv(paste0("disease_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel[i],"second.csv"), header = TRUE)
  disease_plot_data %<>% mutate(condition=rep("disease", nrow(.)))
  
  #model 1 vs model 0 - control
  control_pseudo_bayes_factors[i,1] <- ((control_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "0") %>% .$log_likelihood %>% sum())) 
  #model 1 vs model 0 - control & Sex = 0
  control_pseudo_bayes_factors_sex0[i,1] <- ((control_plot_data %>% filter(markov_order == "1" & Sex == 0) %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "0" & Sex == 0) %>% .$log_likelihood %>% sum())) 
  #model 1 vs model 0 - control & Sex = 1
  control_pseudo_bayes_factors_sex1[i,1] <- ((control_plot_data %>% filter(markov_order == "1" & Sex == 1) %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "0" & Sex == 1) %>% .$log_likelihood %>% sum())) 
  
  #model 1 vs model 0 - disease
  disease_pseudo_bayes_factors[i,1] <- ((disease_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "0") %>% .$log_likelihood %>% sum())) 

  #model 1 vs model 0 - disease & Sex = 0
  disease_pseudo_bayes_factors_sex0[i,1] <- ((disease_plot_data %>% filter(markov_order == "1" & Sex == 0) %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "0" & Sex == 0) %>% .$log_likelihood %>% sum())) 
  #model 1 vs model 0 - disease & Sex = 1
  disease_pseudo_bayes_factors_sex1[i,1] <- ((disease_plot_data %>% filter(markov_order == "1" & Sex == 1) %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "0" & Sex == 1) %>% .$log_likelihood %>% sum())) 
  
  #model 2 vs model 2 - control
  control_pseudo_bayes_factors[i,2] <- ((control_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum())) 
  #model 2 vs model 1 - disease
  disease_pseudo_bayes_factors[i,2] <- ((disease_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum())) 
  
  #model 2 vs model 0 - control
  control_pseudo_bayes_factors[i,3] <- ((control_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "0") %>% .$log_likelihood %>% sum())) 
  #model 2 vs model 0 - disease
  disease_pseudo_bayes_factors[i,3] <- ((disease_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "0") %>% .$log_likelihood %>% sum())) 
  
  #model 3 vs model 2 - control
  control_pseudo_bayes_factors[i,4] <- ((control_plot_data %>% filter(markov_order == "3") %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum())) 
  #model 3 vs model 2 - disease
  disease_pseudo_bayes_factors[i,4] <- ((disease_plot_data %>% filter(markov_order == "3") %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "2") %>% .$log_likelihood %>% sum())) 

  #model 3 vs model 1 - control
  control_pseudo_bayes_factors[i,5] <- ((control_plot_data %>% filter(markov_order == "3") %>% .$log_likelihood %>% sum()) - 
                                          (control_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum())) 
  #model 3 vs model 1 - disease
  disease_pseudo_bayes_factors[i,5] <- ((disease_plot_data %>% filter(markov_order == "3") %>% .$log_likelihood %>% sum()) - 
                                          (disease_plot_data %>% filter(markov_order == "1") %>% .$log_likelihood %>% sum())) 
  
}

control_pseudo_bayes_factors %<>% as.data.frame()
row.names(control_pseudo_bayes_factors) <- paste0("resolution_", resolutionLabel)
colnames(control_pseudo_bayes_factors) <- c("order1_vs_order0", "order2_vs_order1", "order2_vs_order0", "order3_vs_order2", "order3_vs_order1")

disease_pseudo_bayes_factors %<>% as.data.frame()
row.names(disease_pseudo_bayes_factors) <- paste0("resolution_", resolutionLabel)
colnames(disease_pseudo_bayes_factors) <- c("order1_vs_order0", "order2_vs_order1", "order2_vs_order0", "order3_vs_order2", "order3_vs_order1")

write.csv(control_pseudo_bayes_factors, "control_pseudo_bayes_factors_resolution_0_64_seconds.csv", row.names = FALSE)
write.csv(disease_pseudo_bayes_factors, "disease_pseudo_bayes_factors_resolution_0_64_seconds.csv", row.names = FALSE)

control_pseudo_bayes_factors <- matrix(NA, length(time_step), 6)
control_pseudo_bayes_factors_sex0 <- matrix(NA, length(time_step), 6)
control_pseudo_bayes_factors_sex1 <- matrix(NA, length(time_step), 6)

disease_pseudo_bayes_factors <- matrix(NA, length(time_step), 6)
disease_pseudo_bayes_factors_sex0 <- matrix(NA, length(time_step), 6)
disease_pseudo_bayes_factors_sex1 <- matrix(NA, length(time_step), 6)

markov_orders <- matrix(NA,6,2)
for(i in 1:length(time_step)) {
  resolution <- 0.04*time_step[i]
  resolutionLabel[i] <- sub("[.]","_",as.character(resolution))
  
  control_plot_data <- read.csv(paste0("control_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel[i],"second.csv"), header = TRUE)
  control_plot_data %<>% mutate(condition=rep("control", nrow(.)))
  
  disease_plot_data <- read.csv(paste0("disease_mouse_behavior_markov_order_likelihood_sub_communities_time_resolution_",resolutionLabel[i],"second.csv"), header = TRUE)
  disease_plot_data %<>% mutate(condition=rep("disease", nrow(.)))
  
  temp_count <- 1
  for(o1 in 0:2) {
    for(o2 in (o1+1):3) {
      markov_orders[temp_count,1] <- o2
      markov_orders[temp_count, 2] <- o1
      control_pseudo_bayes_factors[i,temp_count] <- ((control_plot_data %>% filter(markov_order == o2) %>% .$log_likelihood %>% sum()) - 
                                              (control_plot_data %>% filter(markov_order == o1) %>% .$log_likelihood %>% sum())) 
      
      disease_pseudo_bayes_factors[i,temp_count] <- ((disease_plot_data %>% filter(markov_order == o2) %>% .$log_likelihood %>% sum()) - 
                                              (disease_plot_data %>% filter(markov_order == o1) %>% .$log_likelihood %>% sum())) 
      
      #model 1 vs model 0 - control & Sex = 0
      control_pseudo_bayes_factors_sex0[i,temp_count] <- ((control_plot_data %>% filter(markov_order == o2 & Sex == 0) %>% .$log_likelihood %>% sum()) - 
                                                   (control_plot_data %>% filter(markov_order == o1 & Sex == 0) %>% .$log_likelihood %>% sum())) 
      #model 1 vs model 0 - control & Sex = 1
      control_pseudo_bayes_factors_sex1[i,temp_count] <- ((control_plot_data %>% filter(markov_order == o2 & Sex == 1) %>% .$log_likelihood %>% sum()) - 
                                                   (control_plot_data %>% filter(markov_order == o1 & Sex == 1) %>% .$log_likelihood %>% sum())) 
      
      #model 1 vs model 0 - disease & Sex = 0
      disease_pseudo_bayes_factors_sex0[i,temp_count] <- ((disease_plot_data %>% filter(markov_order == o2 & Sex == 0) %>% .$log_likelihood %>% sum()) - 
                                                   (disease_plot_data %>% filter(markov_order == o1 & Sex == 0) %>% .$log_likelihood %>% sum())) 
      #model 1 vs model 0 - disease & Sex = 1
      disease_pseudo_bayes_factors_sex1[i,temp_count] <- ((disease_plot_data %>% filter(markov_order == o2 & Sex == 1) %>% .$log_likelihood %>% sum()) - 
                                                   (disease_plot_data %>% filter(markov_order == o1 & Sex == 1) %>% .$log_likelihood %>% sum())) 
      temp_count <- temp_count + 1
    }
  }
 
}

control_pseudo_bayes_factors %<>% as.data.frame()
control_pseudo_bayes_factors_sex0 %<>% as.data.frame()
control_pseudo_bayes_factors_sex1 %<>% as.data.frame()
colnames(control_pseudo_bayes_factors) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")
colnames(control_pseudo_bayes_factors_sex0) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")
colnames(control_pseudo_bayes_factors_sex1) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")

disease_pseudo_bayes_factors %<>% as.data.frame()
disease_pseudo_bayes_factors_sex0 %<>% as.data.frame()
disease_pseudo_bayes_factors_sex1 %<>% as.data.frame()
colnames(disease_pseudo_bayes_factors) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")
colnames(disease_pseudo_bayes_factors_sex0) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")
colnames(disease_pseudo_bayes_factors_sex1) <- c("order1_vs_order0", "order2_vs_order0", "order3_vs_order0", "order2_vs_order1", "order3_vs_order1", "order3_vs_order2")


pseudo_bayes_factor <- rbind(control_pseudo_bayes_factors, control_pseudo_bayes_factors_sex0)
pseudo_bayes_factor %<>% rbind(., control_pseudo_bayes_factors_sex1)
pseudo_bayes_factor %<>% rbind(., disease_pseudo_bayes_factors)
pseudo_bayes_factor %<>% rbind(., disease_pseudo_bayes_factors_sex0)
pseudo_bayes_factor %<>% rbind(., disease_pseudo_bayes_factors_sex1)

group <- c("all_control", "female_control", "male_control", "all_disease", "female_disease", "male_disease")
pseudo_bayes_factor %<>% mutate(group=group) %>% select(c(7,1:6))

pseudo_bayes_factor %>%
  write.csv(., "log_pseudo_bayes_factors_resolution_0_64_seconds.csv")
View(control_pseudo_bayes_factors)
View(disease_pseudo_bayes_factors)

temp_control_fit_data <- control_plot_data %>%
  filter(markov_order == 0 | markov_order == 1) %>%
  mutate(markov_order=as.factor(markov_order), Sex=as.factor(Sex))

lmFit <- lm(log_likelihood ~ mouse_id + Sex + markov_order + Sex:markov_order, data = temp_control_fit_data)
summary(lmFit)


temp_disease_fit_data <- disease_plot_data %>%
  filter(markov_order == 0 | markov_order == 2) %>%
  mutate(markov_order=as.factor(markov_order), Sex=as.factor(Sex))

lmFit <- lm(log_likelihood ~ mouse_id + Sex + markov_order + Sex:markov_order, data = temp_disease_fit_data)
summary(lmFit)
