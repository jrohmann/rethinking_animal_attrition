source("load_packages.R")

# settings for styling visualization
cols <- RColorBrewer::brewer.pal(6, "Dark2")

#############################################################
# setting parameter values 
#############################################################
# A alternate treatments with proportion prob_A of treated
# L is normal with mean_L and sd_L
# Y = b_0 + b_1*A + b_2*L + rnorm(0,sd_Y)
# W = g_0 + g_1*A + g_2*L + rnorm(0, sd_W)


prob_A = 0.5          # proportion of animals assigned to treatment
mean_L = 25           # mean of initial infarct size 
sd_L = 5              # standard deviation of initial infarct size
b_0 = 0               # intercept value for final infarct size 
b_1 = 0               # treatment effect - here null-effect 
b_2 = 8               # effect of initial infarct size on final infarct size 
g_0 = 0               # intercept value of animal welfare W 
g_1 = c(-1,-3,-6)     # effect of treatment on animal welfare
g_2 = -1              # effect of initial infarct size on welfare (L -> W)
sd_Y = 10             # standard deviation of epsilon for final infarct size Y  
sd_W = 2              # standard deviation of epsilon for animal welfare W 
cutoff_W = c(0.1,0.25,0.5)    # attrition rates 
n = c(10, 20, 50)     # total number of animals in the experiment


# combining all possible parameter values and values of n 
report <- expand_grid(prob_A, mean_L, 
                      sd_L, b_0, b_1, 
                      b_2, g_0, g_1, 
                      g_2, sd_Y, sd_W, 
                      cutoff_W, n)
# setting number of random draws
B = 10000

######################################################################


######################################################################
# simulation of B random draws with set parameter values


# The simulation will take considerable time to run 
# If you prefer to run the analysis with the dataset 
# given please load the simulated dataset with 
load("df_simulation_fig3.RData")
  
# proceed then with "data wrangling" 
########################################################


# preparation of empty matrix for model 1 ("oracle"), model 2 ("naive") & model 3 ("adjusted") 
# for later storage of values from replications 
oracle <- matrix(NA_real_, ncol=B, nrow=nrow(report))
naive <- matrix(NA_real_, ncol=B, nrow=nrow(report))
adjusted <- matrix(NA_real_, ncol=B, nrow=nrow(report))


# ensure reproducibility of random draw 
set.seed(10000)
######################################################################
# function for replication draws 

# here starts the counting of the row
for (i in 1:nrow(report)) {
  
  # here starts the loop in b
  for (b in 1:B) {
    
    # create the dataset
    dat <- data.frame(A = rep(0:1,report$prob_A[i]*report$n[i], each=1),
                      L = rnorm(n=report$n[i], mean=report$mean_L[i], 
                                sd=report$sd_L[i]))
    dat$W <- report$g_0[i] + report$g_1[i]*dat$A + report$g_2[i]*dat$L + rnorm(report$n[i], 0, report$sd_W[i])
    dat$Y <- report$b_0[i] + report$b_1[i]*dat$A + report$b_2[i]*dat$L + rnorm(report$n[i], 0, report$sd_Y[i])
    dat$S <- dat$W >= quantile(dat$W, probs=report$cutoff_W[i]) 
    dat_s <- dat %>% filter(S==TRUE)  		# filtering the surviving animals
    naive[i,b] =  mean(dat_s$Y[dat_s$A==1]) - mean(dat_s$Y[dat_s$A==0]) 	#NAIVE, Among selected animals only, 
                                                                          # mean of the outcome among treated minus mean of outcome among untreated
    
    oracle[i,b] = mean(dat$Y[dat$A==1]) - mean(dat$Y[dat$A==0])		# ORACLE, Among all animals, 
                                                                  # mean of the outcome among treated minus mean of outcome among untreated
    
    model<-lm(Y ~ A + L, data = dat_s)					  # regression for outcome only among selected animals, adjusting (L)
    adjusted[i,b]<-model$coefficients[2]					# 3 ADJUSTED, extracted coeff for A (saved)
  }
}

################################################################################

# For 27 scenarios (index = scenario)
# data transformation for stratification and visualization 

oracle_tibble <- as_tibble(oracle) %>% mutate(model = "oracle", index = 1:27)
naive_tibble <- as_tibble(naive) %>% mutate(model = "naive", index = 1:27)
adjusted_tibble <- as_tibble(adjusted) %>% mutate(model = "adjusted", index = 1:27)

df <- rbind(oracle_tibble, naive_tibble, adjusted_tibble)

# save(df, file="df_simulation_fig3.RData")

###############################################################################



#################################
# data wrangling 
#################################


names<-NULL

for (i in 1:B) {
  names[i] <- c(paste("estimate_",i, sep = ""))
}
colnames(df)<-c(names, "model", "index") 


#long format 
df2 <- gather(df, 
              key  = sim_number, 
              value = effect_estimate, 
              estimate_1:estimate_10000)

df3 <- report %>% 
  mutate(index = 1:nrow(report))%>%
  select(c(index,n, cutoff_W, g_1))

df3$g_1 <- factor(df3$g_1,levels = c(-1, -3, -6), 
                  labels = c("minor", "moderate", "major"))
df3$cutoff_W <- factor(df3$cutoff_W, levels = c(0.1,0.25,0.5),
                       labels = c("10", "25", "50"))
df4 <- inner_join(df2, df3)


# missing values 
df2 %>%
  group_by(model) %>%
  summarize(n.sum = n(), 
            n.na = sum(is.na(effect_estimate)),
            pct = n.na/n.sum) %>% ungroup()

na_naive <- 
  df2 %>%
  filter(model == "naive") %>%
  group_by(index) %>%
  summarize(n.sum = n(), 
            n.na = sum(is.na(effect_estimate)), 
            pct = n.na/n.sum) %>% ungroup()

na_adjusted <- 
  df2 %>%
  filter(model == "adjusted") %>%
  group_by(index) %>%
  summarize(n.sum = n(), 
            n.na = sum(is.na(effect_estimate)), 
            pct = n.na/n.sum) %>% ungroup()

na_naive <- 
  inner_join(na_naive, df3, by = "index")

na_adjusted <- 
  inner_join(na_adjusted, df3, by = "index")


###################################################
#
#   effect estimates for different values of n, 
#   attrition frequency and side effects 
#
##################################################



table_1<- df4%>%
  group_by(model, g_1, cutoff_W, n)%>%
  summarize(mean_bias = round(mean(effect_estimate, na.rm = TRUE),1),
            quan_25 = round(quantile(effect_estimate, probs = 0.025, na.rm = TRUE),1),
            quan_975 = round(quantile(effect_estimate, probs = 0.975, na.rm = TRUE),1)) %>%
  ungroup() %>% arrange(n)

colnames(table_1) <- c("approach","side effects", "attrition frequencies", "n(total)", 
                  "mean effect estimate", "2.5% quantile", "97.5% quantile")

kableExtra::kable(table_1, 
                  caption = "Bias in naive estimates by sample size\n, side effects and attrition frequencies") %>% 
  kable_classic()


###########################################################################################
# Visualization 
letter <- theme(plot.title = element_text(face = "bold", 
                                          size = 20, 
                                          margin = margin(2,0,0.5,0, unit = "cm")), 
                axis.text = element_text(face = "bold", 
                                         size = 12),
                axis.title.y = element_text(face = "bold", 
                                         size = 14, 
                                         margin = margin(0,0.5,0,2, unit = "cm")),
                axis.title.x = element_text(face = "bold", 
                                           size = 14, 
                                           margin = margin(0.5,0,0.5,0, unit = "cm")),
                plot.margin = unit(c(2,2,2,0), "cm"),
                panel.grid.major = element_blank(),
                legend.position = "bottom", 
                legend.title = element_blank(),
                legend.text = element_text(size = 12, 
                                           face = "bold"),
                strip.text = element_text(size = 12, 
                                          face = "bold"), 
                strip.background = element_rect(fill = "snow2"))




# Figure 3  violin plot 

df4$n <- as.factor(df4$n)
levels(df4$n) <- c(expression("n = 10"), expression("n = 20"), expression("n = 50"))
levels(df4$cutoff_W) <- c("attrition frequency = 10%", "attrition frequency = 25%", "attrition frequency = 50%")
df4$model <- factor(df4$model, levels = c("oracle", "naive", "adjusted"))

fig3<- ggplot(na.omit(df4))+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_violin(aes(x = g_1, y = effect_estimate, fill = model), alpha = 0.7)+
  geom_vline(xintercept = c(1.5, 2.5), 
             linewidth = 0.8)+
  facet_grid( n ~ cutoff_W)+
  coord_cartesian(ylim = c(-70,70))+
  scale_fill_manual(values = c(cols[1], cols[4], cols[8]), 
                    labels = c("oracle", "naive", "adjusted"))+
  labs(y = "estimated effect", 
       x = "strength of negative side effect of treatment on welfare")+
  theme_bw()+
  letter

fig3






