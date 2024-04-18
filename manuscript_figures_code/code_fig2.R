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
g_1 = -6              # effect of treatment on animal welfare
g_2 = -1              # effect of initial infarct size on welfare (L -> W)
sd_Y = 10             # standard deviation of epsilon for final infarct size Y  
sd_W = 2              # standard deviation of epsilon for animal welfare W 
cutoff_W = 0.25       # attrition rates 
n = 20                # total number of animals in the experiment


# combining all possible parameter values and values of n 
report <- expand_grid(prob_A, mean_L, 
                      sd_L, b_0, b_1, 
                      b_2, g_0, g_1, 
                      g_2, sd_Y, sd_W, 
                      cutoff_W, n)
# setting number of random draws 
B = 1

######################################################################


######################################################################
# simulation of B random draws with set parameter values
########################################################


# ensure reproducibility of random draw 
set.seed(120) 
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
    dat_s <- dat %>% filter(S==TRUE)  			
  }
}



#########################################################
#mutating data to prepare for visualization
dat <- dat %>%
  mutate(
    cols = factor(dat$A, labels = c(cols[1], cols[4])),
    S1 = factor(dat$S, labels = c(0, 1)),
    opacity = recode(S1, `0` = 0.3, `1` = 1),
    X1 = factor(recode(A, `0` = "Control", `1` = "Treatment"))
  )

dat_s <- 
  dat_s %>%
  mutate(
    cols = factor(dat_s$A, labels = c(cols[1], cols[4])),
    X1 = factor(recode(A, `0` = "Control", `1` = "Treatment"))
  )

dat_summary <- 
  dat %>% 
  group_by(A) %>%
  summarize(
    mean_L = mean(L),
    mean_Y = mean(Y),
    sd_L = sd(L),
    sd_Y = sd(Y)
  ) %>%
  mutate(
    cols = c(cols[1], cols[4]),
    xend = c(1.3, 2.3),
    x = c(0.7, 1.7)
  )

dat_s_summary <- 
  dat_s %>% 
  group_by(A) %>%
  summarize(
    mean_L = mean(L),
    mean_Y = mean(Y),
    sd_L = sd(L),
    sd_Y = sd(Y)
  ) %>%
  mutate(
    cols = c(cols[1], cols[4]),
    xend = c(1.3, 2.3),
    x = c(0.7, 1.7)
  )

#######################################################



######################################################
# visualization for Figure 2

letter <- theme(plot.title = element_text(face = "bold", 
                                          size = 14), 
                axis.title = element_text(face = "bold", 
                                          size = 14), 
                axis.text = element_text(face = "bold", 
                                         size = 12)
)


fig2a <- 
  ggplot() +
  geom_quasirandom(
    data = dat,
    aes(y = L, x = X1, fill = cols, alpha = opacity),
    shape = 21,
    color = "darkblue",
    size = 5,
    width = 0.3
  ) +
  scale_fill_manual(values = c(cols[1], cols[4]))+
  geom_segment(
    data = dat_summary,
    aes(
      x = x,
      y = mean_L,
      yend = mean_L,
      xend = xend, 
    ),
    color = dat_summary$cols,
    linewidth = 1,
    alpha = 0.3
  ) +
  geom_segment(
    data = dat_s_summary,
    aes(
      x = x,
      y = mean_L,
      yend = mean_L,
      xend = xend, 
    ),
    color = dat_summary$cols,
    linewidth = 1,
    alpha = 1
  ) +
  ylab(expression("Initial infarct size (mm" ^ 3 * ")")) +
  xlab("") +
  ylim(10, 40) +
  ggtitle("A") +
  theme_classic()+
  letter + 
  theme(legend.position = "none", 
        plot.margin = unit(c(2,0,2,2), "cm"),
  )




plot(fig2a)


fig2b <- 
  ggplot() +
  geom_quasirandom(
    data = dat,
    aes(y = Y, 
        x = X1, 
        group = X1, 
        fill = cols, 
        alpha = opacity),
    shape = 21,
    color = "darkblue",
    size = 5,
    width = 0.3
  ) +
  scale_fill_manual(values = c(cols[1], cols[4]))+
  geom_segment(
    data = dat_summary,
    aes(
      x = x,
      y = mean_Y,
      yend = mean_Y,
      xend = xend
    ),
    color = dat_summary$cols,
    linewidth = 1,
    alpha = 0.3
  ) +
  geom_segment(
    data = dat_s_summary,
    aes(
      x = x,
      y = mean_Y,
      yend = mean_Y,
      xend = xend
    ),
    color = dat_summary$cols,
    linewidth = 1,
    alpha = 1
  ) +
  ylab(expression("Final infarct size (mm" ^ 3 * ")")) +
  xlab("") +
  ylim(120, 300) +
  ggtitle("B") +
  theme_classic()+
  letter + 
  theme(legend.position = "none", 
        plot.margin = unit(c(2,2,2,0), "cm"),
  )


plot(fig2b)

fig2 <- grid.arrange(
  
  fig2a, 
  fig2b,
  
  nrow = 1,
  ncol = 2
)


fig2


