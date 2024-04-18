source("load_packages.R")

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
g_0 = 0               # intercept value of animal welfare W 
g_1 = c(-1,-3,-6)     # effect of treatment on animal welfare
g_2 = -1              # effect of initial infarct size on welfare (L -> W)
sd_W = 2              # standard deviation of epsilon for animal welfare W 
cutoff_W = c(0.1, 0.25,0.5)    # attrition rates 



# combining all possible parameter values and values of n 
report <- expand_grid(prob_A, mean_L, 
                      sd_L, g_0, g_1, 
                      g_2, sd_W, 
                      cutoff_W)

# METHOD 1: analytical values

report$prob_treated_analytical <- NA 
report$prob_untreated_analytical <- NA

for (i in 1:nrow(report)){
  mean1 = report$g_0[i] + report$g_1[i] + report$g_2[i]*report$mean_L[i] 
  #mean welfare score in the treated
  mean2 = report$g_0[i] + report$g_2[i]*report$mean_L[i] 
  #mean welfare score in the non-treated
  sd = sqrt((report$sd_W[i])^2+((report$g_2[i])^2)*(report$sd_L[i])^2)
  # calculated sd
  
  quant_select <- qmixnorm(p = report$cutoff_W[i], mean = c(mean1, mean2), sd = c(sd, sd), pro = c(0.5, 0.5)) 
  # thresd on the welfare score that leaves the specified cutoff_W % of animals out
  report$prob_treated_analytical[i] <- pnorm(quant_select, mean1, sd = sd) 
  # probability of being out if you are in the treated group
  report$prob_untreated_analytical[i] <- pnorm(quant_select, mean2, sd = sd) 
  # probability of being out if you are in the non-treated group
  
  rm(mean1, mean2, sd, quant_select)
}

# METHOD 2: approximation by simulation

# setting sample size and number of draws
n = 100000
B=100
# ensure reproducibility of random draw 
set.seed(10000)

report$out_treated <- NA 
report$out_non_treated <- NA

# here starts the counting of the row
for (i in 1:nrow(report)) {
  
  out_treat <- rep(NA_real_, B)
  out_non_treat <- rep(NA_real_, B)
  
  # here starts the loop in b
  for (b in 1:B) {
    
    # create the dataset
    dat <- data.frame(A = rep(0:1,report$prob_A[i]*n, each=1),
                      L = rnorm(n=n, mean=report$mean_L[i], sd=report$sd_L[i]))
    dat$W <- report$g_0[i] + report$g_1[i]*dat$A + report$g_2[i]*dat$L + rnorm(n, 0, report$sd_W[i])
    dat$S <- dat$W >= quantile(dat$W, probs=report$cutoff_W[i]) 
    
    out_treat[b] <- mean(dat[dat$A==1,"S"]==FALSE)
    out_non_treat[b] <- mean(dat[dat$A==0,"S"]==FALSE)
  }
  report$out_treated[i] <- mean(out_treat)
  report$out_non_treated[i] <- mean(out_non_treat)  
}


kable(report) %>% kable_classic()

