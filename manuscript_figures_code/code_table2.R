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
n = 10000000


# combining all possible parameter values and values of n 
report <- expand_grid(prob_A, mean_L, 
                      sd_L, g_0, g_1, 
                      g_2, sd_W, 
                      cutoff_W, n)

# setting number of random draws
B = 100

# ensure reproducibility of random draw 
set.seed(10000)
######################################################################
# function for replication draws

report$out_treated <- NA 
report$out_non_treated <- NA

# here starts the counting of the row
for (i in 1:nrow(report)) {
  
  out_treat <- rep(NA_real_, B)
  out_non_treat <- rep(NA_real_, B)
  
  # here starts the loop in b
  for (b in 1:B) {
    
    # create the dataset
    dat <- data.frame(A = rep(0:1,report$prob_A[i]*report$n[i], each=1),
                      L = rnorm(n=report$n[i], mean=report$mean_L[i], 
                                sd=report$sd_L[i]))
    dat$W <- report$g_0[i] + report$g_1[i]*dat$A + report$g_2[i]*dat$L + rnorm(report$n[i], 0, report$sd_W[i])
    dat$S <- dat$W >= quantile(dat$W, probs=report$cutoff_W[i]) 
    
    out_treat[b] <- mean(dat[dat$A==1,"S"]==FALSE)
    out_non_treat[b] <- mean(dat[dat$A==0,"S"]==FALSE)
  }
  report$out_treated[i] <- mean(out_treat)
  report$out_non_treated[i] <- mean(out_non_treat)  
}


kable(round(report,2)) %>% kable_classic()

