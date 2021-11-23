####################################################################################
# comparing the traditional t-test with a custom bayes factor test and
# original bayes factor t-test from the BayesFactor package
####################################################################################

# install.packages("BayesFactor")
library(BayesFactor)

trad_tt = c() # traditional t-test
bf_cust = c() # my custom bayesian test
bf_orig = c() # the original bayesian test from BayesFactor package

# change the values of the mean to make the null hypothesis true/false

# From my observations, the bayesian test(s) are pretty robust to false positives
# and provide a good balance of true negatives and true positives (contrary to choosing
# a 5% or 1% alpha, where 5% alphas can have higher false positives and 1% alphas can lack power
# in small-sample-moderate-variance cases like below)
for (i in 1:10000)
{
  x = rnorm(60, mean = 2, sd = 2)
  y = rnorm(60, mean = 3, sd = 2)
  trad_tt[i] = ifelse(t.test(x,y)$p.value < 0.01, 1, 0)
  
  # custom bayes factor test (specific to this standard error (e.g sample size 60 and total variance of 8))
  # if you change the sample size and/or variance, you MUST re-derive the alternate model!
  # perhaps in the future I will derive a general form for this and create a function
  tstat = mean(x) - mean(y)
  null_model = dnorm(tstat, mean = 0, sd = sqrt(8/60))
  alt_model = (sqrt(15)*exp((-15*tstat^2)/8)) / (2^(3/2)*sqrt(3.14159))
  
  bf_cust[i] = ifelse(alt_model/null_model > 3, 1, 0)
  
  # R bayes factor test
  bf_orig[i] = ifelse(extractBF(ttestBF(x,y))$bf > 3, 1, 0)
  
}

mean(trad_tt)
mean(bf_cust)
mean(bf_orig)

#########################################################################################################
# how to derive the custom bayes test? It's a hybrid approach of frequentist and Bayesian statistics!
#
# Let T(x,y) be the test statistic of interest. Then compute its standard error SE(T).
#
# The null model will be assuming the test statistic is centered at 0 with standard error SE(T).
#
# The alternate model will allow the mean to vary according to the standard error,
# i.e, the prior distribution of the mean will have standard deviation of SE(T).
#
# The intuition behind this is, as the standard error grows, the variation of the effect size also 
# grows, so we allow the prior distribution to be determined by the size of the standard error.
#
# To finally derive the alternate model, we take the marginal:
# p(x) = int_{-infty}^{infty} p(x|mu)p(mu)dmu, p(x) is now our alternate model.
#
# We then take the bayes factor given a value of our test statistic T = t:
# BF = alt_model(t) / null_model(t)
#
# and define some decision rule (typically 3) to determine if we "reject" the null model. I.e,
# if BF > 3 then we "reject" the null, otherwise, we lack sufficient evidence that an effect is present.
#########################################################################################################
