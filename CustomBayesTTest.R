####################################################################################
# comparing the traditional t-test with a custom bayes factor test and
# original bayes factor t-test from the BayesFactor package
####################################################################################

#' @title Bayesian Alternative Model
#' 
#' @description The general form of the alternative model, contrary to the null model.
#' 
#' @param t The value of your observed test statistic T = t
#' @param standard_error_variance The variance of the standard error. In the Gaussian context of comparing two means, the standard error would be sqrt((s_x^2 + s_y^2)/n) (assuming an equal sample size of n). Thus, the variance would just be (s_x^2 + s_y^2) / n
#'
#' @returns The likelihood of the observed test statistic T = t under the alternative model. This can be used to obtain the bayes factor by dividing against the likelihood under the null model
alt_model = function(t, standard_error_variance)
{
  return(exp(-(t^2)/(4*standard_error_variance)) / (2*sqrt(3.14159)*sqrt(standard_error_variance)))
}



#' @title Bayesian Null Model
#' 
#' @description The (asymptotic) sampling distribution of the test statistic given the null hypothesis is true.
#' 
#' @param t The value of your observed test statistic T = t
#' @param standard_error_variance The variance of the standard error. In the Gaussian context of comparing two means, the standard error would be sqrt((s_x^2 + s_y^2)/n) (assuming an equal sample size of n). Thus, the variance would just be (s_x^2 + s_y^2) / n.
#' 
#' @returns The likelihood of the observed test statistic T = t under the null model. This is the denominator when computing the bayes factor.
null_model = function(t, standard_error_variance)
{
  return(dnorm(t, mean = 0, sd = sqrt(standard_error_variance)))
}



#' @title Bayesian T-Test
#' 
#' @description This is a bayesian variant of the t-test which mergest frequentist and bayesian ideas together. We establish the null model as the asymptotic distribution of the test statistic given the null hypothesis is true (so far this is the same as frequentist statistics). However, when we establish the alternative model, we allow the mean to vary according to the standard error, i.e, the prior distribution of the mean is centered at zero and varies according to the standard error of the test stistic. The intuition behind this is that when the standard error grows, the variation in the effect size will also grow and vice versa when it shrinks. Thus, we make the variation of the mean to be dependent on the standard error.
#' 
#' @param g1 Group 1
#' @param g2 Group 2 (assuming equal sample size to Group 1)
#' 
#' @returns The bayes factor which can be used to construct a decision rule, such as if bayes factor is > 3, then we "reject" the null model
bayes_ttest = function(g1, g2)
{
  # general variance assuming g1 and g2 are independent
  standard_error_variance = (var(g1) + var(g2)) / length(g1)
  t = mean(g1) - mean(g2)
  return(alt_model(t, standard_error_variance)/null_model(t, standard_error_variance))
}



# simulations (experiment with different means, variances, sample sizes)
library(BayesFactor)

trad_tt = c()
bf_cust = c()
bf_orig = c()
for (i in 1:10000)
{
  x = rnorm(30, mean = 2, sd = 2)
  y = rnorm(30, mean = 3, sd = 2)
  trad_tt[i] = ifelse(t.test(x,y)$p.value < 0.01, 1, 0)
  
  # custom bayes factor test
  bf_cust[i] = ifelse(bayes_ttest(x, y) > 3, 1, 0)
  
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
