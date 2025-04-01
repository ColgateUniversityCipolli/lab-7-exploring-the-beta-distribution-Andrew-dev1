library(tidyverse)
library(e1071)
library(patchwork)
library(cumstats)
library(nleqslv)
####################################################
# Plot the beta distributions and calculate variables
####################################################

calculate <- function(alpha, beta){
  beta.mean <- alpha/ (alpha +beta)
  beta.variance <- (alpha*beta)/ ((alpha +beta)^2 *(alpha+beta+1))
  beta.skewness <- (2*(beta-alpha)* (alpha+beta+1)^(1/2)) /
    ((alpha+beta+2)*(alpha*beta)^(1/2))
  beta.kurtosis <-6*((alpha -beta)^2 *(alpha+beta+1) - ((alpha*beta)*(alpha+beta+2)))/ 
    ((alpha*beta)*(alpha+beta+2)*(alpha+beta+3))
  
  dat.comparison <- tibble() %>%
    summarize(variable = sprintf("Beta(%g,%g)", alpha,beta),
              mean = beta.mean,
              variance = beta.variance,
              skewness = beta.skewness,
              excess.kurtosis = beta.kurtosis
    )  
  return(dat.comparison)
}

betas <- tibble(alpha = c(2,5,5,0.5), beta = c(5,5,2,0.5))
value.dat <- betas %>%
  pmap_dfr(~calculate(..1, ..2))

# generate a grid of points and calculate beta values for different alphas and betas
sequence <- seq(0, 1.0, length.out=500)
figure1.data <- tibble(x = sequence)|>   
  mutate(beta.2.5 = dbeta(x, 2, 5),
         beta.5.5 = dbeta(x, 5, 5),
         beta.5.2 = dbeta(x, 5, 2),
         beta.half = dbeta(x, 0.5, 0.5)
         )

# plot beta distributions
beta.distributions <- ggplot(data= figure1.data)+                     # specify data
  geom_line(aes(x=x, y=beta.2.5, color="Beta(2,5)")) +                
  geom_line(aes(x=x, y=beta.5.5, color="Beta(5,5)")) + 
  geom_line(aes(x=x, y=beta.5.2, color="Beta(5,2)")) + 
  geom_line(aes(x=x, y=beta.half, color="Beta(0.5,0.5)")) + 
  geom_hline(yintercept=0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "lightblue", "purple",
                                    "orange"))+                 # change colors
  theme(legend.position = "bottom")


# Task Two, moment function was created to bypass the unknown x variable 
beta.moment <- function(alpha, beta, k, centered){
  moment <- function(x) {
    x^k * dbeta(x, shape1 = alpha, shape2 = beta)
  }
  value <- integrate(moment, lower = 0, upper = 1)$value
  if(centered){
    center.moment <- function(x) {
      (x-value)^k * dbeta(x, shape1 = alpha, shape2 = beta)
    }
    return(integrate(center.moment, lower = 0, upper = 1)$value)
  }
  return(value)
}

beta.moment(2,3,4,T)


# Task three: create random beta samples to view distributions
set.seed(7272)

sample.size <- 500
beta.sample1 <- rbeta(n = sample.size,  # sample size
                     shape1 = 2,   # alpha parameter
                     shape2 = 5)    # beta parameter

sample.df1 <- as.data.frame(beta.sample1)
betaplot1 <- ggplot(sample.df1, aes(x= beta.sample1))+
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(from = 0, to = 1, by = 0.075))


beta.sample2 <- rbeta(n = sample.size,  # sample size
                      shape1 = 5,   # alpha parameter
                      shape2 = 5)    # beta parameter
sample.df2 <- as.data.frame(beta.sample2)
betaplot2 <- ggplot(sample.df2, aes(x= beta.sample2))+
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(from = 0, to = 1, by = 0.075))


beta.sample3 <- rbeta(n = sample.size,  # sample size
                      shape1 = 5,   # alpha parameter
                      shape2 = 2)    # beta parameter
sample.df3 <- as.data.frame(beta.sample3)
betaplot3 <- ggplot(sample.df3, aes(x= beta.sample3))+
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(from = 0, to = 1, by = 0.075))

beta.sample4 <- rbeta(n = sample.size,  # sample size
                      shape1 = 0.5,   # alpha parameter
                      shape2 = 0.5)    # beta parameter
sample.df4 <- as.data.frame(beta.sample4)
betaplot4 <- ggplot(sample.df4, aes(x= beta.sample4))+
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(from = 0, to = 1, by = 0.075))

(betaplot1 + betaplot2)
(betaplot3 + betaplot4)

calculate2 <- function(beta.sample, num1, num2){
  tibs <- tibble()|>
    summarize(variable = sprintf("Sample Beta(%g,%g)", num1, num2),
            mean = mean(beta.sample),
            variance = var(beta.sample),
            skewness = skewness(beta.sample),
            excess.kurtosis = kurtosis(beta.sample)
    )
  return(tibs)
}
samples <- list(
  list(beta.sample1, 2, 5),
  list(beta.sample2, 5, 5),
  list(beta.sample3, 5, 2),
  list(beta.sample4, 0.5, 0.5)
)

sample.values.dat <- map_dfr(samples, ~ calculate2(.x[[1]], .x[[2]], .x[[3]]))
comparison <- merge(value.dat, sample.values.dat, all.y = T, all.x = T)


## Task four 

dat.beta.2.5 <- tibble(x = beta.sample1)|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, shape1 = 2, shape2 = 5),
         cummulative.mean = cummean(x),
         cummulative.variance = cumvar(x),
         cummulative.skewness = cumskew(x),
         cummulative.kurtosis = cumkurt(x)
         )



mean.plot <- ggplot()+
  geom_line(data = dat.beta.2.5, aes(x = 1:500, y = cummulative.mean))+
  geom_hline(yintercept = comparison$mean[2])+
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Cumulative means")

variance.plot <- ggplot()+
  geom_line(data = dat.beta.2.5, aes(x = 1:500, y = cummulative.variance))+
  geom_hline(yintercept = comparison$variance[2])+
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("variances")

skewness.plot <- ggplot()+
  geom_line(data = dat.beta.2.5, aes(x = 1:500, y = cummulative.skewness))+
  geom_hline(yintercept = comparison$skewness[2])+
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Skewnesses")

kurtosis.plot <- ggplot()+
  geom_line(data = dat.beta.2.5, aes(x = 1:500, y = cummulative.kurtosis))+
  geom_hline(yintercept = comparison$excess.kurtosis[2])+
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Kurts")
 
(mean.plot + variance.plot) / (skewness.plot + kurtosis.plot)



for(i in 2:50){
  set.seed(7272 + i)
  beta.sample <- rbeta(n = sample.size,  # sample size
                        shape1 = 2,   # alpha parameter
                        shape2 = 5)
  loops.data <- tibble(x = beta.sample)|>   # generate a grid of points
    mutate(beta.pdf = dbeta(x, shape1 = 2, shape2 = 5),
           cummulative.mean = cummean(x),
           cummulative.variance = cumvar(x),
           cummulative.skewness = cumskew(x),
           cummulative.kurtosis = cumkurt(x)
    )  
  mean.plot <- mean.plot +
    geom_line(data = loops.data, aes(x = 1:500, y = cummulative.mean), color= i)
  variance.plot <- variance.plot +
    geom_line(data = loops.data, aes(x = 1:500, y = cummulative.variance), color= i)
  skewness.plot <- skewness.plot +
    geom_line(data = loops.data, aes(x = 1:500, y = cummulative.skewness), color= i)
  kurtosis.plot <- kurtosis.plot +
    geom_line(data = loops.data, aes(x = 1:500, y = cummulative.kurtosis), color= i)
  
}

(mean.plot + variance.plot) / (skewness.plot + kurtosis.plot)

##Task 5 
summary <- tibble()
for(i in 1:1000){
  set.seed(7272 + i)
  beta.sample <- rbeta(n = sample.size, 
                       shape1 = 2,   
                       shape2 = 5)
  summary <- rbind(summary, calculate2(beta.sample, 0, i))
  
}
ggplot(summary)+
  geom_histogram(aes(x = mean, y = after_stat(density)))+
  stat_density(aes(x = mean), geom="line")+
  theme_bw()

ggplot(summary)+
  geom_histogram(aes(x = variance, y = after_stat(density)))+
  stat_density(aes(x = variance), geom="line")+
  theme_bw()


# ggplot(summary)+
#   geom_histogram(aes(x = variance, y = after_stat(density)))
# ggplot(summary)+
#   geom_histogram(aes(x = skewness, y = after_stat(density)))
# ggplot(summary)+
#   geom_histogram(aes(x = excess.kurtosis, y = after_stat(density)))
#############################
# Task Six
#############################

deaths.data <- read_csv(file = "Global deaths from World Bank/total_deaths.csv",
                        skip = 4) |>
  select(c(1,2,3,4,67)) |>
  rename("y2022" = "2022") |>
  mutate(y2022 = y2022 /1000) |>
  mutate(y2022 = case_when(is.na(y2022) ~ 0,
                   TRUE ~ y2022))

ggplot(deaths.data)+
  geom_histogram(aes(x =  y2022, y = after_stat(density)),
                 breaks = seq(from = 0.000, to = 0.024, by = 0.001),
                 color = "black")+
  stat_density(aes(x = y2022), geom="line")+
  theme_bw()

MOM.beta <- function(data, par){
  alpha <- par[1]
  beta <- par[2]
  
  EX1 <- alpha/ (alpha +beta)
  EX2 <- (alpha * (alpha +1))/((alpha + beta + 1) * (alpha +beta))
  
  m1 <- mean(data)
  m2 <- mean(data^2)
  
  return( c(EX1 - m1, EX2 - m2) )
}

moms<- nleqslv(x = c(1,1),
               fn = MOM.beta,
               data=deaths.data$y2022)

alpha.moms <- moms$x[1]
beta.moms <- moms$x[2]

moms.internal <- MOM.beta(deaths.data$y2022, c(1,1))

MLE.beta <- function(data, par, neg = F){
  alpha <- par[1]
  beta <- par[2]
  ll <- sum(log(dbeta(x = data, shape1 = alpha, shape2 = beta)))
  
  return(ifelse(neg, -ll, ll))
}

mles <- optim(par = c(1, 1),
              fn = MLE.beta,
              data=deaths.data$y2022,
              neg=T)

alpha.mle <- mles$par[1]
beta.mle <- mles$par[2]

alpha = 8
beta = 950
result <- tibble()
for(i in 1:1000){
  set.seed(7272+i)
  beta.sample <- rbeta(n = 266,      # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta) # beta parameter

  moms.internal <- nleqslv(x = c(1,1),
                           fn = MOM.beta,
                           data= beta.sample)
  mles.internal <- optim(par = c(1,1),
                         fn = MLE.beta,
                         data= beta.sample,
                         neg=T)
  result <- bind_rows(result, tibble(moms.alpha = moms.internal$x[1],
                                     moms.beta = moms.internal$x[2],
                                     mles.alpha = mles.internal$par[1],
                                     mles.beta = mles.internal$par[2]))
}

plot.moms.alpha = ggplot(data = result, aes(x = moms.alpha)) + 
  geom_density(fill = "#FD8A8A")+ 
  theme_bw() +
  geom_hline(yintercept=0)+
  xlab("alpha values from moms")+
  ggtitle("Estimated density for moms")

plot.moms.beta = ggplot(data = result, aes(x = moms.beta)) + 
  geom_density(fill = "#A8D1D1")+ 
  theme_bw() +
  geom_hline(yintercept=0)+
  xlab("beta values from moms")+
  ggtitle("Estimated density for moms")

plot.mles.alpha = ggplot(data = result, aes(x = mles.alpha)) + 
  geom_density(fill = "#d5d1e9")+ 
  theme_bw() +
  geom_hline(yintercept=0)+
  xlab("alpha values from mles")+
  ggtitle("Estimated density for mles")

plot.mles.beta = ggplot(data = result, aes(x = mles.beta)) + 
  geom_density(fill = "#f5cf9f")+ 
  theme_bw() +
  geom_hline(yintercept=0)+
  xlab("beta values from mles")+
  ggtitle("Estimated density for mles")

(plot.moms.alpha + plot.moms.beta) / (plot.mles.alpha + plot.mles.beta)

values = function(theta.hats, theta){
  (bias <- mean(theta.hats) - theta)
  (precision <- 1/var(theta.hats))
  (mse <- var(theta.hats) + bias^2)
  return(tibble(bias, precision, mse))
}

a.moms <- values(result$moms.alpha, alpha)
b.moms <- values(result$moms.beta, beta)
a.mles <- values(result$mles.alpha, alpha)
b.mles <- values(result$mles.beta, beta)

accuracy <- tibble()|> 
  bind_rows(a.moms, b.moms, a.mles, b.mles) |> 
  mutate(names = c("moms alpha", "moms beta", "mles alpha", "mles beta"),
        actual = c(alpha, beta, alpha, beta),
        estimated = c(mean(result$moms.alpha), mean(result$moms.beta), 
                      mean(result$mles.alpha), mean(result$mles.beta))
         )

accuracy
