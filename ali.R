library(tidyverse)
library(e1071)
library(patchwork)
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

value.dat2 <- betas %>%
  pmap_dfr(~calculate(..1, ..2))
figure1.data <- tibble(x = seq(0, 1.0, length.out=500))|>   # generate a grid of points
  mutate(beta.2.5 = dbeta(x, 2, 5),
         beta.5.5 = dbeta(x, 5, 5),
         beta.5.2 = dbeta(x, 5, 2),
         beta.half = dbeta(x, 0.5, 0.5)
         )
beta.distributions <- ggplot(data= figure1.data)+                                            # specify data
  geom_line(aes(x=x, y=beta.2.5, color="Beta(2,5)")) +                 # plot beta dist
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
    x^n * dbeta(x, shape1 = alpha, shape2 = beta)
  }
  value <- integrate(moment, lower = 0, upper = 1)$value
  if(centered){
    center.moment <- function(x) {
      (x-value)^n * dbeta(x, shape1 = alpha, shape2 = beta)
    }
    return(integrate(center.moment, lower = 0, upper = 1)$value)
  }
  return(value)
}

beta.moment(2,3,4,T)

# Task three
set.seed(7272)

sample.size <- 500
beta.sample1 <- rbeta(n = sample.size,  # sample size
                     shape1 = 2,   # alpha parameter
                     shape2 = 5)    # beta parameter
sample.df1 <- as.data.frame(beta.sample1)
betaplot1 <- ggplot(sample.df1, aes(x= beta.sample1))+
  geom_histogram(aes(y = after_stat(density)))


beta.sample2 <- rbeta(n = sample.size,  # sample size
                      shape1 = 5,   # alpha parameter
                      shape2 = 5)    # beta parameter
sample.df2 <- as.data.frame(beta.sample2)
betaplot2 <- ggplot(sample.df2, aes(x= beta.sample2))+
  geom_histogram(aes(y = after_stat(density)))


beta.sample3 <- rbeta(n = sample.size,  # sample size
                      shape1 = 5,   # alpha parameter
                      shape2 = 2)    # beta parameter
sample.df3 <- as.data.frame(beta.sample3)
betaplot3 <- ggplot(sample.df3, aes(x= beta.sample3))+
  geom_histogram(aes(y = after_stat(density)))

beta.sample4 <- rbeta(n = sample.size,  # sample size
                      shape1 = 0.5,   # alpha parameter
                      shape2 = 0.5)    # beta parameter
sample.df4 <- as.data.frame(beta.sample4)
betaplot4 <- ggplot(sample.df4, aes(x= beta.sample4))+
  geom_histogram(aes(y = after_stat(density)))

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