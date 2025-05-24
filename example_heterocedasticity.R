library(dplyr)
library(ggplot2)

set.seed(42)

n <- 100
x <- runif(n, 0, 10)

error <- rnorm(n, mean = 0, sd = 2 * x)

beta0 <- 2
beta1 <- 3

y <- beta0 + beta1 * x + error

plot(x, y)

model_fit <- lm(y ~ x)
summary(model_fit)


plot(predict(model_fit), rstudent(model_fit), 
     xlab = "Fitted values", ylab = "Studentised residuals")

ggplot() + 
  geom_point(aes(x, y)) +
  geom_smooth(aes(x, y), method = "lm", se = FALSE)
