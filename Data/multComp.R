rm(list = ls())

## Demonstrate multiple testing

simple_t_test <- function(){
  x <- runif(100); y <- runif(100)
  res <- t.test(x,y)
  return(c(res$statistic, res$p.value))
}

hist( replicate(n = 10000, simple_t_test()[1]),
      main = "Distribution of Differences",
      xlab = "Difference Between X and Y",
      ylab = "Frequency")

hist( replicate(n = 10000, simple_t_test()[2]),
      main = "Distribution of p-values",
      xlab = "p-value")
abline(v = 0.05, lwd = 10, col = "red", lty = "dashed")


numFalsePos <- c(sum( replicate(n = 1, simple_t_test()[2])<0.05),
  sum( replicate(n = 10, simple_t_test()[2])<0.05),
  sum( replicate(n = 100, simple_t_test()[2])<0.05),
  sum( replicate(n = 1000, simple_t_test()[2])<0.05),
  sum( replicate(n = 10000, simple_t_test()[2])<0.05))

numTests <- c(1,10,100,1000,10000)

 