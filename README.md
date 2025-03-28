# qgm: Quantile Graphical Models
This repository is based on the paper [Quantile Graph Discovery through QuACC: Quantile Association via Conditional Concordance](https://arxiv.org/abs/2411.17033). Take a look at it for details on motivation, underlying proofs, and real world results on large-scale biobanks.

# Motivation

Graphical structure learning is an effective way to assess and visualize cross-biomarker dependencies in biomedical settings. Standard approaches to estimating graphs rely on conditional independence tests that may not be sensitive to associations that manifest at the tails of joint distributions, i.e., they may miss connections among variables that exhibit associations mainly at lower or upper quantiles. 

We propose a novel measure of quantile-specific conditional association called QuACC: Quantile Association via Conditional Concordance. For a pair of variables and a conditioning set, QuACC quantifies agreement between the residuals from two quantile regression models, which may be linear or more complex, e.g., quantile forests. Using this measure as the basis for a test of null (quantile) association, we introduce a new class of quantile-specific graphical models. Through simulation we show our method is powerful for detecting dependencies that manifest at the tails of distributions. 


# Install and Run
If you have the `devtools` package installed, you can install the package from Github using the following commands:
```
devtools::install_github("zainakh/qgm")
library(qgm)
```

# Example

## Additional Dependencies (aka: `R` package managers are not fun) 
A quick demo of how to get started is as follows. Start by installing dependencies that aren't installed automatically by CRAN (as of Mar 2025 you still have to do this on R version 4.4.1). This is required to use the `pcalg` library.

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("graph", "RBGL"))

```

Lets load in dependencies for our data generating process and case study below:
 ```
library(pcalg)
library(igraph)
library(msm)
library(qgm)
library(quantreg)
```
## Create dataset that exhibits effects of interest
Create a dataset that has tail dependent effects.  
```
linear.quantile.dgp <- function(n=1000) {
  umat <- c()

  alpha <- runif(4, min=-0.4, max=0.4)
  alpha <- alpha * 0 # Can change this to have nonzero general trends

  beta.plus_minus <- rbinom(4, 1, 0.5)
  beta <- beta.plus_minus * runif(4, min = 0.3, max = 0.8) + (1 - beta.plus_minus) * runif(4, min = -0.8, max = -0.3)
  gamma.plus_minus <- rbinom(4, 1, 0.5)
  gamma <- gamma.plus_minus * runif(4, min = 0.3, max = 0.8) + (1 - gamma.plus_minus) * runif(4, min = -0.8, max = -0.3)

  z <- rtnorm(n, mean=0, sd=1, lower=-2, upper=2)
  y <- z * alpha[1] + z * beta[1] * ifelse(z >= quantile(z, 0.9), 1, 0) + z * gamma[1] * ifelse(z <= quantile(z, 0.1), 1, 0) + rnorm(n)
  x <- z * alpha[2] + z * beta[2] * ifelse(z >= quantile(z, 0.9), 1, 0) + z * gamma[2] * ifelse(z <= quantile(z, 0.1), 1, 0) + rnorm(n)
  wx <- x * alpha[3] + x * beta[3] * ifelse(x >= quantile(x, 0.9), 1, 0) + x * gamma[3] * ifelse(x <= quantile(x, 0.1), 1, 0)
  wy <- y * alpha[4] + y * beta[4] * ifelse(y >= quantile(y, 0.9), 1, 0) + y * gamma[4] * ifelse(y <= quantile(y, 0.1), 1, 0)
  w <- wx + wy + rnorm(n)

  df <- data.frame(cbind(z, y, x, w))
  names(df) <- c("Z", "Y", "X", "W")
  return(df)
}

df <- linear.quantile.dgp()
head(df)
 ```
 
 If we plot `Y` against `Z`, we can see that at the upper tail of `Z`, there is a slight positive association but otherwise there doesn't appear to be a strong relationship throughout the median area. See how when $Z > 1.2$, then $Y$ has an upwards trajectory whereas this is not true in general.
 
 ![Y plotted against Z!](/demo/dgp-association.png "Example of DGP associations.")

We can plot the true causal graph using igraph. 

 ```
g <- igraph::make_graph(~ Z--Y, Z--X, Y--W, X--W)
op <- par(mar=c(1,1,1,1))
igraph::plot.igraph(g,layout=igraph::layout_in_circle,vertex.size=30,vertex.color="lightblue")
 ```
  ![DGP Graph!](/demo/dgp.png "Data generating process.")
 
 
But note that based on the above plots that variables may not be related at all quantiles the same way! For example, comparing `Y` and `Z` above, we see that at the median of both variables we would not expect to see an edge between the two as their relationship appears noisy at best. However at the upper tail of Z there is a positive association between the two. 


## Running QGM and analyzing results
Let's run the QGM algorithm and see what the resulting graphs are. At the median quantile, $\tau = 0.5$, we have

 ![QGM at median!](/demo/qgm-median.png "QGM at the median.")

When we look at boundary quantiles, $\tau = 0.9$, we have

 ![QGM at upper tails!](/demo/qgm-upper.png "QGM at the upper tail.")

This is working as intended! Because of the relationship between `Y` and `Z`, we don't expect their edge to show up in the median graph but we do expect it at the extreme tails.
 
 
 
 
 
 
 
 
