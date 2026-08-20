






































# **tensortools**

Tensor operations, decompositions, and distributions in R.

``` r
library("tensortools")
library("tidyverse")
library("tictoc")
old_digits <- getOption("digits")
options(digits = 3)
```

# Operations

## n-mode prod

The $n$-mode matrix product of an $O$-th order tensor
$\mathcal{A} \in \mathbb{k}^{I_1 \times I_2 \times \dots \times I_{n} \times \dots \times I_{O}}$
with a matrix $\textbf{U} \in \mathbb{k}^{J \times I_n}$  
is

```math
(\mathcal{A} \times_n \textbf{U}) \in \mathbb{k}^{I_1 \times \dots \times I_{n-1} \times J \times I_{n+1} \times \dots \times I_{O}}
```

with entries

```math
(\mathcal{A} \times_n \textbf{U})_{i_1, \dots, i_{n-1}, j, i_{n+1}, \dots i_{O}} = \sum_{i_n = 1}^{I_{N}} a_{i_{1}, i_{2}, \dots, i_{O}} u_{j, i_n}.
```

The tensor and matrix share one mode in common, denoted here as $I_{n}$.
This operation is also called the tensor times matrix product.

The n-mode product between a tensor and a matrix can be computed with
`n_mode_prod()`.

``` r
a <- array(1:3, dim = c(3, 1, 1))
x <- matrix(4:9, nrow = 2, ncol = 3)

n_prod(a, 1, x)
#> , , 1
#> 
#>      [,1]
#> [1,]   40
#> [2,]   46
```

## nm-mode prod

The $nm$-mode product between a tensor and a tensor can be computed with
`nm_prod()`. The $nm$-mode matrix product of an $O$-th order tensor
$\mathcal{A} \in \mathbb{k}^{I_{1} \times I_{2} \times \dots \times I_{n} \times \dots \times I_{O}}$
with a $P$-th order tensor
$\textbf{U} \in \mathbb{k}^{J_{1} \times {J_2} \times \dots \times J_{m} = I_{n} \dots \times J_{P}}$
is

```math
(\mathcal{A} \sideset{_{n}}{_{m}}{\mathop{\boldsymbol{\times}}} \mathcal{U}) \in \mathbb{k}^{I_1 \times \dots \times I_{n-1} \times I_{n+1} \times \dots \times I_{O} \times J_{1} \times \dots \times J_{m-1} \times J_{m+1} \times \dots \times J_{P}}
```

with entries

```math
(\mathcal{A} \sideset{_{n}}{_{m}}{\mathop{\boldsymbol{\times}}} \mathcal{U})_{i_1, \dots, i_{n-1}, k, i_{n+1}, \dots I_{O}, j_{1}, \dots, j_{m-1}, j_{m+1}, \dots j_{P}} = \sum_{i_n = 1}^{I_n} a_{i_1, i_2, \dots, i_{n}, \dots i_N} u_{j_{1}, j_{2}, \dots i_n, \dots j_{P}}.
```

The $n$-th mode of the first tensor matches the $m$-th mode of the
second tensor, denoted here as $I_{n} = J_{m}$. This operation is also
called the tensor times tensor product.

``` r
A <- matrix(c(1, 2, 3, 4), nrow = 2)
x <- matrix(c(5, 6), nrow = 2)

nm_prod(A, x, 1, 1)
#>      [,1]
#> [1,]   17
#> [2,]   39
```

## [Tensor product](https://en.wikipedia.org/wiki/Tensor_product)

The tensor product between a tensor and a tensor can be computed with
`tensor_prod()`. Take 2 tensors
$\mathcal{A} \in \mathbb{k}^{I_{1} \times I_{2} \times \dots \times I_{O}}$
and
$\mathcal{U} \in \mathbb{k}^{J_{1} \times J_{2} \times \dots \times J_{P}}$

Then the tensor product is
```math
(\mathcal{A} \otimes \mathcal{B}) \in \mathbb{k}^{I_{1} \times I_{2} \dots \times I_{O} \times J_{1} \times \dots \times J_{P}}
```
with entries

```math
(\mathcal{A} \otimes \mathcal{U})_{i_1, \dots, i_{O}, j_{1}, \dots, j_{P}} = a_{i_1, i_2, \dots, i_{o}} u_{j_{1}, j_{2}, \dots \dots j_{p}}.
```

``` r
A <- matrix(c(1, 2, 3, 4), nrow = 2)
x <- matrix(c(5, 6), nrow = 2)

tensor_prod(A, x)
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    5   15
#> [2,]   10   20
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> [1,]    6   18
#> [2,]   12   24
```

# Simulating tensor variate normal draws

The density function of the multilinear normal distribution
is<sup>1</sup>

```math
f(x) =
(2\pi)^{-p^{*}/2}
\biggl(\prod_{i=1}^{k} |\Sigma_i|^{-p^{*}/(2p_i)}\biggr)
\exp\biggl\{ -\frac{1}{2} (x-\mu)^{T} \Sigma_{1:k}^{-1} (x-\mu) \biggr\}
```

where $\Sigma$ is positive definite, $x \in \mathbb{R}^p$,
$\mu \in \mathbb{R}^p$ and $\Sigma_{1:k} \in \mathbb{R}^p.$

The function `rtnorm()` simulates random draws from the tensor variate
normal with a specified mean array mu and a list of covariance matrices
called list_sigmas.

Since the univariate normal and matrix variate normal are simpler cases
of the tensor variate normal, this function can simulate from them as
well. Here, we simulate random draws from a univariate normal with mean
$-2$ and variance $4$.

``` r
univar_draws <- rtnorm(n = 1000, mu = -2, 
                       sigmas = list(matrix(4)))

mean(simplify2array(univar_draws))
#> [1] -2.02

var(simplify2array(univar_draws))
#> [1] 4.28
```

We can also simulate random draws from a multivariate normal
distribution by specifying a mean vector and the covariance matrix. The
dimensions of the draws are dependent on the $\mu$ provided and the
corresponding covariance matrices $\mathbf{\Sigma}$ must conform to the
$\mu$ provided.

``` r
S1 <- crossprod(matrix(data = c(1, 0.5, 0.5, 1), nrow = 2))

multivar_draws <- rtnorm(n = 1000, mu = c(2, 3), sigmas = list(S1))

multivar_draws[1:5]
#> [[1]]
#> [1] 3.268930 3.420696
#> 
#> [[2]]
#> [1] 3.243178 2.705054
#> 
#> [[3]]
#> [1] 1.026441 3.307681
#> 
#> [[4]]
#> [1] 2.235605 3.536821
#> 
#> [[5]]
#> [1] 2.077587 3.024604
```

And now we simulate from the matrix variate normal of size $2 \times 3$.
By default, the covariance matrices will be the identity.

``` r
matrix_draws <- rtnorm(n = 1e3, mu = matrix(1:6, nrow = 2, ncol = 3))
```

When simulating draws, the output is represented a list of $n$ draws.
Each element of the list represents one tensor observation. This makes
it clear that the sample size indexes observations and works well with
other functions in the package. Functions such as `simplify2array()` can
be used to create a higher-order tensor with one mode being the sample
index.

``` r
matrix_draws[[1]]
#>           [,1]     [,2]     [,3]
#> [1,] 1.7391149 1.483627 3.674582
#> [2,] 0.8653698 3.381173 6.263703
```

Below is a simulation from a tensor variate normal of size
$3 \times 4 \times
2$.The covariance structure is specified with a list of matrices.

``` r
mu_true <- array(1:24, dim = c(3, 4, 2))

S1 <- crossprod(matrix(rnorm(9), nrow = 3))
S2 <- crossprod(matrix(rnorm(16), nrow = 4))
S3 <- crossprod(matrix(rnorm(4), nrow = 2))

tvn_draws <- rtnorm(n = 1e3, mu = mu_true,
                    sigmas = list(S1, S2, S3))

tvn_draws[[1]]
#> , , 1
#> 
#>           [,1]      [,2]     [,3]      [,4]
#> [1,] 0.1772103  2.686398  3.74531  8.142328
#> [2,] 1.6455212 10.565258 10.34400 16.043991
#> [3,] 2.1656371  5.491703 11.64316 13.430647
#> 
#> , , 2
#> 
#>          [,1]     [,2]     [,3]     [,4]
#> [1,] 13.96192 16.68971 19.00767 22.71868
#> [2,] 14.12012 15.10580 20.49584 21.57830
#> [3,] 14.41665 20.40647 18.01518 24.32685
```

# Estimation

## MLE estimation for the tensor variate normal

The package supports MLE estimation. Given an array of draws, it will
return a list containing the MLE for the mean and covariance matrices.
Earlier, we generated univariate normal draws with $\mu = 2$ and
$\sigma^{2} = 4.$ Calling the function `tensor_mle()` will return MLEs
similar to the true results.

``` r
(univarnorm_est <- tensor_mle(draws = univar_draws, model = "normal"))
#> $mu
#> [1] -2.023296
#> 
#> $sigmas
#> $sigmas[[1]]
#> [1] 4.284203
```

The multivariate normal draws had a mean vector of
$\begin{pmatrix} 2 \\\\ 3 \end{pmatrix}$ and a covariance matrix
$\begin{pmatrix} 1.25 & 1 \\\\ 1 & 1.25 \end{pmatrix}.$

``` r
(multivarnorm_est <- tensor_mle(draws = multivar_draws, model = "normal"))
#> $mu
#> [1] 1.981819 2.995725
#> 
#> $sigmas
#> $sigmas[[1]]
#>          [,1]     [,2]
#> [1,] 1.350600 1.098575
#> [2,] 1.098575 1.371289
```

For the matrix variate draws, the mean matrix was
$\begin{pmatrix} 1 & 3 & 5 \\\\ 2 & 4 & 6 \end{pmatrix}$ and the
covariance matrices were the identity.

``` r
matrix_est <- tensor_mle(matrix_draws, model = "normal")

matrix_est$mu
#>          [,1]     [,2]     [,3]
#> [1,] 1.016722 2.988299 4.995788
#> [2,] 1.979936 3.979311 5.978654
matrix_est$sigmas |> lapply(round, 3)
#> [[1]]
#>       [,1]  [,2]
#> [1,] 0.996 0.007
#> [2,] 0.007 1.004
#> 
#> [[2]]
#>        [,1]   [,2]   [,3]
#> [1,]  1.028 -0.042  0.007
#> [2,] -0.042  1.018 -0.010
#> [3,]  0.007 -0.010  0.997
```

The `tensor_mle()` function also works for tensor-valued data. Below, we
compare the estimated and true covariance matrices using the mean
squared error (MSE).

Also, note that there is nonidentifiability in the scaling of the
covariance matrices. Given the tensor variate normal draws, we know the
kronecker of all of the covariances, but not the scaling of each one. To
compare the true and estimated covariance matrices, we first trace
normalization. This removes arbitrary scale differences while preserving
the relative covariance structure within each mode.

``` r
tensor_est <- tensor_mle(draws = tvn_draws, model = "normal")

true_sigmas <- list(S1, S2, S3)

for(i in 1:3) {
  true_scaled <- true_sigmas[[i]] / sum(diag(true_sigmas[[i]]))
  est_scaled <- tensor_est$sigmas[[i]] / sum(diag(tensor_est$sigmas[[i]]))
  
  print(mean((est_scaled - true_scaled)^2))
}
#> [1] 5.504582e-05
#> [1] 1.754755e-05
#> [1] 2.100553e-05
```

# Other models

To generate tensor variate skewed distributions, we will use the normal
variance mean mixture model An $r$-dimensional random vector $x$ is a
normal variance-mean mixture with mixing distribution $F$ if, for a
given $u \ge 0$ that follows a probability distribution $F$ on
$[0, \infty)$, $x|u \sim N_r(\mu + u \beta, u \Sigma)$.
$\Sigma \mathbb{R}^{r \times r}$ is a constant, positive-definite matrix
and $\mu \in \mathbb{R}^{r}$ and $\beta \in \mathbb{R}^{r}$ are constant
vectors.

We state that $x$ is a normal variance-mean mixture with position $\mu$,
drift $\beta$, structure matrix $\Sigma$, and mixing distribution $F.$
Different mixing distributions lead to different distributions. For
example, if $F$, the mixing distribution, is the generalized inverse
Gaussian distribution GIG($\lambda, \delta^2, \kappa^2$) then the
resulting distribution of $x$ is a generalized hyperbolic distribution.

In the tensor variate case, the normal variance mean mixture model is an
efficient way to introduce skewness. A random tensor $\mathcal{X}$ can
be written as
$\mathcal{X} = \mathcal{M} + \textbf{W} \mathcal{A} + \sqrt{\textbf{W}} \mathcal{V}.$

## Tensor variate normal inverse Gaussian

For the tensor variate normal inverse Gaussian distribution, we let
$W \sim \mathrm{IG}(1, \kappa)$.

The density function of the normal inverse Gaussian distribution is

```math
\begin{aligned}
f_{\text{TVNIG}}(\mathcal{X}|\mathbf{V})
&= \frac{2 \exp\biggl\{\text{vec}(\mathcal{X} - \mathcal{M})^{T}
\bigotimes_{d=1}^{D} \Sigma_{d}^{-1} \text{vec}(\mathcal{A}) + \kappa \biggr\}}
{(2\pi)^{\frac{n^{*}+1}{2}} \prod_{d=1}^{D} |\Sigma_{d}|^{\frac{n^{*}}{2n_{d}}}} \\
&\quad \times
\biggl(\frac{\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + 1}
{\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \kappa^2} \biggr)^{-\frac{1 + n^{*}}{4}} \\
&\quad \times
K_{-\frac{1 + n^{*}}{2}} \biggl(\sqrt{\Bigl[\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \kappa^{2}\Bigr]
\Bigl[\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + 1\Bigr]} \biggr).
\end{aligned}
```

where $\Sigma$ is positive definite, $x \in \mathbb{R}^p$,
$\mu \in \mathbb{R}^p$ and $\Sigma_{1:k} \in \mathbb{R}^p.$

The function `rtinvgauss()` simulates random draws from the tensor
variate normal inverse Gaussian with a specified mean array mu, a list
of covariance matrices, a skew array, and $\kappa$ which describes the
shape of the inverse Gaussian distribution

``` r
mu_true <- array(1:24, dim = c(3, 4, 2))
skew_true <- array(seq(0, 4, length.out = 24), dim = c(3, 4, 2))
kappa_true <- 2

invgauss_draws <- rtinvgauss(n = 1e3, mu = mu_true, skew = skew_true, 
                             sigmas = list(S1, S2, S3), kappa = kappa_true)
```

``` r
invgauss_est <- tensor_mle(invgauss_draws, model = "invgauss")
```

For the inverse gaussian distribution, the true mean
$E[X] = \mathcal{M} + \mathcal{A}/\kappa$. We can compare the model’s
estimation for the mean with the true values using MSE.

``` r
mean((with(invgauss_est, mu + skew/kappa) -
        (mu_true + skew_true/kappa_true))^2)
#> [1] 0.002919574
```

``` r
mean((invgauss_est$mu - mu_true)^2)
#> [1] 0.01684996
mean((invgauss_est$skew - skew_true)^2)
#> [1] 2.726522
invgauss_est$kappa
#> [1] 0.6375593

for(i in 1:2) {
  true_scaled <- true_sigmas[[i]] / sum(diag(true_sigmas[[i]]))
  est_scaled <- invgauss_est$sigmas[[i]] / sum(diag(invgauss_est$sigmas[[i]]))
  
  print(mean((est_scaled - true_scaled)^2))
}
#> [1] 1.194769e-05
#> [1] 1.17049e-05
```

## Tensor variate generalized hyperbolic

For the tensor variate generalized hyperbolic distribution, we let
$W \sim I(\omega, 1, \lambda)$.

The density function of the generalized hyperbolic distribution is

```math
\begin{aligned}
f_{\text{TVGH}}(\mathcal{X}|\mathbf{V})
&= \frac{\exp\biggl\{\text{vec}(\mathcal{X} - \mathcal{M})^{T}
\bigotimes_{d=1}^{D} \Sigma_{d}^{-1} \text{vec}(\mathcal{A})\biggr\}}
{(2\pi)^{\frac{n^{*}}{2}} \prod_{d=1}^{D} |\Sigma_{d}|^{\frac{n^{*}}{2n_{d}}} K_{\lambda}(\omega)} \\
&\quad \times
\biggl(\frac{\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \omega}
{\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \omega} \biggr)^{\frac{\lambda - n^{*}/2}{2}} \\
&\quad \times
K_{\lambda - n^{*}/2} \biggl(\sqrt{\Bigl[\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \omega\Bigr]
\Bigl[\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \omega\Bigr]}\biggr).
\end{aligned}
```
where $\Sigma$ is positive definite, $x \in \mathbb{R}^p$,
$\mu \in \mathbb{R}^p$ and $\Sigma_{1:k} \in \mathbb{R}^p.$

The function `rtgenhyper()` simulates random draws from the tensor
variate generalized hyperbolic with a specified mean array mu and a list
of covariance matrices called list_sigmas.

Since the univariate normal and matrix variate normal are simpler cases
of the tensor variate normal, this function can simulate from them as
well. Here, we simulate random draws from a univariate normal with mean
$-2$ and variance $4$.

``` r
lambda_true <- 10
omega_true <- 12
 
genhyper_draws <- rtgenhyper(n = 1e3, mu = mu_true, skew = skew_true,
                             sigmas = list(S1, S2, S3),
                             lambda = lambda_true, omega = omega_true)
```

For the generalized hyperbolic distribution, the true mean
$E[X] = \mathcal{M} + \frac{K_{\lambda + 1}(\omega)}{K_{\lambda}(\omega)} \mathcal{A}$.
We can compare the model’s estimation for the mean with the true values
using MSE.

``` r
genhyper_est <- tensor_mle(genhyper_draws, model = "genhyper", quiet = FALSE, tol = 1e-5)
#> Iteration 50: criterion = 1.962e-05
#> Converged at iteration 74

mean((
  with(genhyper_est, mu + besselK(x = omega, nu = lambda + 1)/
                    besselK(x = omega, nu = lambda) * skew) -
    (mu_true + besselK(x = omega_true, nu = lambda_true + 1)/
       besselK(x = omega_true, nu = lambda_true) * skew_true)
  )^2)
#> [1] 0.004171649

mean((genhyper_est$mu - mu_true)^2)
#> [1] 0.8171771
mean((genhyper_est$skew - skew_true)^2)
#> [1] 3.339013

for(i in 1:3) {
  true_scaled <- true_sigmas[[i]] / sum(diag(true_sigmas[[i]]))
  est_scaled <- genhyper_est$sigmas[[i]] / sum(diag(genhyper_est$sigmas[[i]]))
  
  print(mean((est_scaled - true_scaled)^2))
}
#> [1] 2.80842e-05
#> [1] 1.591615e-05
#> [1] 7.089782e-06

genhyper_est$lambda
#> [1] 8.624151
genhyper_est$omega
#> [1] 2.14075
```

## Tensor variate variance gamma

For the tensor variate variance gamma distribution, we let
$W \sim \text{Gamma}(\gamma, \gamma)$.

The density function of the variance gamma distribution is

```math
\begin{aligned}
f_{\text{TVVG}}(\mathcal{X}|\mathbf{V})
&= \frac{2\gamma^{\gamma} \exp\biggl\{\text{vec}(\mathcal{X} - \mathcal{M})^{T}
\bigotimes_{d=1}^{D} \Sigma_{d}^{-1} \text{vec}(\mathcal{A})\biggr\}}
{(2\pi)^{\frac{n^{*}}{2}} \prod_{d=1}^{D} |\Sigma_{d}|^{\frac{n^{*}}{2n_{d}}} \Gamma(\gamma)} \\
&\quad \times
\biggl(\frac{\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1})}
{\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + 2\gamma} \biggr)^{\frac{\gamma - n^{*}/2}{2}} \\
&\quad \times
K_{\gamma - n^{*}/2} \biggl(\sqrt{\Bigl[\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + 2\gamma\Bigr]
\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1})} \biggr).
\end{aligned}
```

``` r
vargamma_draws <- rtvargamma(n = 1e3, mu = mu_true, skew = skew_true,
                             sigmas = list(S1, S2, S3))
```

``` r
vargamma_est <- tensor_mle(vargamma_draws, model = "vargamma", quiet = FALSE)
#> Converged at iteration 6

mean((with(vargamma_est, mu + skew) - (mu_true + skew_true))^2)
#> [1] 0.1187094

mean((vargamma_est$mu - mu_true)^2)
#> [1] 0.003826563
mean((vargamma_est$skew - skew_true)^2)
#> [1] 0.1008045

for(i in 1:3) {
  true_scaled <- true_sigmas[[i]] / sum(diag(true_sigmas[[i]]))
  est_scaled <- vargamma_est$sigmas[[i]] / sum(diag(vargamma_est$sigmas[[i]]))
  
  print(mean((est_scaled - true_scaled)^2))
}
#> [1] 1.333307e-05
#> [1] 8.393338e-06
#> [1] 3.607587e-06

vargamma_est$gamma
#> [1] 1.830932
```

## Tensor variate skewed t

For the tensor variate skewed t distribution, we let
$W \sim \mathrm{Inv\text{-}Gamma}(\nu/2, \nu/2)$.

The density function of the skewed t distribution is

```math
\begin{aligned}
f_{\text{TVST}}(\mathcal{X}|\mathbf{V})
&= \frac{2\bigl(\frac{\nu}{2}\bigr)^{\nu/2} \exp\biggl\{\text{vec}(\mathcal{X} - \mathcal{M})^{T}
\bigotimes_{d=1}^{D} \Sigma_{d}^{-1} \text{vec}(\mathcal{A})\biggr\}}
{(2\pi)^{\frac{n^{*}}{2}} \prod_{d=1}^{D} |\Sigma_{d}|^{\frac{n^{*}}{2n_{d}}} \Gamma\bigl(\frac{\nu}{2}\bigr)} \\
&\quad \times
\biggl(\frac{\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \nu}
{\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1})} \biggr)^{-\frac{\nu + n^{*}}{4}} \\
&\quad \times
K_{-\frac{\nu + n^{*}}{2}} \biggl(\sqrt{\rho(\mathcal{A}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1})
\Bigl[\delta(\mathcal{X}; \mathcal{M}, \bigotimes_{d=1}^{D} \Sigma_{d}^{-1}) + \nu\Bigr]} \biggr).
\end{aligned}
```

``` r
nu_true <- 20
skew_true <- array(c(
   0.6,  0.3, -0.2,  0.5,
   0.4, -0.5,  0.2, -0.3,
  -0.4,  0.7,  0.1, -0.6,

   0.5, -0.2,  0.3,  0.4,
  -0.3,  0.6, -0.4,  0.2,
   0.1, -0.5,  0.7, -0.1
), dim = c(3, 4, 2))

skewt_draws <- rtskewt(n = 1e3, mu = mu_true, skew = skew_true,
                       sigmas = list(S1, S2, S3), nu = nu_true)
```

``` r
skewt_est <- tensor_mle(skewt_draws, model = "skewt",
                        quiet = FALSE)
#> Converged at iteration 28

mean((with(skewt_est, mu + (nu)/(nu-2) * skew) -
        (mu_true + nu_true / (nu_true - 2) * skew_true))^2)
#> [1] 0.002083219

for(i in 1:3) {
  true_scaled <- true_sigmas[[i]] / sum(diag(true_sigmas[[i]]))
  est_scaled <- skewt_est$sigmas[[i]] / sum(diag(skewt_est$sigmas[[i]]))
  
  print(mean((est_scaled - true_scaled)^2))
}
#> [1] 6.403674e-06
#> [1] 1.334612e-06
#> [1] 1.931086e-06
```

## Assessing tensor variate normality

We assess tensor variate normality using the Mahalanobis distance. If
$\mathcal{X}$ follows an order-$O$ tensor variate normal distribution,
then $\text{vec}(\mathcal{X})$ follows a multivariate normal
distribution with mean $\text{vec}(\mathcal{M})$ and covariance matrix
$\boldsymbol{\Sigma}_{O} \otimes \cdots \otimes \boldsymbol{\Sigma}_{1}$.
This is why we can vectorize the tensor draws and compare the
multivariate Mahalanobis distances to the tensor variate Mahalanobis
distances. If we have properly reconstructed the structure using the
MLEs $\hat{\mathcal{M}}$ and
$\hat{\Sigma}_{O}, \cdots \hat{\Sigma}_{1}$, then the standardized draws
should be from a standard normal.

If a tensor normal covariance structure is present, then the
multivariate and tensor variate distances should be close to one
another. Following the matrix variate normality paper, this can be
visualized with a distance-distance plot: when the tensor normal
structure is appropriate, the distances lie close to the reference line;
when the Kronecker product structure is absent, the distances diverge.

The function `mahalanobis_dist()` computes both distances using fitted
MLEs, and `mahalanobis_test()` performs a two-sample Kolmogorov-Smirnov
test to compare their empirical distributions.

``` r
A <- matrix(rnorm(6^2), 6, 6)
s_vec <- crossprod(A)   # positive definite, general covariance

vec_draws <- rtnorm(n = 1e3, mu = 1:6, sigmas = list(s_vec))

matrix_nocovar_draws <- lapply(1:1e3, function(i) {
  array(vec_draws[[i]], dim = c(2, 3))
})
```

``` r
plot_malanobis <- function(draws, title = "Dist") {
  distances <- mahalanobis_dist(draws)
  test <- mahalanobis_test(distances)
  plot_dist(distances)
  title(main = title)
  with(test, mtext(text = sprintf("D = %.3f, p = %.3f", statistic, p.value)), 
       side = 3, line = 1)
}
```

``` r
par(mfcol = c(1, 4))
plot_malanobis(multivar_draws, title = "Multivar")
#> Warning in ks.test.default(vec, tensor): p-value will be approximate in the
#> presence of ties
plot_malanobis(matrix_draws, title = "Mat norm")
plot_malanobis(matrix_nocovar_draws, title = "Mat no covar norm")
plot_malanobis(tvn_draws, title = "Tensor norm")
```

<img src="man/figures/README-unnamed-chunk-27-1.png" width="100%" />

We are also interested in performing a likelihood ratio test. We are
curious if one of the covariance matrices might be the identity matrix.
Thus, we run the algorithm and restrict these matrices.

For a restriction set $\mathcal{R} \subset \{1,\ldots,D\}$, the
restricted model fixes $\Sigma_r = I_{n_r}$ for every
$r \in \mathcal{R}$, while all other covariance matrices are estimated.
The likelihood ratio test compares

```math
H_0: \Sigma_r = I_{n_r} \ \forall r \in \mathcal{R} \quad \text{vs.} \quad H_A: \Sigma_r \neq I_{n_r} \ \text{for at least one } r \in \mathcal{R}.
```

The test statistic is

```math
2\{\ell(\hat\theta_{\mathrm{full}}) - \ell(\hat\theta_{\mathcal{R}})\}
\sim \chi^2_{\nu_{\mathcal{R}}}, \quad \nu_{\mathcal{R}} = \sum_{r \in \mathcal{R}} \biggl\{\frac{n_r(n_r+1)}{2} - 1\biggr\}
```

If the $p$-value is large, then we fail to reject the null hypothesis
and conclude the restricted model fits the data as well as the
unrestricted model. If the $p$-value is small, we reject the null
hypothesis, which states the unrestricted model has a significantly
better fit than the restricted model.

``` r
true_str <- kronecker(S3, S2) |> kronecker(S1)
  
test_restrict <- function(draws, restrict) {
  curr_restrict <- tensor_mle(draws, restrict = restrict, model = "normal")
  
  if(is.null(restrict)) restrict <- "None"

  curr_str <- with(curr_restrict, kronecker(sigmas[[3]], sigmas[[2]]) |>
                                  kronecker(sigmas[[1]]))

  tibble(mse = mean((curr_str - true_str)^2),
         loglik = curr_restrict$loglik,
         k = curr_restrict$k,
         BIC = curr_restrict$BIC,
         restrict = paste(restrict, collapse = ", "))
}

(res <- map_dfr(list(NULL, c(1), c(2), c(3), c(1, 2), c(1, 3), c(2, 3)),
        draws = tvn_draws, test_restrict))
#> # A tibble: 7 × 5
#>       mse  loglik     k     BIC restrict
#>     <dbl>   <dbl> <dbl>   <dbl> <chr>   
#> 1 0.00724 -42243.    41  84769. None    
#> 2 1.80    -51553.    36 103354. 1       
#> 3 2.47    -56939.    32 114099. 2       
#> 4 3.91    -60844.    39 121958. 3       
#> 5 3.39    -65847.    27 131881. 1, 2    
#> 6 4.34    -69754.    34 139742. 1, 3    
#> 7 4.49    -75248.    30 150703. 2, 3
```

``` r
lrt_test <- function(draws, restrict) {
  full_model <- tensor_mle(draws, model = "normal")
  restrict_model <- tensor_mle(draws, model = "normal", restrict = restrict)
  
  pval <- pchisq(2 * (full_model$loglik - restrict_model$loglik), 
                 df = full_model$k - restrict_model$k, lower.tail = FALSE)
  
  tibble(restrict = paste(restrict, collapse = ", "), 
         pval = pval)
}

map_dfr(list(NULL, c(1), c(2), c(3), c(1, 2), c(1, 3), c(2, 3)),
        draws = tvn_draws, lrt_test)
#> # A tibble: 7 × 2
#>   restrict  pval
#>   <chr>    <dbl>
#> 1 ""           1
#> 2 "1"          0
#> 3 "2"          0
#> 4 "3"          0
#> 5 "1, 2"       0
#> 6 "1, 3"       0
#> 7 "2, 3"       0
```

``` r
true_str <- kronecker(S3, diag(4)) |> kronecker(S1)

second_identity <- rtnorm(n = 1e3, mu = mu_true, 
                          sigmas = list(S1, diag(4), S3))

map_dfr(list(NULL, c(1), c(2), c(3), c(1, 2), c(1, 3), c(2, 3)), 
        draws = second_identity, test_restrict)
#> # A tibble: 7 × 5
#>         mse  loglik     k    BIC restrict
#>       <dbl>   <dbl> <dbl>  <dbl> <chr>   
#> 1 0.000206  -25784.    41 51852. None    
#> 2 0.0673    -34465.    36 69180. 1       
#> 3 0.0000499 -25794.    32 51809. 2       
#> 4 0.0742    -30418.    39 61105. 3       
#> 5 0.0672    -34480.    27 69146. 1, 2    
#> 6 0.115     -39055.    34 78345. 1, 3    
#> 7 0.0741    -30435.    30 61077. 2, 3

map_dfr(list(NULL, c(1), c(2), c(3), c(1, 2), c(1, 3), c(2, 3)),
        draws = second_identity, lrt_test)
#> # A tibble: 7 × 2
#>   restrict   pval
#>   <chr>     <dbl>
#> 1 ""       1     
#> 2 "1"      0     
#> 3 "2"      0.0228
#> 4 "3"      0     
#> 5 "1, 2"   0     
#> 6 "1, 3"   0     
#> 7 "2, 3"   0
```

## Tensor Variate Skewed T

``` r
tibble(skew25 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                               skew = matrix(0.25, nrow = 2, ncol = 3), nu = 6)),
       skew50 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                               skew = matrix(0.5, nrow = 2, ncol = 3), nu = 6)),
       skew0 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                              skew = matrix(0, nrow = 2, ncol = 3), nu = 6)),
       skew1 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                              skew = matrix(1, nrow = 2, ncol = 3), nu = 6)),
       skew5 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                              skew = matrix(5, nrow = 2, ncol = 3), nu = 6)),
       skew10 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                               skew = matrix(10, nrow = 2, ncol = 3), nu = 6))) |>
  pivot_longer(cols = everything()) |> 
  mutate(name = factor(name, levels = c("skew0", "skew1", 
                                     "skew5", "skew10", "skew25", "skew50"))) |> 
  ggplot() +
  geom_histogram(aes(x = value), color = "black", fill = "pink") +
  facet_wrap(~name, nrow = 2, scales = "free_x") +
  theme_minimal()
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

<img src="man/figures/README-unnamed-chunk-31-1.png" width="100%" />

``` r

tibble(nu05 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                             skew = matrix(1, nrow = 2, ncol = 3), nu = 0.5)),
       nu1 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                            skew = matrix(1, nrow = 2, ncol = 3), nu = 1)),
       nu5 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                            skew = matrix(1, nrow = 2, ncol = 3), nu = 5)),
       nu10 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                             skew = matrix(1, nrow = 2, ncol = 3), nu = 10)),
       nu20 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                             skew = matrix(1, nrow = 2, ncol = 3), nu = 20)),
       nu50 = unlist(rtskewt(n = 1e3, mu = matrix(0, nrow = 2, ncol = 3), 
                             skew = matrix(1, nrow = 2, ncol = 3), nu = 50))) |>
  pivot_longer(cols = everything()) |>
  mutate(name = factor(name, levels = c("nu05", "nu1", "nu5", "nu10", "nu20", "nu50"))) |> 
  ggplot() +
  geom_histogram(aes(x = value), color = "black", fill = "pink") +
  facet_wrap(~name, nrow = 2, scales = "free_x") +
  theme_minimal()
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

<img src="man/figures/README-unnamed-chunk-31-2.png" width="100%" />

Suppose $X_{n \times p} \sim N_{np}(0, \Sigma_1, \Sigma_2)$. A natural
question is whether the usual sample covariance of the rows can estimate
$\Sigma_2$, and whether the sample covariance of the columns can
estimate $\Sigma_1$.

If we compute the usual row covariance

```math
S_{\text{row}} = \frac{1}{n-1} X^T \biggl(I_n - \frac{1}{n}\mathbf{1}_n \mathbf{1}_n^T\biggr) X,
```

then under the matrix normal model,

```math
\mathbb{E}(S_{\text{row}})
= \frac{\mathrm{tr}\left[\left(I_n - \frac{1}{n}\mathbf{1}_n \mathbf{1}_n^T\right)\Sigma_1\right]}{n-1}\Sigma_2.
```

So the usual sample covariance does recover the of $\Sigma_2$, but only
up to a multiplicative constant determined by $\Sigma_1$. The same
argument with $X^T$ shows that the usual column covariance recovers the
shape of $\Sigma_1$, but only up to a multiplicative constant determined
by $\Sigma_2$.

The simulation below shows exactly that. The trace-normalized naive
covariance matrices are close to the truth in all three cases, but the
raw traces are distorted by the covariance matrix on the other mode.

``` r
mse <- function(mat1, mat2) {
  mat1_scale <- mat1/diag(mat1)
  mat2_scale <- mat2/diag(mat2)
  
  mean((mat1_scale - mat2_scale)^2)
}

compare_naive_covs <- function(sigmas, n_rep = 200,
                               n_row = 30, n_col = 4) {
  
  draws <- rtnorm(n = n_rep, mu = matrix(0, nrow = n_row, ncol = n_col),
                  sigmas = sigmas)

  draws_mat <- matrix(aperm(simplify2array(draws), c(1, 3, 2)), ncol = 4)
  transpose_draws_mat <- matrix(aperm(simplify2array(draws), c(2, 1, 3)), ncol = 30)
  
  s2_naive <- var(draws_mat)
  s1_naive <- var(transpose_draws_mat)
  
  mle_est <- tensor_mle(draws, model = "normal")
  
  s2_hat <- mle_est$sigmas[[2]]
  s1_hat <- mle_est$sigmas[[1]]

  tibble(name = c("s2_naive", "s1_naive", "s2_hat", "s1_hat"),
         error = c(mse(s2_naive, sigmas[[2]]),
                   mse(s1_naive, sigmas[[1]]),
                   mse(s2_hat, sigmas[[2]]),
                   mse(s1_hat, sigmas[[1]])),
         trace_ratio = c(sum(diag(s2_naive)) / sum(diag(sigmas[[2]])),
                         sum(diag(s1_naive)) / sum(diag(sigmas[[1]])),
                         sum(diag(s2_hat)) / sum(diag(sigmas[[2]])),
                         sum(diag(s1_hat)) / sum(diag(sigmas[[1]]))))
}

set.seed(123)

s1 <- crossprod(matrix(rnorm(30^2), nrow = 30))
s2 <- crossprod(matrix(rnorm(4^2), nrow = 4))

rbind(
  compare_naive_covs(list(diag(30), s2)) |> mutate(setting = "Sigma1 = I"),
  compare_naive_covs(list(s1, diag(4))) |> mutate(setting = "Sigma2 = I"),
  compare_naive_covs(list(s1, s2)) |> mutate(setting = "both non-identity")
) |> View()
```

In this table, `s2_error` measures how close `var(X)` is to the shape of
$\Sigma_2$ after trace normalization, and `s1_error` does the same for
`var(t(X))` and $\Sigma_1$. The trace ratios show the scale distortion.
For example, when $\Sigma_2 = I$, the usual covariance of the rows still
recovers the shape of $\Sigma_2$ well, but its trace is inflated by a
factor of about $29$, coming from $\Sigma_1$.

- `var(X)` can estimate the shape of $\Sigma_2$, but not its absolute
  scale, unless the scale contribution from $\Sigma_1$ is fixed.
- `var(t(X))` can estimate the shape of $\Sigma_1$, but not its absolute
  scale, unless the scale contribution from $\Sigma_2$ is fixed.
- This is exactly why the tensor and matrix normal models need an
  identifying constraint such as fixing the trace of one mode
  covariance.

``` r
make_skew <- function(dims, strength = 1) {
  A <- array(rnorm(prod(dims)), dim = dims)
  A / sqrt(sum(A^2)) * strength
}
```

``` r
model_compare <- function(iter) {
  s1 <- crossprod(matrix(rnorm(4), nrow = 2))
  s2 <- crossprod(matrix(rnorm(9), nrow = 3))
  s3 <- crossprod(matrix(rnorm(16), nrow = 4))

  skew_true <- make_skew(dims = c(2, 3, 4), strength = 1.5)
  
  draws <- rtskewt(n = 1e3, mu = array(1:24, dim = c(2, 3, 4)), 
                   sigmas = list(s1, s2, s3), skew = skew_true,
                   nu = 40)
  
  res_normal <- tensor_mle(draws, model = "normal", quiet = FALSE)
  res_skewt <- tensor_mle(draws, model = "skewt", quiet = FALSE)
  res_vargam <- tensor_mle(draws, model = "vargamma", quiet = FALSE)
  res_invgauss <- tensor_mle(draws, model = "invgauss", quiet = FALSE)
  res_genhyper <- tensor_mle(draws, model = "genhyper", quiet = FALSE)

  all_mod <- 
    tibble(model = c("Normal", "Skewt", "Vargamma", "Invgauss", "Genhyper"),
    loglik = c(res_normal$loglik, res_skewt$loglik,
               res_vargam$loglik, res_invgauss$loglik, res_genhyper$loglik),
    k = c(res_normal$k, res_skewt$k,
          res_vargam$k, res_invgauss$k, res_genhyper$k),
    AIC = c(res_normal$AIC, res_skewt$AIC,
            res_vargam$AIC, res_invgauss$AIC, res_genhyper$AIC),
    BIC = c(res_normal$BIC, res_skewt$BIC,
            res_vargam$BIC, res_invgauss$BIC, res_genhyper$BIC)
  )
  
  all_mod
  # return best model by BIC
  # tibble(iter = iter, 
  #        res = all_mod |> arrange(BIC) |> slice_head(n=1) |> pull(model))
}
```

## Naive classification

A natural supervised setup is to keep the tensor-valued observations and
their class labels as parallel objects:

``` r
train <- list(
  x = imgs,      # list of tensor-valued draws
  y = labels     # vector of class labels, same length as x
)
```

The labels do not need to be attached to each tensor. Keeping `x` and `y`
separate makes it easy to pass class-specific subsets directly to
`tensor_mle()`.

For a simple generative classifier, fit one class-conditional tensor
distribution per class. The distribution family can differ by class, as
long as all classes use the same tensor dimensions and preprocessing.
The prediction rule is the class with the largest posterior score,
`log(prior) + log(density)`.

``` r
candidate_models <- c("normal", "skewt", "vargamma", "invgauss", "genhyper")

fit_class_model <- function(x, model) {
  fit <- tensor_mle(x, model = model, quiet = TRUE)
  fit$model <- model
  fit
}

fit_naive_tensor_classifier <- function(x, y,
                                        models = candidate_models) {
  classes <- sort(unique(y))

  fits <- lapply(classes, function(cls) {
    class_x <- x[y == cls]

    model_fits <- lapply(models, function(model) {
      tryCatch(fit_class_model(class_x, model), error = function(e) NULL)
    })
    names(model_fits) <- models
    model_fits <- Filter(Negate(is.null), model_fits)

    best_model <- names(which.min(sapply(model_fits, `[[`, "BIC")))

    list(
      fit = model_fits[[best_model]],
      prior = length(class_x) / length(x)
    )
  })

  names(fits) <- classes
  fits
}

tensor_log_density <- function(x, fit) {
  switch(
    fit$model,
    normal = dtnorm(x, mu = fit$mu, sigmas = fit$sigmas, log = TRUE),
    skewt = dtskewt(x, mu = fit$mu, skew = fit$skew,
                    sigmas = fit$sigmas, nu = fit$nu, log = TRUE),
    vargamma = dtvargamma(x, mu = fit$mu, skew = fit$skew,
                          sigmas = fit$sigmas, scale = fit$gamma,
                          log = TRUE),
    invgauss = dtinvgauss(x, mu = fit$mu, skew = fit$skew,
                          sigmas = fit$sigmas, kappa = fit$kappa,
                          log = TRUE),
    genhyper = dtgenhyper(x, mu = fit$mu, skew = fit$skew,
                          sigmas = fit$sigmas, lambda = fit$lambda,
                          omega = fit$omega, log = TRUE)
  )
}

predict_naive_tensor_classifier <- function(x_new, classifier) {
  scores <- sapply(classifier, function(class_fit) {
    log(class_fit$prior) + tensor_log_density(x_new, class_fit$fit)
  })

  names(which.max(scores))
}
```

For example, if the normal model has the lowest BIC for cats and the
generalized hyperbolic model has the lowest BIC for maple trees, those
two class densities can still be compared on a new image because both
return log densities on the same input space.

## Installation

You can install the development version of tensortools from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("ankan12/tensortools")
```

Reset digits

``` r
options(digits = old_digits)
```

<div id="refs" class="references csl-bib-body" entry-spacing="0"
line-spacing="2">

<div id="ref-tensornormprop" class="csl-entry">

<span class="csl-left-margin">1.
</span><span class="csl-right-inline"><span class="nocase">Ohlson, M.,
Rauf Ahmad, M. & von Rosen, D.</span> [The multilinear normal
distribution: Introduction and some basic
properties](https://doi.org/10.1016/j.jmva.2011.05.015). *Journal of
Multivariate Analysis* **113**, 37–47 (2013).</span>

</div>

</div>
