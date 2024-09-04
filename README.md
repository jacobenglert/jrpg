
<!-- README.md is generated from README.Rmd. Please edit that file -->

# jrpg

<!-- badges: start -->
<!-- badges: end -->

`jrpg` is an R package that is used to facilitate the sampling of latent
Pólya-Gamma random variables in Bayesian models with logistic
likelihoods. This data augmentation strategy, originally proposed in
[Polson et al.
(2013)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001)
and [Windle
(2013)](https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1)
is very useful for writing custom MCMC algorithms (see also this
[post](https://gregorygundersen.com/blog/2019/09/20/polya-gamma/) by
Gregory Gunderson and this
[post](https://tiao.io/post/polya-gamma-bayesian-logistic-regression/)
by Louis Tiao for a blog-style overview).

This package seeks to compile existing implementations of $PG(b,z)$
samplers into one hybrid sampler. Recently, a Julia language
[package](https://github.com/wzhorton/PolyaGammaHybridSamplers.jl/tree/main)
has been created to address the same problem.

# Previous Implementations

The Pólya-Gamma packages I have found to be most useful in my work thus
far are:

- [BayesLogit](https://github.com/jwindle/BayesLogit/tree/master) (TMSA
  Lab @ UIUC)
- [pg](https://github.com/tmsalab/pg/tree/main) (Jesse Windle)
- [pgdraw](https://cran.r-project.org/web//packages/pgdraw/index.html)
  (Enes Makalic and Daniel F Schmidt)

There are no doubt other implementations out there, but these are the
most popular in my experience. All of the above implement some form of
hybrid $PG(b,z)$ sampler, which automatically selects the fastest
sampling strategy for the supplied parameters. For example, `BayesLogit`
has the `rpg` function and `pg` has the `rpg_hybrid` function. However,
the hybrid samplers differ between packages. This leads to certain
implementations being “better” for different scenarios.

In general, there are 4 reliable sampling techniques for $PG(b,z)$
random variables used throughout these packages:

1.  Sum-of-gammas approximation: best for very small values of $b$. If
    $X \sim PG(b,z)$, then the density of $X$ can be expressed as an
    infinite sum of independent $\text{Gamma}(b,1)$ random variables.
    This method simply truncates the infinite sum to, say, 200 terms.
2.  Devroye Method: best for integer values of $b$. Samples $b$
    $PG(1,z)$ random variables using a method outlined in [Polson et al.
    (2013)](https://www.tandfonline.com/doi/abs/10.1080/01621459.2013.829001).
    This method is exact, however for large enough $b$ it often suffices
    to use one of the following approximations.
3.  Saddlepoint Approximation: best for intermediate values of $b$.
    Samples $X$ from by first creating an envelope over the $PG(b,z)$
    density consisting of a mixture of right-truncated inverse-Gaussian
    and left-truncated gamma distributions.
4.  Normal Approximation: best for large values of $b$. Samples $X$ from
    a normal distribution with mean and variance determined by the first
    and second moments of the $PG(b,z)$ distribution.

Some of these methods are faster than others, and the implementation in
each packages differs slightly, which also affects speed. The goal of
`jrpg` is to compile the methods that work best (read: fastest) in each
scenario, per my personal experience. For parallel sampling, the `pgR`
package appears to be the only option. `jrpg` does not use parallel code
by personal preference, as I prefer to reserve this for running multiple
chains when needed. Perhaps it will be added in the future.

# Hybrid sampler for the `jrpg` R package

The hybrid sampler for this package uses the following strategy:

<table style="width:88%;">
<caption>Hybrid Sampler for <code>jrpg</code> R package</caption>
<colgroup>
<col style="width: 29%" />
<col style="width: 29%" />
<col style="width: 29%" />
</colgroup>
<thead>
<tr class="header">
<th>Case</th>
<th><code>jrpg</code> Sampler</th>
<th>Other Implementations</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><span class="math inline">0 &lt; <em>b</em> &lt; 1</span></td>
<td>Truncated Sum-of-Gammas</td>
<td><ul>
<li><p><code>BayesLogit</code>: Truncated Sum-of-Gammas</p></li>
<li><p><code>pg</code>: Truncated Sum-of-Gammas</p></li>
<li><p><code>pgdraw</code>: Devroye Method</p></li>
<li><p><code>pgR</code>: Devroye Method</p></li>
</ul></td>
</tr>
<tr class="even">
<td><span
class="math inline"><em>b</em> <em>i</em><em>n</em>1, 2, …, 13</span></td>
<td>Devroye Method</td>
<td><ul>
<li><p><code>BayesLogit</code>: Truncated Sum-of-Gammas for <span
class="math inline"><em>b</em> ∈ 1, 2</span>; Devroye Method
otherwise</p></li>
<li><p><code>pg</code>: Truncated Sum-of-Gammas for <span
class="math inline"><em>b</em> ∈ 1, 2</span>; Devroye Method
otherwise</p></li>
<li><p><code>pgdraw</code>: Devroye Method</p></li>
<li><p><code>pgR</code>: Devroye Method</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><span class="math inline">1 &lt; <em>b</em> &lt; 13</span>
(non-integer)</td>
<td>Devroye Method + Truncated Sum-of-Gammas</td>
<td><ul>
<li><p><code>BayesLogit</code>: Truncated Sum-of-Gammas</p></li>
<li><p><code>pg</code>: Truncated Sum-of-Gammas</p></li>
<li><p><code>pgdraw</code>: N/A</p></li>
<li><p><code>pgR</code>: Devroye Method</p></li>
</ul></td>
</tr>
<tr class="even">
<td><span class="math inline">13 ≤ <em>b</em> &lt; 170</span></td>
<td>Saddlepoint Approximation</td>
<td><ul>
<li><p><code>BayesLogit</code>: Saddlepoint Approximation</p></li>
<li><p><code>pg</code>: Saddlepoint Approximation</p></li>
<li><p><code>pgdraw</code>: Devroye Method (integers only)</p></li>
<li><p><code>pgR</code>: Devroye Method</p></li>
</ul></td>
</tr>
<tr class="odd">
<td><span class="math inline"><em>b</em> ≥ 170</span></td>
<td>Normal Approximation</td>
<td><ul>
<li><p><code>BayesLogit</code>: Normal Approximation</p></li>
<li><p><code>pg</code>: Normal Approximation</p></li>
<li><p><code>pgdraw</code>: Devroye Method (integers only)</p></li>
<li><p><code>pgR</code>: Normal Approximation</p></li>
</ul></td>
</tr>
</tbody>
</table>

Hybrid Sampler for `jrpg` R package

The approach is identical to that used in the
[PolyaGammaHybridSamplers.jl](https://github.com/wzhorton/PolyaGammaHybridSamplers.jl/tree/main)
package. The function corresponding to the hybrid sampler is this
package is `jrpg()`.

## Examples

For some examples on how to use `jrpg` to facilitate Gibbs sampling in
custom MCMC routines, see the package vignette [Custom MCMC Samplers
with jrpg](inst/doc/mcmc-with-jrpg.html).

## Installation

You can install the development version of `jrpg` like so:

``` r
install.packages('devtools')
devtools::install_github('jacobenglert/jrpg')
```

Note that installation will require compilation of `Rcpp` and
`RcppArmadillo` code. If you have used packages with `Rcpp` dependence,
you will probably not have issues. If not, you will need to install a
working C++ compiler. This will be either
[Xcode](https://developer.apple.com/xcode/) (Mac) or
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) (Windows).

## References

1.  Balamuta, J. J. (2021). Bayesian estimation of restricted latent
    class models: Extending priors, link functions, and structural
    models. University of Illinois Urbana-Champaign. Available at:
    <https://www.ideals.illinois.edu/items/121209>.

2.  Enes Makalic and Daniel Schmidt (2016). High-Dimensional Bayesian
    Regularised Regression with the BayesReg Package. arXiv. Available
    at: <https://arxiv.org/abs/1611.06649>.

3.  Polson, N.G., Scott, J.G. and Windle, J. (2013) ‘Bayesian Inference
    for Logistic Models Using Pólya–Gamma Latent Variables’, *Journal of
    the American Statistical Association*, 108(504), pp. 1339–1349.
    Available at: <https://doi.org/10.1080/01621459.2013.829001>.

4.  Windle, J. (2013) *Forecasting High-Dimensional, Time-Varying
    Variance-Covariance Matrices with High-Frequency Data and Sampling
    Pólya-Gamma Random Variates for Posterior Distributions Derived from
    Logistic Likelihoods*. The University of Texas at Austin.

5.  Windle, J., Polson, N.G. and Scott, J.G. (2014) ‘Sampling
    Polya-Gamma random variates: alternate and approximate techniques’.
    arXiv. Available at: <http://arxiv.org/abs/1405.0506>.
