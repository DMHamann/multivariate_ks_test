# Multivariate Kolmogorov–Smirnov Test

Implementation and simulation study for the bivariate extension of the Kolmogorov–Smirnov goodness-of-fit test, based on Justel, Peña & Zamar (1997).

## Background

The classical Kolmogorov–Smirnov test makes use of the ordering of the
real line to define an empirical CDF and compare it against a reference CDF (depending on $H_0$). 
Considering the multivariate case, or even $\mathbb{R}^p$ for $p \geq 2$, has no such ordering.
Thus, a multivariate extension must be constructed differently. In addition, critical values have to be obtained via simulation and cannot be found in a lookup table.

## The statistic

Let the sample $x_1, \dots, x_n$ of i.i.d. $p$-dimensional random variables with distribution $F$ be given. For the test $H_0 : F = F_0$ against $H_1 : F \neq F_0$,

$$D_n = \max_{j=1,2,\dots} d_n^j$$

is the multivariate version of the Kolmogorov-Smirnov test statistic. Here $d_n^j = \sup_{y^j}|G_n(y^j) - y_1^j \cdots y_p^j|$ and

$$y_1^j = F(z_1^j)$$

$$y_i^j = F(z_i^j \mid z_{i-1}^j, \dots, z_1^j)$$

is the sequence of transformations (Rosenblatt transformation) that generates all $(y_1^j, \dots, y_p^j)$. $(z_1^j, \dots, z_p^j)$ for $j = 1, \dots, p!$ is the $j$-th permutation of $(x_1, \dots, x_p)$.

## Further information

The full write-up and description of the statistic and the algorithm can be found [`here`](docs/write_up.md). This is a translation of my [`German version`](docs/ks_test_write_up_ger.pdf).

## Usage

See the notebooks for critical-value simulation and full power-study
examples.

## References
The full citation of the paper my simulation is based on can be found [`here`](docs/references.md).
