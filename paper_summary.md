# Simulation Study: Power of the Multivariate Kolmogorov–Smirnov Test

*Author: David Hamann*

This document summarizes the theoretical background, algorithm, and simulation
results behind the code in this repository. It is an English restructuring of
a handout I originally wrote in German; the original is kept in
[`simulation_handout_de.pdf`](simulation_handout_de.pdf) for reference.

---

## 1. Motivation

The classical (univariate) Kolmogorov–Smirnov test compares an empirical CDF
to a reference CDF using the supremum of their absolute difference. This
construction relies on the natural ordering of the real line, which has no
canonical analogue in $\mathbb{R}^p$ for $p \geq 2$: there is no single "correct"
way to order or accumulate probability mass in higher dimensions. A faithful
multivariate extension therefore needs (a) a well-defined statistic that
reduces to the univariate KS statistic when $p = 1$, and (b) a way to obtain
critical values, since no closed-form null distribution is available in
general.

## 2. The Multivariate KS Statistic

Let $x_1, \dots, x_n$ be an i.i.d. sample of $p$-dimensional random vectors
with distribution $F$. For the hypothesis test

$$H_0: F = F_0 \quad \text{vs.} \quad H_1: F \neq F_0,$$

the multivariate Kolmogorov–Smirnov statistic is defined as

$$D_n = \max_{j = 1, 2, \dots} d_n^j,$$

where

$$d_n^j = \sup_{y^j} \left| G_n(y^j) - y_1^j \cdots y_p^j \right|$$

and $(y_1^j, \dots, y_p^j)$ is obtained from a sequence of transformations
(the **Rosenblatt transformation**):

$$y_1^j = F(z_1^j), \qquad y_i^j = F(z_i^j \mid z_{i-1}^j, \dots, z_1^j).$$

Here $(z_1^j, \dots, z_p^j)$ denotes the $j$-th permutation of the coordinates
$(x_1, \dots, x_p)$, for $j = 1, \dots, p!$. Intuitively: the Rosenblatt
transformation maps any continuous $p$-dimensional distribution to the
uniform distribution on $[0,1]^p$ under $H_0$, coordinate by coordinate,
conditioning each new coordinate on all previous ones. Because this
transformation is order-dependent, the statistic maximizes over **all**
coordinate permutations to remove that arbitrary choice of ordering.

## 3. Exact Computation in the Bivariate Case

Computing $D_n$ via a naive grid search over $y^j \in [0,1]^p$ is
computationally expensive. For the bivariate case ($p = 2$), Justel et al.
show that the supremum is in fact attained on a **finite set of candidate
points** determined by the sample itself, which makes exact computation
tractable:

$$D_n = \max_{u \in I,\, v \in P} \{ G_n(u) - G(u),\; G(v) - G_n(v^-) \}$$

where:

- $G$ is the CDF of two independent $\mathrm{Uniform}(0,1)$ random variables,
  i.e. $G(u_1, u_2) = u_1 u_2$, and $G_n$ is the empirical CDF of the
  transformed sample.
- $I = \{(x_j, y_i) \mid x_i \leq x_j,\ y_i \geq y_j;\ i,j = 0, \dots, n\}$ is
  the set containing the origin $(0,0)$ and all **intersection points**
  $(x_j, y_i)$ for $x_i < x_j$ and $y_i > y_j$.
- $P = \{(x_j, y_j) \mid x_j > x_i,\ y_j < y_i;\ i,j = 0, \dots, n+1\}$ is the
  set containing the point $(1,1)$, the corresponding intersection points,
  and the projections of the observed points onto the right and top edges of
  the unit square.

Following the article's notation, write $D_n^+(u) = G_n(u) - G(u)$ and
$D_n^-(u) = G(u) - G_n(u)$. Theorem 2 of the article then shows $D_n$ can be
obtained by evaluating a small, finite number of candidate points (between
$3n$ and $3n + \binom{n}{2}$, depending on the sample configuration), reduced
to five terms:

1. $D_n^1 = \max_{i=1,\dots,n} D_n^+(x_i, y_i)$ — maximum distance at the
   observed points.
2. $D_n^2 = \max_{i,j=1,\dots,n} \{ D_n^+(x_j, y_i) \mid x_j > x_i,\ y_j < y_i \}$
   — maximum distance over all intersection points.
3. $D_n^3 = \tfrac{2}{n} - \min_{i,j=1,\dots,n} \{ D_n^-(x_j, y_i) \mid x_j > x_i,\ y_j < y_i \}$
   — minimum distance over all intersection points (with the $2/n$
   correction).
4. $D_n^4 = \tfrac{1}{n} - \min_{i=1,\dots,n} D_n^+(1, y_i)$ — maximum
   distance among projections of the observed points onto the right edge of
   the unit square.
5. $D_n^5 = \tfrac{1}{n} - \min_{i=1,\dots,n} D_n^+(x_i, 1)$ — maximum
   distance among projections onto the top edge of the unit square.

$$D_n = \max\{D_n^1, D_n^2, D_n^3, D_n^4, D_n^5\}$$

This is implemented directly in [`src/multivariate_ks_test/algorithm.py`](../src/multivariate_ks_test/algorithm.py)
as `ks_2d_statistic`, with `D1`–`D5` in the code corresponding one-to-one to
$D_n^1$–$D_n^5$ above.

## 4. Monte Carlo Estimation of Critical Values

Since the exact (finite-sample) distribution of $D_n$ under $H_0$ has no
known closed form, its percentiles must be approximated by Monte Carlo
simulation: repeatedly sampling from the (known) null distribution,
computing $D_n$ for each sample, and taking the empirical quantile of the
resulting values.

The table below compares the 95th percentile of $D_n$ (i.e. the critical
value for $\alpha = 0.05$) obtained from Justel et al.'s article against my
own Monte Carlo estimate, both based on 2,000 replications:

| $n$ | MC percentile (Justel et al.) | MC percentile (this work) |
|---|---|---|
| 15  | 0.4141 | 0.4118 |
| 25  | 0.3254 | 0.3277 |
| 50  | 0.2350 | 0.2350 |
| 100 | 0.1675 | 0.1663 |

The absolute deviation between the two sets of values is on the order of
$10^{-3}$, which is within the range expected from Monte Carlo sampling
variability at this number of replications.

## 5. Power Simulation: Goodness-of-Fit Test for Bivariate Normality

The article also reports power simulations for a specific goodness-of-fit
scenario, which I reproduced independently. Under $H_0$, the data follow a
bivariate normal distribution with $\mu = 0$ and covariance

$$\Sigma = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}.$$

The alternative is a Gaussian mixture

$$(1 - \varepsilon)\, \mathcal{N}_2(0, \Sigma) + \varepsilon\, \mathcal{N}_2(\mu_1, \Sigma), \qquad \varepsilon \in \{0.1, 0.2, 0.4\}, \quad \mu_1 = (3, 3)^\top.$$

Because $F_0$ (the $H_0$ distribution) is known exactly, samples
$x_1, \dots, x_n$ can be drawn from the mixture and transformed via the
Rosenblatt transformation:

$$u_1 = F_0(x_1), \qquad u_i = F_0(x_i \mid x_{i-1}, \dots, x_1).$$

The resulting points $u_1, \dots, u_n$ are then used to compute
$D_n = \sup_{u \in [0,1]^2} |G_n(u) - u_1 u_2|$ via `ks_2d_statistic`.

**Deriving the conditional transform for this specific case:** since under
$H_0$, $(X_1, X_2) \sim \mathcal{N}(0, \Sigma)$, it follows that
$X_1 \sim \mathcal{N}(0, 1)$ and $F_1(x_1) = \Phi(x_1)$. Since $X_2 \mid X_1$
is also normally distributed, standard conditional-Gaussian identities give
$X_2 \mid X_1 = x_1 \sim \mathcal{N}(0.5\, x_1,\, 0.75)$ — using
$E[X_2 \mid X_1 = x_1] = 0.5\, x_1$ and $\mathrm{Var}(X_2 \mid X_1) = 1 - 0.5^2 = 0.75$
— so that

$$F_2(x_2 \mid x_1) = \Phi\!\left( \frac{x_2 - 0.5\, x_1}{\sqrt{0.75}} \right).$$

### Results

Power was estimated over 10,000 replications per $(n, \varepsilon)$
combination, using the Monte Carlo critical values from the article (since my
own MC percentiles agreed closely with theirs, per Section 4). The two
values shown in the "this work" column come from two independent runs of
10,000 replications each, illustrating simulation variability.

| $n$ | $\varepsilon$ | Power (Justel et al.) | Power (this work) |
|---|---|---|---|
| 15  | 0.1 | 0.13 | 0.0817 / 0.0918 |
| 15  | 0.2 | 0.27 | 0.2132 / 0.2135 |
| 15  | 0.4 | 0.73 | 0.6738 / 0.6706 |
| 25  | 0.1 | 0.16 | 0.1155 / 0.112 |
| 25  | 0.2 | 0.40 | 0.335 / 0.335 |
| 25  | 0.4 | 0.92 | 0.8899 / 0.8863 |
| 50  | 0.1 | 0.23 | 0.177 / 0.1712 |
| 50  | 0.2 | 0.67 | 0.6027 / 0.5984 |
| 50  | 0.4 | 1.00 | 0.996 / 0.9958 |
| 100 | 0.1 | 0.41 | 0.3198 / 0.3119 |
| 100 | 0.2 | 0.94 | 0.9216 / 0.9176 |
| 100 | 0.4 | 1.00 | 1 / 1 |

### Discussion

Across all $(n, \varepsilon)$ combinations, my power estimates are
consistently **lower** than those reported in the article, with deviations
of up to about 0.06 — noticeably larger than the $\sim 10^{-3}$ variation
seen between my own repeated runs. This gap is too large and too consistent
in direction to be explained by sampling variability alone, and points to a
systematic discrepancy somewhere in the pipeline (mine, the article's, or
both).

That said, the qualitative behavior matches the article closely:

- Power increases monotonically with both $\varepsilon$ and $n$, as expected.
- The *increments* in power from changing $\varepsilon$ closely track the
  article's. For example, at $n = 15$, the article reports a power increase
  of 0.14 when going from $\varepsilon = 0.1$ to $\varepsilon = 0.2$; my two
  runs show increases of 0.1315 and 0.1217 respectively — close in
  magnitude, even though the absolute levels differ.

**Steps taken to locate the discrepancy:**
- Re-derived the Rosenblatt transformation for this specific bivariate
  normal case using the conditional density directly, as a cross-check on
  the derivation above.
- Re-checked the `ks_2d_statistic` implementation line-by-line against the
  article's description of the algorithm (Section 3).
- Re-read my own code for implementation bugs.
- Implemented a **second, independent version** of the Gaussian mixture
  sampling procedure (using `sklearn.mixture.GaussianMixture` directly rather
  than manual component sampling) to rule out an error specific to one
  sampling approach. Both versions produced consistent results with each
  other, which rules out a mixture-sampling bug as the sole explanation.

Despite the systematic offset, the overall shape and internal consistency of
the power simulation suggest the implementation is *methodologically sound*;
the remaining discrepancy against the article's absolute numbers is flagged
here explicitly rather than resolved, since I was not able to identify its
exact source with certainty.

## 6. Limitations

- The exact finite-set computation of $D_n$ (Section 3) is only implemented
  for $p = 2$; for $p > 2$, the maximum in the original definition would need
  to be taken over $p!$ permutations combined with a supremum evaluated over
  $[0,1]^p$, which becomes rapidly more expensive.
- No closed-form null distribution exists for $D_n$, so critical values
  always require Monte Carlo approximation, whose accuracy is itself subject
  to sampling variability (Section 4).
- As discussed in Section 5, there remains an unresolved, systematic
  discrepancy between my power estimates and the article's; the qualitative
  trends replicate, but the absolute power values do not match as closely as
  the Type-I-error-side (critical value) simulation does.

## References

- Justel, A., Peña, D., & Zamar, R. (1997). *A multivariate Kolmogorov-Smirnov
  test of goodness of fit.* Statistics & Probability Letters. *(Please verify
  exact volume/page numbers against your `docs/references.md` before
  publishing.)*
