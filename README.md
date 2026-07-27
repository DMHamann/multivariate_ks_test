# Multivariate Kolmogorov–Smirnov Test

Implementation, exact-computation algorithm, and simulation study for the
bivariate extension of the Kolmogorov–Smirnov goodness-of-fit test, based on
the finite-representation result of Justel, Peña & Zamar (1997), plus my own
Monte Carlo validation of critical values and test power.

## Why this is a nontrivial problem

The classical Kolmogorov–Smirnov test relies on the natural ordering of the
real line to define an empirical vs. reference CDF comparison. That ordering
has no canonical analogue in $\mathbb{R}^p$ for $p \geq 2$, so a faithful
multivariate extension needs a different construction entirely and,
because no closed-form null distribution exists in the multivariate case,
critical values must be obtained via simulation rather than a lookup table.

## The statistic

For an i.i.d. sample $x_1, \dots, x_n$ from a $p$-dimensional distribution
$F$, testing $H_0: F = F_0$ against $H_1: F \neq F_0$, the multivariate KS
statistic is

$$D_n = \max_{j=1,2,\dots} \; \sup_{y^j} \left| G_n(y^j) - y_1^j \cdots y_p^j \right|,$$

maximized over all $p!$ coordinate permutations combined with the Rosenblatt
transformation. For the bivariate case, this reduces to a maximum over a
finite, sample-determined set of candidate points, making exact computation
tractable — see [`docs/paper_summary.md`](docs/paper_summary.md) for the
full derivation, the five-term computational procedure, and my Monte Carlo
validation against the original article's simulation results.

## What's in this repository

```
├── docs/
│   ├── paper_summary.md        # full write-up: theory, algorithm, simulation results
│   ├── references.md           # citations
│   └── simulation_handout_de.pdf  # original handout (German)
├── notebooks/
│   ├── 01_ks_test_simulation.ipynb    # statistic + power/Type-I-error study
│   └── 02_MC_quantiles_simulation.ipynb  # Monte Carlo critical-value estimation
├── src/multivariate_ks_test/
│   ├── __init__.py
│   └── algorithm.py            # ks_2d_statistic implementation
├── requirements.txt
└── pyproject.toml
```

## Results at a glance

- **Critical values (Type I error):** Monte Carlo estimates of the 95th
  percentile of $D_n$ agree with Justel et al.'s reported values to within
  $\sim 10^{-3}$ across $n \in \{15, 25, 50, 100\}$.
- **Power:** power increases with both sample size and mixture weight
  $\varepsilon$, as expected, and increments closely track the article's;
  absolute power values run consistently lower than the article's by up to
  ~0.06 — an unresolved, systematic discrepancy discussed openly in
  [`docs/paper_summary.md`](docs/paper_summary.md#discussion).

## Usage

```python
from multivariate_ks_test.algorithm import ks_2d_statistic, G_uniform

Dn = ks_2d_statistic(X, Y, G_uniform)
```

See the notebooks for critical-value simulation and full power-study
examples.

## References

Full citation and related work in [`docs/references.md`](docs/references.md).
