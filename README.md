# Multivariate Kolmogorov–Smirnov Test


Implementation and simulation study of the multivariate Kolmogorov–Smirnov goodness-of-fit test proposed by Justel, Peña, and Zamar (1997).

This repository contains a Python implementation of the test statistic together with Monte Carlo procedures for estimating critical values and evaluating the statistical power under alternative distributions.

The project was developed as part of a Master's thesis in mathematical statistics, focusing on goodness-of-fit testing, empirical processes, and simulation-based inference.

---

# Repository Structure
.
├── docs/                 # Handout in German and Reference to the paper
├── figures/              # Result tables
├── notebooks/            # Jupyter notebooks for experiments
├── src/multivariate_ks_test                
│   └── algorithm.py/     # Main Algorithm
├── requirements.txt      # Python dependencies
├── pyproject.toml        # Project configuration and package metadata
├── README.md             # Project documentation
└── .gitignore            # Files excluded from version control


# Installation

Clone the repository:

git clone <repository-url>
cd <repository-name>

Create and activate a virtual environment:

python -m venv .venv
source .venv/bin/activate      # Linux/macOS
.venv\Scripts\activate         # Windows


Install the required dependencies:

pip install -r requirements.txt

For development, install the package in editable mode:

pip install -e .


# Reproducibility

To reproduce the experiments:

Install the dependencies.
Run the notebooks with according specification of n and critical values.
