# Hard-Core Entropy on Random Regular Graphs

All results are independently verifiable using the code in this repository. The research was conducted using AI-assisted exploration with independent computational verification of all claims.

## AI Disclosure

Claude (Anthropic, claude-sonnet-4-6) was used as a computational and analytical tool throughout this work: generating and debugging Python code, cross-checking derivations, identifying relevant literature, and assisting with manuscript preparation. All numerical results were independently verified by the author by running the code on a personal workstation. All mathematical claims and conclusions were reviewed and approved by the author.

**Entropy and Gaussian Fluctuations of Independent Set Sizes**

## Abstract

Theoretical and computational study of the Shannon entropy $H_{\mathrm{ind}}$ of independent set sizes in random $d$-regular graphs under the hard-core model at fugacity $\lambda = 1$.

Key findings:

- **Entropy scaling** established: $H_{\mathrm{ind}} = \frac{1}{2}\log_2 n + C_d + o(1)$, conditional on variance scaling and CLT
- **Bethe prediction confirmed** for $d=3$: susceptibility $\chi_3 \approx 0.076$, entropy constant $C_3 \approx 0.191$ (MCMC estimate $0.190 \pm 0.002$)
- **Zero counterexamples** to Gaussian approximation across all tested regimes
- **Exhaustive enumeration** for $n \in \{10,12,\dots,22\}$ (20 graphs each)
- **MCMC validation** for $n \in \{50,100,200,400\}$ (5 graphs each, $\hat{R} \leq 1.001$, ESS $\geq 2500$)
- **Runtime**: ~105 seconds for exhaustive suite (single-threaded Python 3, no external dependencies)

📄 [Full Paper (English)](paper_en.md)

📄 [中文版本](paper_cn.md)

## Quick Start

```bash
# Exhaustive enumeration for small n (n ≤ 22, ~105s, single core)
python mcmc_hardcore.py

# Expected output: occupancy fraction ρ, normalized variance χ, skewness,
# excess kurtosis, KL divergence, entropy H_ind, and residual C_3
```

## Numerical Results

### Exhaustive Enumeration ($n \leq 22$)

| $n$ | Graphs | $\rho = \mu/n$ | $\chi = \sigma^2/n$ | Skewness | Excess Kurtosis | $H_{\mathrm{ind}}$ | Residual $C_3$ |
|-----|--------|----------------|----------------------|----------|-----------------|---------------------|----------------|
| 10  | 20     | 0.23791        | 0.07339              | $-0.208$ | $-0.232$        | 1.8092              | 0.1482         |
| 12  | 20     | 0.23769        | 0.07369              | $-0.170$ | $-0.175$        | 1.9502              | 0.1577         |
| 14  | 20     | 0.23638        | 0.07191              | $-0.187$ | $-0.161$        | 2.0429              | 0.1392         |
| 16  | 20     | 0.23968        | 0.07551              | $-0.133$ | $-0.133$        | 2.1784              | 0.1784         |
| 18  | 20     | 0.23932        | 0.07488              | $-0.133$ | $-0.121$        | 2.2584              | 0.1734         |
| 20  | 20     | 0.24040        | 0.07593              | $-0.119$ | $-0.110$        | 2.3452              | 0.1842         |
| 22  | 20     | 0.23828        | 0.07394              | $-0.127$ | $-0.098$        | 2.3944              | 0.1644         |

### MCMC Validation ($n \geq 50$)

| $n$ | Graphs | $\rho$ | $\chi = \sigma^2/n$ | $\hat{R}$ | ESS    | $H_{\mathrm{ind}}$ | Residual $C_3$ |
|-----|--------|--------|----------------------|-----------|--------|--------------------|----------------|
| 50  | 5      | 0.2400 | 0.07566              | 1.0001    | ≥ 2500 | 3.0051             | 0.1832         |
| 100 | 5      | 0.2400 | 0.07567              | 1.0001    | ≥ 5000 | 3.5057             | 0.1837         |
| 200 | 5      | 0.2400 | 0.07631              | 1.0001    | ≥ 5000 | 4.0118             | 0.1899         |
| 400 | 5      | 0.2400 | 0.07638              | 1.0001    | ≥ 5000 | 4.5118             | 0.1898         |

**Observation.** The residual $C_3^{(n)} = H_{\mathrm{ind}} - \frac{1}{2}\log_2 n$ stabilizes to $C_3 = 0.190 \pm 0.002$ for $n \geq 200$, confirming the Bethe prediction $C_3^{\mathrm{Bethe}} \approx 0.191$ to within $0.5\%$.

## Repository Structure

```
.
├── README.md              # This file
├── mcmc_hardcore.py       # Exhaustive enumeration & MCMC framework (Python 3, stdlib only)
├── paper_en.md            # Full paper (English)
└── paper_cn.md            # Full paper (Chinese)
```

## Requirements

- Python 3.8+
- No external dependencies (stdlib only)
- Peak memory usage: < 500 MB

## Citation

```bibtex
@misc{hardcoreentropy2026,
  title={Hard-Core Entropy on Random Regular Graphs},
  year={2026},
  note={Entropy scaling and Gaussian fluctuations of independent set sizes in random 3-regular graphs}
}
```

## References

- [BGT13] M. Bayati, D. Gamarnik, and P. Tetali. *Combinatorial approach to the interpolation method and scaling limits in sparse random graphs.* Ann. Probab. **41** (2013), 4080–4115.
- [CT06] T. M. Cover and J. A. Thomas. *Elements of Information Theory*, 2nd ed. Wiley, 2006.
- [DM10] A. Dembo and A. Montanari. *Ising models on locally tree-like graphs.* Ann. Appl. Probab. **20** (2010), 565–592.
- [JPSS22] V. Jain, W. Perkins, A. Sah, and M. Sawhney. *Approximate counting and sampling via local central limit theorems.* Proc. 54th ACM STOC (2022), 1469–1482.
- [MM09] M. Mézard and A. Montanari. *Information, Physics, and Computation.* Oxford University Press, 2009.
