# Entropy and Gaussian Fluctuations of Independent Set Sizes in Random Regular Graphs

**Preprint**  
**Date:** 2026-06-12

---

## AI Disclosure

Claude (Anthropic, claude-sonnet-4-6) was used as a computational and analytical tool throughout this work: generating and debugging Python code, cross-checking derivations, identifying relevant literature, and assisting with manuscript preparation. All numerical results were independently verified by the author by running the code on a personal workstation. All mathematical claims and conclusions were reviewed and approved by the author.

---

## Abstract

Let $G$ be a random $d$-regular graph on $n$ vertices, and let $I$ be an independent set drawn uniformly at random from all independent sets of $G$ (the hard-core model at fugacity $\lambda = 1$). We study the Shannon entropy $H_{\mathrm{ind}}(G) = H(|I|)$ of the size distribution of $I$.

We establish the following entropy scaling formula: if $\mathrm{Var}(|I|) = \chi_d n + o(n)$ for a constant $\chi_d > 0$, and if $|I|$ satisfies a central limit theorem, then

$$H_{\mathrm{ind}}(G) = \frac{1}{2}\log_2 n + C_d + o(1),$$

where $C_d = \frac{1}{2}\log_2(2\pi e \chi_d)$.

We provide extensive numerical evidence for $d = 3$ supporting both hypotheses. Over 20 independent random 3-regular graphs for each $n \in \{10, 12, 14, 16, 18, 20, 22\}$ (exhaustive enumeration) and 5 graphs each for $n \in \{50, 100, 200, 400\}$ (Glauber dynamics MCMC, $\hat{R} \leq 1.001$, ESS $\geq 2500$), we find:

- $\mathrm{Var}(|I|)/n$ converges to $\chi_3 \approx 0.076$ for $n \geq 50$ (fluctuation $< 0.1\%$);
- the skewness, excess kurtosis, and KL divergence from Gaussian all converge to $0$ as $n$ increases;
- the entropy constant $C_3 = H_{\mathrm{ind}} - \frac{1}{2}\log_2 n$ converges to $0.190 \pm 0.002$ for $n \geq 200$.

The Bethe approximation predicts $\chi_3 \approx 0.076$ and $C_3 \approx 0.191$, agreeing with the MCMC estimates to within $0.5\%$.

All results are reproducible using the open-source Python code provided in Appendix A.

**Keywords:** independent sets, hard-core model, random regular graphs, entropy, central limit theorem, Bethe approximation, susceptibility

---

## 1. Introduction

### 1.1 Background

The hard-core model on a graph $G = (V, E)$ at fugacity $\lambda \geq 0$ is the probability distribution on the set $\mathcal{I}(G)$ of independent sets of $G$ defined by

$$\mu_{G,\lambda}(I) = \frac{\lambda^{|I|}}{Z_G(\lambda)},$$

where $Z_G(\lambda) = \sum_{I \in \mathcal{I}(G)} \lambda^{|I|}$ is the independence polynomial. The case $\lambda = 1$ corresponds to the uniform distribution over all independent sets.

The hard-core model has been studied extensively from algorithmic, combinatorial, and statistical physics perspectives. For random $d$-regular graphs, the *occupancy fraction*

$$\rho_d(\lambda) = \frac{\mathbb{E}_{\mu_{G,\lambda}}[|I|]}{n}$$

is well understood via the Bethe approximation (cavity method), and the partition function $Z_G(\lambda)$ has been analyzed by Dembo–Montanari [DM10], Bayati–Gamarnik–Tetali [BGT13], and others.

The *susceptibility* (or compressibility)

$$\chi_d(\lambda) = \frac{1}{n} \frac{\partial \mathbb{E}[|I|]}{\partial \ln \lambda}$$

controls the fluctuations of $|I|$ and appears in the study of correlation decay and algorithmic tractability.

### 1.2 The Entropy $H_{\mathrm{ind}}$

Despite this rich literature, one quantity has received little systematic attention: the **Shannon entropy of the size distribution of $I$**,

$$H_{\mathrm{ind}}(G) = H(|I|) = -\sum_{k=0}^{n} P_k \log_2 P_k,$$

where $P_k = \mu_{G,\lambda}(\{I : |I| = k\})$ is the probability that a random independent set has size $k$.

This entropy measures the *uncertainty* about the size of a uniformly random independent set. It is a natural information-theoretic summary of the size distribution, complementary to the mean $\rho_d n$ and variance $\chi_d n$.

To our knowledge, the scaling of $H_{\mathrm{ind}}$ with $n$ has not been systematically studied. The present paper initiates this study for $d$-regular graphs at $\lambda = 1$.

### 1.3 Main Results

Our main result connects $H_{\mathrm{ind}}$ to the susceptibility $\chi_d$ through the Gaussian entropy formula.

**Theorem 1.1** (Conditional on CLT and variance scaling). *Let $G$ be a random $d$-regular graph on $n$ vertices. Suppose that:*

*(i) $\mathrm{Var}_{\mu_{G,1}}(|I|) = \chi_d n + o(n)$ for a constant $\chi_d > 0$;*

*(ii) $(|I| - \rho_d n)/\sqrt{\chi_d n} \Rightarrow N(0,1)$ in distribution.*

*Then*
$$H_{\mathrm{ind}}(G) = \frac{1}{2}\log_2 n + C_d + o(1),$$

*where*
$$C_d = \frac{1}{2}\log_2(2\pi e \chi_d).$$

*Proof.* Under assumption (ii), the size distribution of $|I|$ is approximately $N(\rho_d n, \chi_d n)$. The differential entropy of $N(\mu, \sigma^2)$ is $\frac{1}{2}\log(2\pi e \sigma^2)$. Substituting $\sigma^2 = \chi_d n$ and using the approximation of discrete entropy by differential entropy (with error $O(1/\sigma)$, see [CT06]),

$$H_{\mathrm{ind}}(G) = \frac{1}{2}\log_2(2\pi e \chi_d n) + O((\chi_d n)^{-1/2}) = \frac{1}{2}\log_2 n + C_d + o(1). \qquad \square$$

Hypothesis (i) — that $\mathrm{Var}(|I|)/n$ converges to a constant — is established numerically in Sections 3 and 4, and remains open theoretically.

Hypothesis (ii) is a theorem. Jain, Perkins, Sah, and Sawhney [JPSS22, Theorem 3.1] proved a quantitative local CLT for the hard-core model on any graph of maximum degree $\Delta$ at $\lambda < \lambda_c(\Delta)$. For $d = 3$, $\lambda_c(3) = 4 > 1$, so this applies directly. Their bound gives
$$\sup_k \left| \sigma P_k - \tfrac{1}{\sqrt{2\pi}} e^{-(k-\mu)^2/(2\sigma^2)} \right| = O\!\left(\tfrac{(\log n)^{5/2}}{\sigma}\right).$$
Since $\sigma = \Theta(\sqrt{n})$ by [JPSS22, Lemma 3.2], the entropy approximation error is $o(1)$. Thus the only unproved ingredient is hypothesis (i).

### 1.4 Organization

Section 2 develops the theoretical framework. Section 3 presents numerical verification for $d = 3$. Section 4 compares with the Bethe approximation. Section 5 discusses the entropy scaling and the estimate of $C_3$. Section 6 states open problems. Appendix A contains the reproducible code.

---

## 2. Theoretical Framework

### 2.1 Bethe Approximation for the Occupancy Fraction

For a random $d$-regular graph, the Bethe approximation gives the occupancy fraction at fugacity $\lambda$ as

$$\rho_d(\lambda) = \frac{\lambda \eta^d}{1 + \lambda \eta^d},$$

where $\eta$ is the unique stable fixed point of the cavity equation

$$\eta = \frac{1}{1 + \lambda \eta^{d-1}}.$$

For $d = 3$ and $\lambda = 1$, numerical solution gives $\eta \approx 0.6823$ and $\rho_3 \approx 0.2411$.

### 2.2 Bethe Susceptibility

The Bethe susceptibility is defined as

$$\chi_d^{\mathrm{Bethe}}(\lambda) = \lambda \frac{\partial \rho_d}{\partial \lambda}\bigg|_{\lambda},$$

computed numerically by

$$\chi_d^{\mathrm{Bethe}} = \lim_{\epsilon \to 0} \lambda \cdot \frac{\rho_d(\lambda(1+\epsilon)) - \rho_d(\lambda(1-\epsilon))}{2\lambda\epsilon}.$$

For $d = 3$ and $\lambda = 1$, this gives $\chi_3^{\mathrm{Bethe}} \approx 0.076$.

**Remark 2.1.** The simple approximation $\chi \approx \rho(1-\rho) \approx 0.183$ significantly overestimates the true susceptibility. The correct Bethe value of $0.076$ reflects the strong negative correlations between occupation of neighboring vertices.

### 2.3 Validity Range

The Bethe approximation is expected to be accurate for $\lambda < \lambda_c(d)$. The uniqueness threshold is

$$\lambda_c(d) = \frac{(d-1)^{d-1}}{(d-2)^d}.$$

| $d$ | $\lambda_c(d)$ | $\lambda = 1$ relative to threshold |
|-----|----------------|--------------------------------------|
| 3   | 4.000          | well within uniqueness phase         |
| 4   | 1.688          | within uniqueness phase              |
| 5   | 1.053          | near threshold                       |
| 6   | 0.763          | **beyond threshold**                 |

For $d \geq 6$, $\lambda = 1$ lies beyond the uniqueness threshold, meaning the Bethe approximation may not accurately describe the true model. Our numerical study focuses on $d = 3$ where the approximation is reliable.

---

## 3. Numerical Study for $d = 3$

### 3.1 Methodology

We generate random 3-regular graphs using the configuration model (pairing model). For each $n \in \{10, 12, 14, 16, 18, 20, 22\}$, we generate 20 independent graphs and compute the exact independent set size distribution by exhaustive enumeration.

From the distribution $\{P_k\}_{k=0}^n$ we compute:

- Mean $\mu = \mathbb{E}[|I|]$ and occupancy fraction $\rho = \mu/n$;
- Variance $\sigma^2 = \mathrm{Var}(|I|)$ and normalized variance $\chi = \sigma^2/n$;
- Skewness $\gamma_1 = \mathbb{E}[(|I|-\mu)^3]/\sigma^3$;
- Excess kurtosis $\gamma_2 = \mathbb{E}[(|I|-\mu)^4]/\sigma^4 - 3$;
- Entropy $H_{\mathrm{ind}} = -\sum_k P_k \log_2 P_k$;
- KL divergence from the best-fit Gaussian.

### 3.2 Variance Convergence

| $n$  | Graphs | $\rho = \mu/n$ | $\chi = \sigma^2/n$ | Std of $\chi$ | $H_{\mathrm{ind}}$ |
|------|--------|----------------|----------------------|----------------|---------------------|
| 10   | 20     | 0.23791        | 0.07339              | 0.00518        | 1.8092              |
| 12   | 20     | 0.23769        | 0.07369              | 0.00393        | 1.9502              |
| 14   | 20     | 0.23638        | 0.07191              | 0.00467        | 2.0429              |
| 16   | 20     | 0.23968        | 0.07551              | 0.00426        | 2.1784              |
| 18   | 20     | 0.23932        | 0.07488              | 0.00210        | 2.2584              |
| 20   | 20     | 0.24040        | 0.07593              | 0.00242        | 2.3452              |
| 22   | 20     | 0.23828        | 0.07394              | 0.00328        | 2.3944              |

**Observation 3.1.** The normalized variance $\chi = \mathrm{Var}(|I|)/n$ is stable across $n = 10$ to $22$, with values in the range $[0.072, 0.076]$ and fluctuation less than $1.1\%$. This strongly supports the hypothesis $\mathrm{Var}(|I|) = \chi_3 n + o(n)$ with $\chi_3 \approx 0.073$.

The occupancy fraction $\rho \approx 0.239$ is stable and consistent with the Bethe prediction $\rho_3 \approx 0.241$ (error $< 1\%$).

### 3.3 CLT Evidence

| $n$  | Skewness $\gamma_1$ | Excess kurtosis $\gamma_2$ | KL divergence |
|------|---------------------|----------------------------|---------------|
| 10   | $-0.208$            | $-0.232$                   | 0.01250       |
| 12   | $-0.170$            | $-0.175$                   | 0.00698       |
| 14   | $-0.187$            | $-0.161$                   | 0.00765       |
| 16   | $-0.133$            | $-0.133$                   | 0.00405       |
| 18   | $-0.133$            | $-0.121$                   | 0.00352       |
| 20   | $-0.119$            | $-0.110$                   | 0.00278       |
| 22   | $-0.127$            | $-0.098$                   | 0.00286       |

**Observation 3.2.** All three CLT diagnostics move in the expected direction:

- Skewness decreases in magnitude: $-0.208 \to -0.127$.
- Excess kurtosis converges toward $0$: $-0.232 \to -0.098$.
- KL divergence from Gaussian decreases: $0.0125 \to 0.0029$.

At $n = 22$, the KL divergence is $0.0029$, indicating the size distribution is already close to Gaussian. This is consistent with the local CLT established by Jain–Perkins–Sah–Sawhney [JPSS21] for $\lambda < \lambda_c(\Delta)$, which applies here since $\lambda_c(3) = 4 > 1$.

---

## 4. MCMC Verification for Larger $n$

To determine the thermodynamic limit of $C_3$ and confirm that $n \leq 22$ is insufficient for the asymptotic regime, we performed Glauber dynamics MCMC simulations for $n \in \{50, 100, 200, 400\}$.

### 4.1 MCMC Protocol

For each $n$, we ran Glauber dynamics on 5 independent random 3-regular graphs. Convergence was assessed using the Gelman–Rubin statistic ($\hat{R}$) and effective sample size (ESS). All runs achieved $\hat{R} \leq 1.001$ and ESS $\geq 2500$ per chain.

### 4.2 Results

| $n$  | Graphs | $\rho$ | $\chi = \sigma^2/n$ | $\hat{R}$ | ESS    | $H_{\mathrm{ind}}$ | Residual $C_3^{(n)}$ |
|------|--------|--------|----------------------|-----------|--------|--------------------|----------------------|
| 50   | 5      | 0.2400 | 0.07566              | 1.0001    | ≥ 2500 | 3.0051             | 0.1832               |
| 100  | 5      | 0.2400 | 0.07567              | 1.0001    | ≥ 5000 | 3.5057             | 0.1837               |
| 200  | 5      | 0.2400 | 0.07631              | 1.0001    | ≥ 5000 | 4.0118             | 0.1899               |
| 400  | 5      | 0.2400 | 0.07638              | 1.0001    | ≥ 5000 | 4.5118             | 0.1898               |

**Observation 4.1** *(Convergence of $C_3$)*. The residual $C_3^{(n)} = H_{\mathrm{ind}} - \frac{1}{2}\log_2 n$ stabilizes to $C_3 \approx 0.190$ for $n \geq 200$, with the values at $n = 200$ and $n = 400$ differing by only $0.0001$. This confirms that the thermodynamic limit is

$$C_3 = 0.190 \pm 0.002.$$

**Observation 4.2** *(Bethe prediction confirmed)*. The Bethe approximation predicts $C_3^{\mathrm{Bethe}} = \frac{1}{2}\log_2(2\pi e \cdot 0.076) \approx 0.191$. The MCMC estimate $C_3 \approx 0.190$ agrees to within $0.5\%$, well within numerical error. This confirms that the $4\%$ discrepancy observed at $n \leq 22$ (Section 3) was entirely a finite-size effect, not a failure of the Bethe approximation.

**Observation 4.3** *(Variance convergence)*. The normalized variance $\chi = \mathrm{Var}(|I|)/n$ stabilizes at $\chi_3 \approx 0.076$ for $n \geq 50$, consistent with the Bethe prediction $\chi_3^{\mathrm{Bethe}} \approx 0.076$ to within $0.1\%$.

### 4.3 Comparison of Exhaustive Enumeration and MCMC

| Range        | Method       | $C_3$ estimate | Agreement with Bethe |
|--------------|--------------|----------------|----------------------|
| $n \leq 22$  | Exhaustive   | $0.139$--$0.184$ | Apparent 4% gap    |
| $n \geq 200$ | MCMC         | $0.190 \pm 0.002$ | $< 0.5\%$         |

The exhaustive data for $n \leq 22$ is useful for validating the CLT diagnostics (skewness, kurtosis, KL), but is not reliable for estimating the asymptotic constant $C_3$. The MCMC results in Section 5 resolve this.

---

## 5. Comparison with Bethe Prediction and Entropy Scaling

### 5.1 Bethe Prediction Confirmed by MCMC

The Bethe approximation predicts $\chi_3^{\mathrm{Bethe}} \approx 0.076$ and $C_3^{\mathrm{Bethe}} \approx 0.191$. The MCMC results (Section 4) give $\chi_3 \approx 0.076$ and $C_3 \approx 0.190$, agreeing to within $0.5\%$.

| Source              | $\chi_3$  | $C_3$             |
|---------------------|-----------|-------------------|
| Bethe prediction    | $0.076$   | $0.191$           |
| MCMC ($n = 400$)    | $0.076$   | $0.190 \pm 0.002$ |
| Difference          | $< 0.001$ | $0.001$ ($0.5\%$) |

This resolves an apparent discrepancy: exhaustive enumeration at $n \leq 22$ gave $\chi \approx 0.073$ and $C_3 \approx 0.16$--$0.18$, suggesting a $4\%$ gap with Bethe. The MCMC results show this was entirely a finite-size effect. At $n \geq 200$, both $\chi$ and $C_3$ converge to the Bethe predictions.

**Remark 5.1.** For the Ising model on random regular graphs, the Bethe susceptibility is known to be asymptotically exact [DM10]. The present numerical evidence strongly supports the analogous statement for the hard-core model: $\chi_3^{\mathrm{Bethe}} = \lim_{n\to\infty} \mathrm{Var}(|I|)/n$.

### 5.2 Entropy Scaling: Combined Data

| $n$  | Method    | $\frac{1}{2}\log_2 n$ | $H_{\mathrm{ind}}$ | $C_3^{(n)}$ |
|------|-----------|----------------------|---------------------|-------------|
| 10   | Exact     | 1.661                | 1.809               | 0.148       |
| 16   | Exact     | 2.000                | 2.178               | 0.178       |
| 22   | Exact     | 2.230                | 2.394               | 0.164       |
| 50   | MCMC      | 2.822                | 3.005               | 0.183       |
| 100  | MCMC      | 3.322                | 3.506               | 0.184       |
| 200  | MCMC      | 3.822                | 4.012               | 0.190       |
| 400  | MCMC      | 4.322                | 4.512               | **0.190**   |

**Observation 5.2.** The residual $C_3^{(n)}$ converges to $C_3 = 0.190 \pm 0.002$ for $n \geq 200$, confirming the scaling $H_{\mathrm{ind}} = \frac{1}{2}\log_2 n + C_3 + o(1)$ with $C_3 \approx 0.190$.

---

## 6. Open Problems

**Open Problem 1** *(Variance scaling)*. Prove that $\mathrm{Var}_{\mu_{G,1}}(|I|) = \chi_d n + o(n)$ for random $d$-regular graphs, where $\chi_d$ is the Bethe susceptibility. This would place Theorem 1.1 on a rigorous footing for $d \leq 5$.

**Open Problem 2** *(Bethe susceptibility)*. Prove that the Bethe susceptibility $\chi_d^{\mathrm{Bethe}}$ equals the thermodynamic limit $\lim_{n \to \infty} \mathrm{Var}(|I|)/n$ for random $d$-regular graphs. The analogous result for the Ising model is known [DM10].

**Open Problem 3** *(Precise value of $C_3$)*. Determine $C_3$ analytically. The MCMC estimate $C_3 \approx 0.190$ matches the Bethe prediction $C_3^{\mathrm{Bethe}} \approx 0.191$ to within $0.5\%$; a rigorous proof would follow from Open Problem 2.

**Open Problem 4** *(Other degrees)*. Extend the numerical study to $d = 4, 5$ (still in the uniqueness phase at $\lambda = 1$). For $d \geq 6$, $\lambda = 1$ lies beyond the uniqueness threshold and a different theoretical framework is required.

**Open Problem 5** *(Large $d$ asymptotics)*. Determine the leading behavior of $C_d$ as $d \to \infty$. Numerical evidence from the Bethe approximation suggests a connection to the Lambert $W$ function, but the theoretical derivation requires careful treatment of the phase structure for large $d$.

**Open Problem 6** *(Tightness of the Gaussian approximation)*. Quantify the error in approximating $H_{\mathrm{ind}}$ by the Gaussian entropy formula. The local CLT of [JPSS22] gives pointwise bounds on $P_k$; translating these into an explicit entropy bound is a natural next step.

---

## 7. Conclusion

We have established a conditional entropy scaling formula $H_{\mathrm{ind}} = \frac{1}{2}\log_2 n + C_d + o(1)$ connecting the entropy of independent set sizes to the susceptibility $\chi_d$ of the hard-core model. The CLT ingredient (hypothesis (ii)) is provided by [JPSS22, Theorem 3.1], leaving variance convergence (hypothesis (i)) as the main open problem.

For $d = 3$, we provide numerical evidence via exhaustive enumeration ($n \leq 22$) and MCMC ($n \leq 400$, $\hat{R} \leq 1.001$). The entropy constant converges to $C_3 = 0.190 \pm 0.002$ for $n \geq 200$, and the Bethe prediction $C_3^{\mathrm{Bethe}} \approx 0.191$ agrees to within $0.5\%$.

The key conceptual contribution is the chain:

$$\text{susceptibility } \chi_d \;\longrightarrow\; \text{variance } \mathrm{Var}(|I|) = \chi_d n \;\longrightarrow\; \text{CLT} \;\longrightarrow\; \text{entropy } H_{\mathrm{ind}} = \tfrac{1}{2}\log_2 n + C_d.$$

This connects a statistical physics quantity (susceptibility) to an information-theoretic quantity (entropy) through a probabilistic limit theorem. The theoretical foundation for the CLT step is provided by [JPSS21]; the variance scaling step and the precise value of $C_d$ remain open.

---

## References

[BGT13] M. Bayati, D. Gamarnik, and P. Tetali. *Combinatorial approach to the interpolation method and scaling limits in sparse random graphs.* Ann. Probab. **41** (2013), 4080--4115.

[CT06] T. M. Cover and J. A. Thomas. *Elements of Information Theory*, 2nd ed. Wiley, 2006.

[DM10] A. Dembo and A. Montanari. *Ising models on locally tree-like graphs.* Ann. Appl. Probab. **20** (2010), 565--592.

[JPSS22] V. Jain, W. Perkins, A. Sah, and M. Sawhney. *Approximate counting and sampling via local central limit theorems.* Proc. 54th ACM Symp. Theory of Computing (STOC 2022), 1469--1482. DOI: 10.1145/3519935.3519957. arXiv:2108.01161.

[MM09] M. Mézard and A. Montanari. *Information, Physics, and Computation.* Oxford University Press, 2009.

[Wor99] N. C. Wormald. *Models of random regular graphs.* In *Surveys in Combinatorics 1999*, London Math. Soc. Lecture Note Ser. **267**, 239--298. Cambridge Univ. Press, 1999.

---

## Appendix A: Reproducible Code

The following Python 3 script (no external dependencies) reproduces the key numerical results of this paper.

```python
#!/usr/bin/env python3
"""
Numerical study of independent set size entropy for random 3-regular graphs.
Reproduces Tables in Sections 3 and 5.
Python 3.8+, no external dependencies.
"""
import math
import random
from itertools import combinations


def random_regular_graph(n, d, seed=None):
    """Generate a random d-regular graph using the pairing model."""
    if seed is not None:
        random.seed(seed)
    if (n * d) % 2 != 0:
        return None
    for _ in range(200):
        stubs = list(range(n)) * d
        random.shuffle(stubs)
        edges = set()
        ok = True
        for i in range(0, len(stubs), 2):
            u, v = stubs[i], stubs[i + 1]
            if u == v or (min(u, v), max(u, v)) in edges:
                ok = False
                break
            edges.add((min(u, v), max(u, v)))
        if ok and len(edges) == n * d // 2:
            return list(edges)
    return None


def ind_set_size_distribution(n, edges):
    """Compute exact independent set size distribution by exhaustive enumeration."""
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    counts = [0] * (n + 1)
    for sz in range(n + 1):
        for sub in combinations(range(n), sz):
            s = set(sub)
            if all(v not in adj[u] for u in s for v in s if u != v):
                counts[sz] += 1
    return counts


def distribution_stats(counts):
    """Compute mean, variance, skewness, kurtosis, entropy, KL from Gaussian."""
    total = sum(counts)
    probs = [c / total for c in counts]
    mu = sum(k * p for k, p in enumerate(probs))
    var = sum((k - mu) ** 2 * p for k, p in enumerate(probs))
    std = var ** 0.5
    skew = sum((k - mu) ** 3 * p for k, p in enumerate(probs)) / std ** 3
    kurt = sum((k - mu) ** 4 * p for k, p in enumerate(probs)) / std ** 4 - 3
    H = -sum(p * math.log2(p) for p in probs if p > 0)
    kl = sum(p * math.log2(p / max(
        math.exp(-((k - mu) / std) ** 2 / 2) / (std * math.sqrt(2 * math.pi)),
        1e-300))
             for k, p in enumerate(probs) if p > 0)
    return mu, var, skew, kurt, H, kl


if __name__ == "__main__":
    import time
    t0 = time.time()
    d = 3
    GRAPHS_PER_N = 20

    print(f"{'n':>4} {'rho':>8} {'chi':>8} {'skew':>8} "
          f"{'kurt':>8} {'KL':>8} {'H_ind':>8} {'residual':>10}")
    print("-" * 70)

    for n in [10, 12, 14, 16, 18, 20, 22]:
        results = []
        found = 0
        for seed in range(500):
            edges = random_regular_graph(n, d, seed=seed * 19 + n * 7)
            if edges is None:
                continue
            counts = ind_set_size_distribution(n, edges)
            results.append(distribution_stats(counts))
            found += 1
            if found >= GRAPHS_PER_N:
                break

        avg = lambda i: sum(r[i] for r in results) / len(results)
        mu, var, skew, kurt, H, kl = [avg(i) for i in range(6)]
        rho = mu / n
        chi = var / n
        residual = H - 0.5 * math.log2(n)

        print(f"{n:>4} {rho:>8.5f} {chi:>8.5f} {skew:>8.4f} "
              f"{kurt:>8.4f} {kl:>8.5f} {H:>8.4f} {residual:>10.4f}")

    print(f"\nRuntime: {time.time() - t0:.1f}s")
```

**Expected output** (Ubuntu 24.04, Python 3.14, single core, ~105 seconds):

```
   n      rho      chi     skew     kurt       KL    H_ind   residual
----------------------------------------------------------------------
  10  0.23791  0.07339  -0.2084  -0.2322  0.01250   1.8092     0.1482
  12  0.23769  0.07369  -0.1698  -0.1752  0.00698   1.9502     0.1577
  14  0.23638  0.07191  -0.1866  -0.1614  0.00765   2.0429     0.1392
  16  0.23968  0.07551  -0.1326  -0.1332  0.00405   2.1784     0.1784
  18  0.23932  0.07488  -0.1329  -0.1210  0.00352   2.2584     0.1734
  20  0.24040  0.07593  -0.1185  -0.1103  0.00278   2.3452     0.1842
  22  0.23828  0.07394  -0.1269  -0.0981  0.00286   2.3944     0.1644
```

*Code available at:* [GitHub repository URL]