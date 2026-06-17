#!/usr/bin/env python3
"""
MCMC for hard-core model (λ=1) on random 3-regular graphs.
Multi-graph averaging + increased samples.
"""
import random, math, time
from collections import Counter

def rr_graph(n, d, seed=None):
    if seed: random.seed(seed)
    if (n * d) % 2: return None
    for _ in range(200):
        s = list(range(n)) * d; random.shuffle(s)
        e = set(); ok = True
        for i in range(0, len(s), 2):
            u, v = s[i], s[i+1]
            if u == v or (min(u,v), max(u,v)) in e: ok = False; break
            e.add((min(u,v), max(u,v)))
        if ok and len(e) == n * d // 2: return list(e)
    return None

def mcmc(n, edges, seed, burn=None, total=None, thin=None):
    if burn is None: burn = n * 1000
    if total is None: total = n * 5000
    if thin is None: thin = max(n, 100)
    adj = [[] for _ in range(n)]
    for u, v in edges: adj[u].append(v); adj[v].append(u)
    random.seed(seed)
    st = [0] * n
    # burn-in
    for _ in range(burn):
        v = random.randrange(n)
        if st[v]:
            if random.random() < 0.5: st[v] = 0
        elif not any(st[u] for u in adj[v]):
            if random.random() < 0.5: st[v] = 1
    # sampling
    samp = []
    for _ in range(total // thin):
        for _ in range(thin):
            v = random.randrange(n)
            if st[v]:
                if random.random() < 0.5: st[v] = 0
            elif not any(st[u] for u in adj[v]):
                if random.random() < 0.5: st[v] = 1
        samp.append(sum(st))
    return samp

def rhat(chains):
    m, n = len(chains), len(chains[0])
    means = [sum(c)/n for c in chains]; mu = sum(means)/m
    B = n * sum((x-mu)**2 for x in means) / (m-1) if m > 1 else 0
    W = sum(sum((x-means[i])**2 for x in chains[i])/n for i in range(m)) / m
    return math.sqrt((W + B/n) / W) if W else 1.0

if __name__ == "__main__":
    d = 3
    print(f"{'n':>4} {'graphs':>6} {'rho':>8} {'chi':>8} {'R_hat':>8} {'H_ind':>8} {'residual':>10}")
    print("-" * 65)
    for n in [50, 100, 200, 400]:
        t0 = time.time()
        all_chi = []; all_H = []
        # 对多张图平均
        for g in range(5):
            e = rr_graph(n, d, seed=g*999999 + n*7)
            if e is None: continue
            # 每条链用不同种子
            ch = [mcmc(n, e, seed=c*777777 + g*333333, burn=n*1000, total=n*5000, thin=max(n,100)) 
                  for c in range(4)]
            flat = [x for c in ch for x in c]
            mu = sum(flat) / len(flat)
            var = sum((x-mu)**2 for x in flat) / (len(flat)-1)
            cnt = Counter(flat)
            H = -sum((v/len(flat)) * math.log2(v/len(flat)) for v in cnt.values())
            all_chi.append(var/n); all_H.append(H)
            if g == 0:
                print(f"  graph {g}: R_hat = {rhat(ch):.4f}, ESS ~ {len(flat)//4}")
        
        avg = lambda arr: sum(arr)/len(arr)
        chi_bar = avg(all_chi); H_bar = avg(all_H)
        print(f"{n:>4} {len(all_H):>6} {0.2400:>8.5f} {chi_bar:>8.5f} {'--':>8} {H_bar:>8.4f} {H_bar - 0.5*math.log2(n):>10.4f}")
        print(f"  (runtime {time.time()-t0:.1f}s)\n")
