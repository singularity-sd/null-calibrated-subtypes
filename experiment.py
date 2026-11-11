"""Are cognitive subtypes real? Null-calibrated consensus clustering for test batteries.

From the Singularity v3 forecast T3v3-2026-0038 in the view v3-classic-U: consensus clustering ->
Experimental and Cognitive Psychology. Subtyping studies cluster test scores and report profiles. We ask how often
common procedures report subtypes in data that have none (one continuous ability factor, skew, ceiling effects), and
whether consensus clustering calibrated against a correlation-preserving null (the M3C idea from genomics) avoids it.

Fixed before the single run: n = 300, 8 tests, factor loading 0.6, scenarios below, 40 datasets per scenario,
k = 2..5, 25 resamples at 80%, 10 null references (subtypes declared only if the observed PAC is below all of them), seed.
Run: python experiment.py -> results.json, fig_rates.png
"""

import json
from pathlib import Path

import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(20261105)
N, P, LOAD, R = 300, 8, 0.6, 40
KS = [2, 3, 4, 5]
B, NREF, ALPHA = 25, 10, 0.05


def one_factor(n, skew=False, ceiling=False):
    f = RNG.standard_normal(n)
    X = LOAD * f[:, None] + np.sqrt(1 - LOAD ** 2) * RNG.standard_normal((n, P))
    if skew:
        X[:, :4] = np.exp(0.8 * X[:, :4])
    if ceiling:
        X[:, 4:] = np.minimum(X[:, 4:], 0.5)
    return X


def mixture(n, k, d):
    X = one_factor(n)
    g = RNG.integers(0, k, n)
    shift = np.zeros((k, P))
    for c in range(1, k):
        shift[c, (c - 1) * 3 % P:(c - 1) * 3 % P + 3] = d     # each subtype stronger on a different block of tests
    return X + shift[g]


def standardise(X):
    return (X - X.mean(0)) / X.std(0)


def silhouette_rule(X):
    best = max(KS, key=lambda k: silhouette_score(X, KMeans(k, n_init=3, random_state=0).fit_predict(X)))
    s = silhouette_score(X, KMeans(best, n_init=3, random_state=0).fit_predict(X))
    return s > 0.25, best          # Kaufman-Rousseeuw: above 0.25, "weak structure" at least


def lpa_bic(X, cov):
    bics = [GaussianMixture(k, covariance_type=cov, n_init=2, random_state=0).fit(X).bic(X) for k in [1] + KS]
    k = ([1] + KS)[int(np.argmin(bics))]
    return k > 1, k


def pac(X, k):
    n = len(X)
    M = np.zeros((n, n)); I = np.zeros((n, n))
    for _ in range(B):
        idx = RNG.choice(n, int(0.8 * n), replace=False)
        lab = KMeans(k, n_init=1, random_state=int(RNG.integers(1e9))).fit_predict(X[idx])
        same = lab[:, None] == lab[None, :]
        M[np.ix_(idx, idx)] += same
        I[np.ix_(idx, idx)] += 1
    C = M[I > 0] / I[I > 0]
    return float(np.mean((C > 0.1) & (C < 0.9)))


def consensus_test(X):
    real = {k: pac(X, k) for k in KS}
    cov = np.cov(X, rowvar=False)
    null = {k: [] for k in KS}
    for _ in range(NREF):
        Z = standardise(RNG.multivariate_normal(np.zeros(P), cov, len(X)))
        for k in KS:
            null[k].append(pac(Z, k))
    best = min(KS, key=lambda k: real[k] - np.mean(null[k]))
    p = (1 + sum(z <= real[best] for z in null[best])) / (NREF + 1)
    return p < ALPHA + 1e-9 if NREF >= 19 else bool(real[best] < min(null[best])), best


SCEN = {"null: one factor": lambda: one_factor(N), "null: skewed tests": lambda: one_factor(N, skew=True),
        "null: ceiling effect": lambda: one_factor(N, ceiling=True),
        "2 subtypes, d = 1.0": lambda: mixture(N, 2, 1.0), "2 subtypes, d = 1.5": lambda: mixture(N, 2, 1.5),
        "2 subtypes, d = 2.0": lambda: mixture(N, 2, 2.0), "3 subtypes, d = 1.5": lambda: mixture(N, 3, 1.5),
        "3 subtypes, d = 2.0": lambda: mixture(N, 3, 2.0)}
METHODS = {"k-means + silhouette > 0.25": silhouette_rule, "latent profiles (diagonal), BIC": lambda X: lpa_bic(X, "diag"),
           "Gaussian mixture (full), BIC": lambda X: lpa_bic(X, "full"), "consensus + correlated null": consensus_test}


def main():
    res = {"protocol": {"N": N, "P": P, "LOAD": LOAD, "R": R, "KS": KS, "B": B, "NREF": NREF}, "rates": {}}
    for sname, gen in SCEN.items():
        hits = {m: 0 for m in METHODS}
        for _ in range(R):
            X = standardise(gen())
            for m, f in METHODS.items():
                hits[m] += bool(f(X)[0])
        res["rates"][sname] = {m: h / R for m, h in hits.items()}
        print(sname, res["rates"][sname], flush=True)
        (HERE / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    figure(res)


def figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sc = list(res["rates"])
    fig, ax = plt.subplots(figsize=(11, 4.6))
    w = 0.2
    for j, m in enumerate(METHODS):
        ax.bar(np.arange(len(sc)) + (j - 1.5) * w, [res["rates"][s][m] for s in sc], w, label=m)
    ax.axvspan(-0.5, 2.5, color="0.92", zorder=0)
    ax.text(1, 1.02, "no subtypes: any bar is a false discovery", ha="center", fontsize=9)
    ax.text(5, 1.02, "real subtypes: higher is better", ha="center", fontsize=9)
    ax.set_xticks(range(len(sc)))
    ax.set_xticklabels([s.replace(": ", ":\n").replace(", ", ",\n") for s in sc], fontsize=8)
    ax.set_ylabel("share of datasets where subtypes are reported")
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4)
    fig.tight_layout()
    fig.savefig(HERE / "fig_rates.png", dpi=160)


if __name__ == "__main__":
    main()
