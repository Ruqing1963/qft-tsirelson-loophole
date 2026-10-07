"""
Effective Hilbert-space dimension of local measurements in the vacuum of a free lattice field.

Part A  Harmonic chain (lattice free scalar field) H = 1/2 sum [pi_i^2 + m^2 phi_i^2 + (phi_{i+1} - phi_i)^2] on a ring of
        N sites, in its ground state.  For a block of l sites the reduced state is Gaussian, a product of thermal modes
        with symplectic eigenvalues nu_k (eigenvalues of sqrt(X_A P_A)).  Its spectrum is the set of products
        prod_k (1 - q_k) q_k^{n_k} with q_k = (nu_k - 1/2)/(nu_k + 1/2).
Part B  Smooth max-entropy H_0^eps: d(eps) is the smallest number of eigenvalues of rho_A whose sum is at least 1 - eps.
        By the gentle-measurement lemma, every correlation between the block and its complement is reproduced to
        within 2 sqrt(eps) by local measurements compressed to d(eps) dimensions.  We compute d(eps) exactly (best-first
        enumeration of occupation vectors) and compare log2 d(eps) with the entanglement entropy and with the naive
        count of the block's degrees of freedom.

Usage:  python qft_effective_dimension.py [--show]
Figures go to ./figures (PNG, PDF), data to ./data (CSV); override with QE_FIG_DIR, QE_DATA_DIR.
"""
import os
import sys
import heapq
import numpy as np
import matplotlib.pyplot as plt

FIG_DIR = os.environ.get("QE_FIG_DIR", "figures")
DATA_DIR = os.environ.get("QE_DATA_DIR", "data")


def savefig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150)


def savedata(name, columns, cols, meta=()):
    os.makedirs(DATA_DIR, exist_ok=True)
    arr = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    with open(os.path.join(DATA_DIR, f"{name}.csv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in meta:
            fh.write(f"# {line}\n")
        fh.write(",".join(columns) + "\n")
        np.savetxt(fh, arr, delimiter=",", fmt="%.10g")


# ================================================================ Part A: Gaussian block states
def vacuum_correlators(N, m):
    k = 2 * np.pi * np.arange(N) / N
    w = np.sqrt(m ** 2 + 4 * np.sin(k / 2) ** 2)
    r = np.arange(N)
    X = np.array([np.mean(np.cos(k * d) / (2 * w)) for d in r])     # <phi_0 phi_d>
    P = np.array([np.mean(np.cos(k * d) * w / 2) for d in r])       # <pi_0 pi_d>
    return X, P


def block_modes(X, P, l):
    idx = np.abs(np.subtract.outer(np.arange(l), np.arange(l)))
    XA, PA = X[idx], P[idx]
    ev = np.linalg.eigvals(XA @ PA).real
    nu = np.sqrt(np.clip(ev, 0.25, None))
    q = (nu - 0.5) / (nu + 0.5)
    return np.sort(q)[::-1]


def entropy_bits(q):
    q = q[q > 1e-300]
    s = -np.log(1 - q) - q / (1 - q) * np.log(q)
    return float(np.sum(s) / np.log(2))


def smooth_dimension(q, eps, cap=2_000_000):
    """Smallest number d of eigenvalues of prod_k thermal(q_k) with total weight >= 1 - eps (best-first search)."""
    base = float(np.sum(np.log1p(-q)))                       # log of the largest eigenvalue (all modes)
    q = q[q > 1e-3 * eps / max(len(q), 1)]                   # excitations of the dropped modes carry < 1e-3 eps in total
    nlq = -np.log(q)
    # each multiset of excitations is generated once: a child increments a mode index >= the last one incremented
    heap = [(-base, 0)]
    tot, d = 0.0, 0
    while heap and tot < 1 - eps:
        nl, last = heapq.heappop(heap)
        tot += np.exp(-nl)
        d += 1
        if d >= cap:
            return None, tot
        for i in range(last, len(q)):
            heapq.heappush(heap, (nl + nlq[i], i))
    return d, tot


def main():
    plt.rcParams["font.family"] = "DejaVu Sans"
    N = 4000
    masses = (1.0, 0.1, 0.01)
    ls = (2, 4, 8, 16, 32, 64, 128)
    epss = (1e-2, 1e-4, 1e-6)
    rows = []
    print(f"[A,B] harmonic chain vacuum, ring of N = {N} sites; block of l sites")
    for m in masses:
        X, P = vacuum_correlators(N, m)
        for l in ls:
            q = block_modes(X, P, l)
            S = entropy_bits(q)
            ds = []
            for e in epss:
                d, w = smooth_dimension(q, e)
                ds.append(np.nan if d is None else d)
            rows.append([m, l, S] + [np.log2(x) if x == x else np.nan for x in ds])
            print(f"    m = {m:5.2f}, l = {l:4d}: S = {S:7.3f} bits; log2 d(eps) = "
                  + ", ".join(f"{np.log2(x):6.2f}" if x == x else "  >cap" for x in ds)
                  + f"  for eps = {', '.join(f'{e:.0e}' for e in epss)}")
    print("    single-particle entanglement energies e_k = ln(1/q_k) (two towers, one per block endpoint):")
    for m in masses:
        X, P = vacuum_correlators(N, m)
        for l in (64, 128):
            q = block_modes(X, P, l)
            e = -np.log(q[:6])
            act = [int(np.sum(q > 1e-3 * eps / l)) for eps in (1e-2, 1e-6)]
            print(f"      m = {m:5.2f}, l = {l:3d}: e_1..e_6 = {np.round(e, 2)}; tower spacing e_3 - e_1 = {e[2] - e[0]:.2f}; "
                  f"modes with q_k > 1e-3 eps/l: {act[0]} (eps = 1e-2), {act[1]} (eps = 1e-6)")
    R = np.array(rows)
    savedata("block_effective_dimension", ["mass", "l", "S_bits"] + [f"log2_d_eps_{e:.0e}" for e in epss],
             [R[:, i] for i in range(R.shape[1])],
             [f"harmonic chain vacuum, ring N={N}; d(eps) = smallest number of eigenvalues of rho_A with weight >= 1-eps"])

    # dependence on the precision at fixed block: exact d(eps) and the product (mode-truncation) bound
    X, P = vacuum_correlators(N, 0.01)
    q = block_modes(X, P, 64)
    S64 = entropy_bits(q)
    e_scan = np.logspace(-1, -12, 23)
    d_scan, b_scan, K_scan = [], [], []
    for e in e_scan:
        d, w = smooth_dimension(q, e)
        d_scan.append(np.nan if d is None else d)
        qq = q[q > 1e-3 * e / len(q)]
        K = len(qq)
        nk = np.ceil(np.log(K / e) / np.log(1 / qq))                 # keep occupations 0..n_k-1 in mode k
        b_scan.append(np.sum(np.log2(nk)))
        K_scan.append(K)
    d_scan = np.array(d_scan, float)
    b_scan = np.array(b_scan)
    lg = np.log2(d_scan)
    loc = np.diff(lg) / np.diff(np.log2(1 / e_scan))
    print(f"    m = 0.01, l = 64 (S = {S64:.3f} bits): exact log2 d(eps) and the product bound sum_k log2 n_k(eps)")
    for e, x, b, K in list(zip(e_scan, lg, b_scan, K_scan))[::2]:
        print(f"      eps = {e:.1e}: log2 d = {x:6.2f}, bound = {b:6.2f} ({K} modes)")
    print(f"    local slope d log2 d / d log2(1/eps): {loc[0]:.3f} at eps ~ 1e-1, {loc[len(loc)//2]:.3f} at eps ~ 1e-6, "
          f"{loc[-1]:.3f} at eps ~ 1e-12 (decreasing: growth slower than any power)")
    savedata("dimension_vs_eps", ["eps", "d", "log2_product_bound", "modes_kept"], [e_scan, d_scan, b_scan, K_scan],
             [f"m=0.01, l=64, S={S64:.4f} bits; bound: sum over kept modes of log2 ceil(ln(K/eps)/ln(1/q_k))"])

    fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    a = ax[0]
    for m in masses:
        sel = R[:, 0] == m
        a.semilogx(R[sel, 1], R[sel, 2], "o-", label=f"S(rho_A), m = {m}")
        a.semilogx(R[sel, 1], R[sel, 5], "s--", label=f"log2 d(1e-6), m = {m}")
    a.semilogx(ls, np.array(ls) * np.log2(20), ":", color="grey", label="l log2(20): 20 levels per site")
    a.set_ylim(0, 60)
    a.set_xlabel("block length l (sites)"); a.set_ylabel("bits")
    a.set_title("Effective dimension of a vacuum block: area law, not volume law")
    a.legend(fontsize=7, ncol=2)
    a = ax[1]
    a.semilogx(1 / e_scan, lg, "o-", label="exact log2 d(eps)")
    a.semilogx(1 / e_scan, b_scan, "s--", label="product bound over truncated modes")
    a.set_xlabel("1/eps (a Bell test resolves eps ~ 1/N)"); a.set_ylabel("bits")
    a.set_title("m = 0.01, l = 64: effective dimension vs precision")
    a.legend(fontsize=8)
    savefig(fig, "qft_effective_dimension")
    print(f"figures written to {FIG_DIR}/, data to {DATA_DIR}/")
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
