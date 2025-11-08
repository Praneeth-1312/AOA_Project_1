"""
experiments.py
Generates all artifacts for the paper:
- Greedy vs DP verification table
- D&C vs brute-force verification and crossover
- Timing with repeats + error bars
- log-log linear regression vs n log n (slope & R^2)
- Plots (PNG), CSVs, JSON summary, LaTeX table & snippet files
Outputs go to ./artifacts (change with --out).
"""

from __future__ import annotations
import argparse, csv, json, math, os, platform, random, statistics as stats, sys, time
from pathlib import Path
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt

from interval_scheduling import (
    select_max_nonoverlap, generate_random_intervals, max_cardinality_via_dp
)
from closest_pair import (
    closest_pair_distance, brute_force, generate_random_points
)

# -------------------- helpers --------------------
def ensure_outdir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True); return p

def loglog_regression(xs: Iterable[float], ys: Iterable[float]) -> Tuple[float, float]:
    """Return slope and R^2 for regression of log2(ys) ~ a + b*log2(xs).
    (Simple least squares on logs.)"""
    lx = [math.log2(x) for x in xs]
    ly = [math.log2(y) for y in ys]
    n = len(lx)
    mx, my = sum(lx)/n, sum(ly)/n
    num = sum((lx[i]-mx)*(ly[i]-my) for i in range(n))
    den = sum((lx[i]-mx)**2 for i in range(n))
    b = num/den if den > 0 else float("nan")
    a = my - b*mx
    # R^2
    ss_tot = sum((ly[i]-my)**2 for i in range(n))
    ss_res = sum((ly[i] - (a + b*lx[i]))**2 for i in range(n))
    r2 = 1 - (ss_res/ss_tot if ss_tot > 0 else 0.0)
    return b, r2

def write_latex_table(outfile: Path, caption: str, header: List[str], rows: List[List[str]], label: str):
    with outfile.open("w", encoding="utf-8") as f:
        f.write("\\begin{table}[h]\n\\centering\n")
        f.write("\\begin{tabular}{%s}\n" % ("|" + "|".join(["l"]*len(header)) + "|"))
        f.write("\\hline\n")
        f.write(" & ".join(header) + " \\\\\n\\hline\n")
        for r in rows:
            f.write(" & ".join(map(str, r)) + " \\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        f.write("\\caption{%s}\\label{%s}\n" % (caption, label))
        f.write("\\end{table}\n")

def env_report() -> dict:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "matplotlib": getattr(plt, "__version__", "unknown")
    }

# -------------------- experiments --------------------
def greedy_timings(ns: Iterable[int], repeats: int, seed: int) -> Tuple[List[Tuple[int,float,float,float]], dict]:
    """Return list of (n, mean, std, selected_mean) and sanity dict."""
    rng = random.Random(seed)
    rows = []
    # sanity: greedy == DP on small instances
    sanity_trials = 100
    sanity_ok = 0
    for _ in range(sanity_trials):
        iv = generate_random_intervals(40, seed=rng.randrange(1<<30))
        if len(select_max_nonoverlap(iv)) == max_cardinality_via_dp(iv):
            sanity_ok += 1
    sanity = {"trials": sanity_trials, "passed": sanity_ok, "pass_rate": sanity_ok/sanity_trials}

    for n in ns:
        times = []
        selected_counts = []
        for _ in range(repeats):
            ivs = generate_random_intervals(n, seed=rng.randrange(1<<30))
            t0 = time.perf_counter()
            sel = select_max_nonoverlap(ivs)
            t1 = time.perf_counter()
            times.append(t1 - t0)
            selected_counts.append(len(sel))
        rows.append((n, stats.mean(times), (stats.stdev(times) if len(times)>1 else 0.0),
                     stats.mean(selected_counts)))
    return rows, sanity

def dc_timings(ns: Iterable[int], repeats: int, seed: int) -> List[Tuple[int,float,float]]:
    rng = random.Random(seed)
    rows = []
    for n in ns:
        times = []
        for _ in range(repeats):
            pts = generate_random_points(n, seed=rng.randrange(1<<30))
            t0 = time.perf_counter()
            _ = closest_pair_distance(pts)
            t1 = time.perf_counter()
            times.append(t1 - t0)
        rows.append((n, stats.mean(times), (stats.stdev(times) if len(times)>1 else 0.0)))
    return rows

def dc_crossover(ns_small: Iterable[int], repeats: int, seed: int) -> List[Tuple[int,float,float,float,float]]:
    """Return list of (n, mean_dc, std_dc, mean_bf, std_bf)."""
    rng = random.Random(seed)
    rows = []
    for n in ns_small:
        tdc, tbf = [], []
        for _ in range(repeats):
            pts = generate_random_points(n, seed=rng.randrange(1<<30))
            t0 = time.perf_counter(); _ = closest_pair_distance(pts); t1 = time.perf_counter()
            t2 = time.perf_counter(); _ = brute_force(pts);          t3 = time.perf_counter()
            tdc.append(t1 - t0); tbf.append(t3 - t2)
        rows.append((n,
                     stats.mean(tdc), (stats.stdev(tdc) if len(tdc)>1 else 0.0),
                     stats.mean(tbf), (stats.stdev(tbf) if len(tbf)>1 else 0.0)))
    return rows

# -------------------- plotting --------------------
def plot_with_errorbars(rows, title, outfile: Path, ylabel="time (s)", nlogn=True):
    # rows: list of tuples where first item is n, second is mean time, third std
    ns  = [r[0] for r in rows]
    mu  = [r[1] for r in rows]
    std = [r[2] for r in rows]
    plt.figure(figsize=(7,5))
    plt.errorbar(ns, mu, yerr=std, marker="o", capsize=4, label="empirical (mean±sd)")
    if nlogn and mu[0] > 0:
        scale = mu[0] / (ns[0] * math.log2(ns[0]))
        ref = [scale * n * math.log2(n) for n in ns]
        plt.plot(ns, ref, "--", label="scaled $n\\log n$")
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("n"); plt.ylabel(ylabel); plt.title(title)
    plt.grid(True, which="both", ls=":"); plt.legend(); plt.tight_layout()
    plt.savefig(outfile); plt.close()

def plot_crossover(rows, outfile: Path):
    ns   = [r[0] for r in rows]
    mu_d = [r[1] for r in rows]
    sd_d = [r[2] for r in rows]
    mu_b = [r[3] for r in rows]
    sd_b = [r[4] for r in rows]
    plt.figure(figsize=(7,5))
    plt.errorbar(ns, mu_d, yerr=sd_d, marker="o", capsize=4, label="D&C")
    plt.errorbar(ns, mu_b, yerr=sd_b, marker="s", capsize=4, label="Brute force")
    plt.xscale("log"); plt.yscale("log")
    plt.xlabel("n"); plt.ylabel("time (s)"); plt.title("Closest Pair: D&C vs Brute Force")
    plt.grid(True, which="both", ls=":"); plt.legend(); plt.tight_layout()
    plt.savefig(outfile); plt.close()

# -------------------- CSV writers --------------------
def write_rows_csv(path: Path, header: List[str], rows: List[Tuple]):
    with path.open("w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

# -------------------- LaTeX snippet generators --------------------
def write_verification_tables(outdir: Path, sanity_g: dict, cross_rows: List[Tuple[int,float,float,float,float]]):
    # Greedy sanity table
    write_latex_table(
        outdir/"greedy_sanity.tex",
        "Greedy vs DP verification on 40-interval random instances.",
        ["Trials","Passed","Pass rate"],
        [[sanity_g["trials"], sanity_g["passed"], f'{100*sanity_g["pass_rate"]:.1f}\\%']],
        "tab:greedy-sanity"
    )
    # Crossover table (first 5 rows)
    head = ["n","D&C mean (s)","D&C sd","Brute mean (s)","Brute sd"]
    rows = []
    for r in cross_rows:
        rows.append([r[0], f"{r[1]:.6g}", f"{r[2]:.3g}", f"{r[3]:.6g}", f"{r[4]:.3g}"])
    write_latex_table(
        outdir/"dc_crossover_table.tex",
        "Closest-pair D\\&C vs brute-force crossover (averaged).",
        head, rows, "tab:dc-crossover"
    )

def write_fig_include_snippet(outdir: Path):
    with (outdir/"figures_snippet.tex").open("w") as f:
        f.write(
r"""\begin{figure}[h]
\centering
\includegraphics[width=0.48\textwidth]{greedy_runtime.png}
\includegraphics[width=0.48\textwidth]{dc_runtime.png}
\caption{Empirical runtimes with $n\log n$ overlays (mean$\pm$sd, log--log).}
\label{fig:runtime-both}
\end{figure}

\begin{figure}[h]
\centering
\includegraphics[width=0.6\textwidth]{dc_crossover.png}
\caption{Closest pair: divide-and-conquer vs brute-force crossover (log--log).}
\label{fig:dc-crossover}
\end{figure}
""")

def write_regression_snippet(outdir: Path, title: str, slope: float, r2: float, label: str):
    with (outdir / f"{label}_regression.tex").open("w") as f:
        f.write(
            f"\\noindent\\textbf{{{title}.}} log--log slope $= {slope:.3f}$, $R^2 = {r2:.4f}$.\\\\\n"
        )

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("artifacts"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--repeats", type=int, default=5, help="Repeats per n for timing (error bars).")
    ap.add_argument("--greedy_ns", type=int, nargs="*", default=[1000,2000,4000,8000,16000,32000])
    ap.add_argument("--dc_ns", type=int, nargs="*", default=[256,512,1024,2048,4096,8192])
    ap.add_argument("--cross_ns", type=int, nargs="*", default=[32,64,128,256])
    args = ap.parse_args()

    outdir = ensure_outdir(args.out)

    # --- GREEDY ---
    g_rows, g_sanity = greedy_timings(args.greedy_ns, args.repeats, args.seed)
    write_rows_csv(outdir/"greedy_runtimes.csv",
                   ["n","time_mean_sec","time_sd_sec","selected_mean"], g_rows)
    plot_with_errorbars(g_rows, "Interval Scheduling (Greedy) Runtime", outdir/"greedy_runtime.png")
    # regression vs n log n
    ns = [r[0] for r in g_rows]; ts = [r[1] for r in g_rows]
    slope_g, r2_g = loglog_regression([n*math.log2(n) for n in ns], ts)
    write_regression_snippet(outdir, "Greedy runtime vs $n\\log n$", slope_g, r2_g, "greedy")

    # --- D&C ---
    d_rows = dc_timings(args.dc_ns, args.repeats, args.seed)
    write_rows_csv(outdir/"dc_runtimes.csv", ["n","time_mean_sec","time_sd_sec"], d_rows)
    plot_with_errorbars(d_rows, "Closest Pair (D&C) Runtime", outdir/"dc_runtime.png")
    slope_d, r2_d = loglog_regression([n*math.log2(n) for n in [r[0] for r in d_rows]], [r[1] for r in d_rows])
    write_regression_snippet(outdir, "D\\&C runtime vs $n\\log n$", slope_d, r2_d, "dc")

    # --- Crossover D&C vs Brute ---
    cross_rows = dc_crossover(args.cross_ns, args.repeats, args.seed)
    write_rows_csv(outdir/"dc_crossover.csv",
                   ["n","dc_mean_sec","dc_sd_sec","brute_mean_sec","brute_sd_sec"], cross_rows)
    plot_crossover(cross_rows, outdir/"dc_crossover.png")

    # --- Tables & snippets ---
    write_verification_tables(outdir, g_sanity, cross_rows)
    write_fig_include_snippet(outdir)

    # --- Env & summary ---
    summary = {
        "seed": args.seed,
        "repeats": args.repeats,
        "env": env_report(),
        "greedy_sanity": g_sanity,
        "regression": {
            "greedy": {"slope_loglog_vs_nlogn": slope_g, "R2": r2_g},
            "dc": {"slope_loglog_vs_nlogn": slope_d, "R2": r2_d}
        },
        "artifacts": {
            "greedy_csv": str(outdir/"greedy_runtimes.csv"),
            "dc_csv": str(outdir/"dc_runtimes.csv"),
            "cross_csv": str(outdir/"dc_crossover.csv"),
            "greedy_plot": str(outdir/"greedy_runtime.png"),
            "dc_plot": str(outdir/"dc_runtime.png"),
            "cross_plot": str(outdir/"dc_crossover.png"),
            "greedy_regression_tex": str(outdir/"greedy_regression.tex"),
            "dc_regression_tex": str(outdir/"dc_regression.tex"),
            "greedy_sanity_tex": str(outdir/"greedy_sanity.tex"),
            "dc_crossover_table_tex": str(outdir/"dc_crossover_table.tex"),
            "figures_snippet_tex": str(outdir/"figures_snippet.tex"),
        }
    }
    with (outdir/"summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    print("Artifacts written to:", outdir.resolve())
    for k, v in summary["artifacts"].items():
        print(f" - {k}: {v}")
    print("Greedy sanity:", g_sanity)
    print("Env:", summary["env"])

if __name__ == "__main__":
    main()
