"""
interval_scheduling.py
Greedy earliest-finish-time algorithm for maximum-cardinality interval scheduling.
Includes: DP verifier, random generator, and a quick self-test CLI.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Iterable, Optional
import bisect, random

@dataclass(frozen=True, order=True)
class Interval:
    start: float
    end: float
    def __post_init__(self):
        if self.end <= self.start:
            raise ValueError("Interval must have end > start")

def select_max_nonoverlap(intervals: Iterable[Tuple[float,float] | Interval]) -> List[Interval]:
    """Greedy: sort by end time; pick first compatible."""
    ivals = [i if isinstance(i, Interval) else Interval(*i) for i in intervals]
    ivals.sort(key=lambda x: x.end)
    selected: List[Interval] = []
    last_end = float("-inf")
    for iv in ivals:
        if iv.start >= last_end:
            selected.append(iv)
            last_end = iv.end
    return selected

# ---------- DP verifier (unit weights) ----------
def _preprocess_p(sorted_by_end: List[Interval]) -> List[int]:
    """p[i] = rightmost j < i with end[j] <= start[i]."""
    ends = [iv.end for iv in sorted_by_end]
    p = []
    for i, iv in enumerate(sorted_by_end):
        j = bisect.bisect_right(ends, iv.start) - 1
        p.append(j)
    return p

def max_cardinality_via_dp(intervals: Iterable[Tuple[float,float] | Interval]) -> int:
    ivals = [i if isinstance(i, Interval) else Interval(*i) for i in intervals]
    ivals.sort(key=lambda x: x.end)
    p = _preprocess_p(ivals)
    n = len(ivals)
    M = [0]*n
    for i in range(n):
        take = 1 + (M[p[i]] if p[i] >= 0 else 0)
        skip = M[i-1] if i-1 >= 0 else 0
        M[i] = max(take, skip)
    return M[-1] if n > 0 else 0

# ---------- Random generator ----------
def generate_random_intervals(n:int, day_length:float=24.0, max_duration:float=2.0,
                              seed:Optional[int]=None) -> List[Interval]:
    if seed is not None:
        random.seed(seed)
    out: List[Interval] = []
    for _ in range(n):
        s = random.random() * day_length
        dur = random.random() * max_duration + 1e-6
        e = min(day_length, s + dur)
        if e > s:
            out.append(Interval(s, e))
    return out

if __name__ == "__main__":
    # tiny smoke test
    ivals = generate_random_intervals(20, seed=1)
    g = len(select_max_nonoverlap(ivals))
    d = max_cardinality_via_dp(ivals)
    print("greedy =", g, " dp =", d, " match =", g == d)
