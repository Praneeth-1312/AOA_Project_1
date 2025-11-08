"""
closest_pair.py
Divide-and-conquer closest pair (2D).
Includes: brute-force verifier, random generator, and a quick self-test CLI.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Iterable, Optional
import math, random

@dataclass(frozen=True)
class Point:
    x: float
    y: float

def dist(a: Point, b: Point) -> float:
    return math.hypot(a.x - b.x, a.y - b.y)

def brute_force(points: List[Point]) -> float:
    best = float("inf")
    n = len(points)
    for i in range(n):
        for j in range(i+1, n):
            d = dist(points[i], points[j])
            if d < best:
                best = d
    return best

def closest_pair_distance(points: Iterable[Tuple[float,float] | Point]) -> float:
    pts = [p if isinstance(p, Point) else Point(*p) for p in points]
    px = sorted(pts, key=lambda p: p.x)
    py = sorted(pts, key=lambda p: p.y)
    return _closest_rec(px, py)

def _closest_rec(px: List[Point], py: List[Point]) -> float:
    n = len(px)
    if n <= 3:
        return brute_force(px)
    mid = n // 2
    Qx, Rx = px[:mid], px[mid:]
    midx = px[mid].x
    Qy, Ry = [], []
    for p in py:
        (Qy if p.x <= midx else Ry).append(p)
    dl = _closest_rec(Qx, Qy)
    dr = _closest_rec(Rx, Ry)
    d  = dl if dl < dr else dr
    strip = [p for p in py if abs(p.x - midx) < d]
    best = d
    for i in range(len(strip)):
        for j in range(i+1, min(i+8, len(strip))):
            dij = dist(strip[i], strip[j])
            if dij < best:
                best = dij
    return best

def generate_random_points(n:int, seed:Optional[int]=None) -> List[Point]:
    if seed is not None:
        random.seed(seed)
    return [Point(random.random(), random.random()) for _ in range(n)]

if __name__ == "__main__":
    pts = generate_random_points(30, seed=2)
    dc = closest_pair_distance(pts)
    bf = brute_force(pts)
    print("dc =", dc, " bf =", bf, " match =", abs(dc-bf) < 1e-12)
