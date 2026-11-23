"""Reference (numpy) implementation of the serial-compartment sinusoid transport for ONE field.
Backward-Euler step for tissue voxels (50x50 lattice, conductance ge between neighbours, conductance gw to its blood node) and K serial blood nodes
(6 sectors x 8 segments, flow Fs through each sector chain, inlet concentration Cin at the first node of every sector). Amounts in pmol, volumes nL, conc mM.
The step is verified against a dense solve of the full (N+K) system; mass ledger is certified. Used to validate the C++ operator."""
import math
import numpy as np
from scipy.linalg import solve_banded, lu_factor, lu_solve
N, NX, K, SECT, SEG = 2500, 50, 48, 6, 8
VOX = 0.008


def node_map(radius=490.0):
    """voxel -> node (sector*8+segment) with the candidate rule of lobule_vascular_geometry_bridge_2026_11_17/geometry.py (outer ring = segment 0, centre = segment 7)."""
    idx = np.arange(N); ix, iy = idx % NX, idx // NX
    x = -500 + 20 * (ix + .5); y = -500 + 20 * (iy + .5)
    r = np.hypot(x, y); theta = np.mod(np.arctan2(y, x) + math.pi / 6, 2 * math.pi)
    sector = np.floor(theta / (math.pi / 3)).astype(int)
    segment = 7 - np.searchsorted(np.arange(1, 8) * radius / 8, r, side="left")
    return sector * 8 + segment


class SerialTransport:
    def __init__(self, dt, ge, gw, Fs, Vnode, node_of):
        self.dt, self.ge, self.gw, self.Fs, self.Vnode, self.node_of = dt, ge, gw, Fs, Vnode, node_of
        self.m = np.bincount(node_of, minlength=K).astype(float)
        assert (self.m > 0).all()
        # tissue matrix A = V/dt + L + gw I (banded, half bandwidth 50)
        A = np.zeros((N, N)); d = np.full(N, VOX / dt + gw)
        for y in range(NX):
            for x in range(NX):
                i = x + NX * y
                for j in ((i + 1) if x < NX - 1 else None, (i + NX) if y < NX - 1 else None):
                    if j is not None: A[i, j] -= ge; A[j, i] -= ge; d[i] += ge; d[j] += ge
        A[np.arange(N), np.arange(N)] = d; self.A = A
        # band storage for solve_banded (u=l=50)
        ab = np.zeros((101, N))
        for off in range(-50, 51):
            diag = np.diagonal(A, offset=off); ab[50 - off, max(0, off): max(0, off) + len(diag)] = diag if off >= 0 else np.r_[np.zeros(0), diag]
        # simpler: build via scipy.sparse for exactness
        import scipy.sparse as sp, scipy.sparse.linalg as spl
        self.As = sp.csc_matrix(A); self.lu = spl.splu(self.As)
        E = np.zeros((N, K)); E[np.arange(N), node_of] = 1.0; self.E = E
        self.W = self.lu.solve(E)                                           # A^-1 E
        B = np.zeros((K, K)); first = (np.arange(K) % SEG) == 0
        for k in range(K):
            B[k, k] = Vnode / dt + Fs + gw * self.m[k]
            if not first[k]: B[k, k - 1] = -Fs
        self.S = B - gw * gw * (E.T @ self.W); self.Slu = lu_factor(self.S); self.B = B; self.first = first

    def step(self, c_old, b_old, Cin):
        dt, gw, Fs, V = self.dt, self.gw, self.Fs, self.Vnode
        rc = VOX / dt * c_old; rb = V / dt * b_old + np.where(self.first, Fs * Cin, 0.0)
        y = self.lu.solve(rc); r = rb + gw * (self.E.T @ y); b = lu_solve(self.Slu, r); c = y + gw * (self.W @ b)
        return c, b

    def dense_step(self, c_old, b_old, Cin):
        n = N + K; M = np.zeros((n, n)); rhs = np.zeros(n)
        M[:N, :N] = self.A; M[:N, N:] = -self.gw * self.E; M[N:, :N] = -self.gw * self.E.T; M[N:, N:] = self.B
        rhs[:N] = VOX / self.dt * c_old; rhs[N:] = self.Vnode / self.dt * b_old + np.where(self.first, self.Fs * Cin, 0.0)
        x = np.linalg.solve(M, rhs); return x[:N], x[N:]

    def ledger(self, c_old, b_old, c, b, Cin):
        """pmol: change in stored amount equals inflow - outflow over the step (outflow evaluated at the new last-node concentration)."""
        before = VOX * c_old.sum() + self.Vnode * b_old.sum(); after = VOX * c.sum() + self.Vnode * b.sum()
        inflow = self.dt * self.Fs * Cin * SECT; last = (np.arange(K) % SEG) == SEG - 1; outflow = self.dt * self.Fs * b[last].sum()
        return (after - before) - inflow + outflow


if __name__ == "__main__":
    nm = node_map(); print("voxels per node: min %d max %d, per segment (all sectors):" % (np.bincount(nm).min(), np.bincount(nm).max()), np.bincount(nm % 8, minlength=8))
    for dt, gw in ((0.05, 3.504), (0.05, 0.00350), (0.5, 3.504)):
        T = SerialTransport(dt, 3.504, gw, 8.06 / 6, 2.0 / 48, nm); rng = np.random.default_rng(1)
        c0 = rng.random(N) * 0.5; b0 = rng.random(K) * 0.5
        c, b = T.step(c0, b0, 0.8); cd, bd = T.dense_step(c0, b0, 0.8)
        print("dt=%s gw=%s  max|schur-dense| tissue %.2e node %.2e  ledger %.2e  min c %.3e min b %.3e" % (dt, gw, abs(c - cd).max(), abs(b - bd).max(), T.ledger(c0, b0, c, b, 0.8), c.min(), b.min()))
