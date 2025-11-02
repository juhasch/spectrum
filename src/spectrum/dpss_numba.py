import numpy as np
from numba import njit


def _build_tridiagonal(N: int, NW: float):
    """Build the symmetric tridiagonal matrix for DPSS eigenproblem.

    This follows the construction used in the original C code (mydpss.c):
      diag[i] = -cos(2*pi*W) * (( (N-1)/2 - i )^2)
      offdiag[i] = -i * (N - i) / 2 for i>=1 (offdiag[0] unused)
    The full matrix T has diag on diagonal and offdiag on first sub/super diagonals.
    """
    W = float(NW) / float(N)
    twopi = 2.0 * np.pi
    cs = np.cos(twopi * W)

    diag = np.empty(N, dtype=np.float64)
    offdiag = np.empty(N, dtype=np.float64)
    for i in range(N):
        ai = float(i)
        diag[i] = -cs * (( (float(N) - 1.0) / 2.0 - ai) * ((float(N) - 1.0) / 2.0 - ai))
        offdiag[i] = -ai * (float(N) - ai) / 2.0
    # Build dense symmetric matrix from tridiagonal components
    T = np.zeros((N, N), dtype=np.float64)
    T[np.arange(N), np.arange(N)] = diag
    for i in range(1, N):
        T[i, i - 1] = offdiag[i]
        T[i - 1, i] = offdiag[i]
    return T


@njit(cache=True)
def _normalize_and_pack(evecs: np.ndarray) -> tuple:
    """Normalize eigenvectors to RMS=1 and compute tapsum; pack to 1D blocks.

    evecs: shape (N, K), column k is the k-th taper before normalization.
    Returns (tapers_flat, tapsum) where tapers_flat is 1D with K contiguous blocks
    of length N (same layout as original C code).
    """
    N, K = evecs.shape
    tapsum = np.zeros(K, dtype=np.float64)
    tapers_flat = np.empty(N * K, dtype=np.float64)
    for k in range(K):
        # Copy column
        base = k * N
        tapsq = 0.0
        for i in range(N):
            val = evecs[i, k]
            tapers_flat[base + i] = val
            tapsum[k] += val
            tapsq += val * val
        # Normalize to RMS=1 (sum(x^2)/N == 1)
        aa = (tapsq / float(N)) ** 0.5
        if aa == 0.0:
            aa = 1.0
        tapsum[k] = tapsum[k] / aa
        for i in range(N):
            tapers_flat[base + i] = tapers_flat[base + i] / aa
    return tapers_flat, tapsum


def dpss_tapers(N: int, NW: float, Kmax: int):
    """Compute DPSS tapers using a Python/NumPy core and numba-accelerated packing.

    Parameters
    ----------
    N : int
        Window length.
    NW : float
        Time-half bandwidth product.
    Kmax : int
        Number of tapers to return (first Kmax tapers corresponding to the
        largest eigenvalues).

    Returns
    -------
    tapers : ndarray, shape (N, Kmax)
        Normalized tapers with L2 norm sqrt(N) (RMS=1), matching the C routine
        prior to the additional 1/sqrt(N) scaling applied in mtm.dpss.
    tapsum : ndarray, shape (Kmax,)
        Sum of each taper (used for sign convention handling).
    """
    # Build the tridiagonal operator and solve for eigenpairs
    T = _build_tridiagonal(N, NW)

    # Solve dense symmetric eigendecomposition (ascending eigenvalues)
    # We select the Kmax eigenvectors with the largest eigenvalues
    w, v = np.linalg.eigh(T)
    idx = np.argsort(w)[-Kmax:]
    # Ensure ascending order across the selected set to match expectations
    idx.sort()
    evecs = v[:, idx]

    # The eigenvectors from eigh have unit L2 norm. Normalize to RMS=1 (L2=sqrt(N))
    tapers_flat, tapsum = _normalize_and_pack(evecs)

    # Return as (N, K) matrix for easier downstream handling
    tapers = tapers_flat.reshape(Kmax, N).T
    return tapers, tapsum


