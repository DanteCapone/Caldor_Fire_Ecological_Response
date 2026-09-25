"""Lake Number from savalbuena/Lake-Tools, upstream commit f9cb02b."""

import numpy as np


def lakenumber(st, uSt, Ao, rho, zm, zt, zg):
    """Compute dimensionless Lake Number using the upstream equation."""
    g = 9.81
    st = st * Ao / g
    zt = zm - zt
    zg = zm - zg

    uSt = np.asarray(uSt, dtype=float).copy()
    uSt[uSt == 0] = np.nan
    return (g * st * (1 - zt / zm)) / (
        rho * uSt**2 * Ao ** (3 / 2) * (1 - zg / zm)
    )
