"""Schmidt stability from savalbuena/Lake-Tools, commit f9cb02b."""

import numpy as np


def schmidtstability(z, rho, A, A0, z_bar_in, z_bar_input, g=9.806):
    """Compute Schmidt stability in J m-2 using the upstream equation."""
    z = np.asarray(z)
    rho = np.asarray(rho)
    A = np.asarray(A)
    volume = np.trapezoid(A, z)
    zcv = z_bar_input if z_bar_in else np.trapezoid(z * A, z) / volume
    integrand = (z - zcv) * rho * A
    return g / A0 * np.trapezoid(integrand, z)
