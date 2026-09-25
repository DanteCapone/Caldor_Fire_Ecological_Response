"""Water-side shear velocity from savalbuena/Lake-Tools, commit f9cb02b."""

import numpy as np

RHO_AIR = 1.2
VON_KARMAN = 0.4


def u_star(wspd, windHeight, EpiRho):
    """Calculate water-side friction velocity from unaveraged wind speed."""
    wspd = np.asarray(wspd, dtype=float).copy()
    cd = 0.0013 * np.ones(len(wspd))
    cd[wspd < 5] = 0.001
    if windHeight != 10:
        wspd = wspd / (
            1 - np.sqrt(cd) / VON_KARMAN * np.log(10 / windHeight)
        )
    tau = cd * RHO_AIR * wspd**2
    return np.sqrt(tau / EpiRho)
