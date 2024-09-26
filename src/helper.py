
import numpy as np
from scipy.optimize import minimize

def prandtl_meyer_from_mach(Mach, gamma=1.4):

    if Mach <= 1:
        return 0

    alpha = np.sqrt((gamma+1) / (gamma-1))
    beta = np.sqrt(Mach**2 - 1)

    return alpha * np.arctan(beta/alpha) - np.arctan(beta)

def mach_from_prandtl_meyer(pm_radians, gamma=1.4):
    res = minimize(lambda mach, gamma: abs(prandtl_meyer_from_mach(mach, gamma) - pm_radians),
                   x0=2.0,
                   args=(gamma))

    return res.x[0]

if __name__ == "__main__":
    for M in np.linspace(1.01, 2.5, 100):
        print('test for M={}'.format(M))
        print(prandtl_meyer_from_mach(M))
        # print(mach_from_prandtl_meyer(prandtl_meyer_from_mach(M)))
        print('error: {}%'.format((mach_from_prandtl_meyer(prandtl_meyer_from_mach(M)) - M)/M * 100))