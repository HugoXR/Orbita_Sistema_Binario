import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

G = 6.674184E-11 # m^3 kg^-1 s^-2
M_1 = 1E30 # Star Mass 1 in Kg
M_2 = 1E30 # Star Mass 2 in Kg
M_p = 1E15 # Planet Mass in Kg

X = np.zeros(18)

X[0] = 1E11 # Position x of Star 1 in m
X[1] = 2E11 # Position y of Star 1 in m
X[2] = 0 # Position z of Star 1 in m
X[3] = -1E11 # Position x of Star 2 in m
X[4] = -2E11 # Position y of Star 2 in m
X[5] = 0 # Position z of Star 2 in m
X[6] = 10E11 # Position x of Planet in m
X[7] = 0 # Position y of Planet in m
X[8] = 0 # Position z of Planet in m

X[9] = 2E3
X[10] = 3E3
X[11] = 0
X[12] = 5E3
X[13] = 3E3
X[14] = 0
X[15] = 2E3
X[16] = 3E3
X[17] = 0

def dxdt(X, t):
    r12 = X[0:3] - X[3:6] # Position between Star 1 and 2
    rp1 = X[6:9] - X[0:3] # Position between Planet and Star 1
    rp2 = X[6:9] - X[3:6] # Position between Planet and Star 2
    mod_r12 = np.linalg.norm(r12) # Module of r12
    mod_rp1 = np.linalg.norm(rp1) # Module of rp1
    mod_rp2 = np.linalg.norm(rp2) # Module of rp2
    return np.hstack((np.array(X[9:]), -(G*M_2*r12)/(mod_r12**3),
                      -(G*M_1*r12)/(mod_r12**3), ((-(G*M_1*rp1)/(mod_rp1**3)) +
                     (-(G*M_2*rp2)/(mod_rp2**3)))))
