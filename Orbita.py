import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from astropy import constants
from astropy import units


def dxdt(X, t):
    r12 = X[0:3] - X[3:6] # Position between Star 1 and 2
    rp1 = X[6:9] - X[0:3] # Position between Planet and Star 1
    rp2 = X[6:9] - X[3:6] # Position between Planet and Star 2
    mod_r12 = np.linalg.norm(r12) # Module of r12
    mod_rp1 = np.linalg.norm(rp1) # Module of rp1
    mod_rp2 = np.linalg.norm(rp2) # Module of rp2
    return np.hstack((np.array(X[9:]), -(G*M_2*r12)/(mod_r12**3),
                      (G*M_1*r12)/(mod_r12**3), ((-(G*M_1*rp1)/(mod_rp1**3)) +
                     (-(G*M_2*rp2)/(mod_rp2**3)))))


M_1 = (1*units.Msun).decompose().value # Massa da primeira estrela = 1 Massa Solar
M_2 = (5*units.Msun).decompose().value # Massa da segunda estrela = 5 Massas Solares
M_P = (10*units.Mjup).decompose().value # Massa do planeta = 10 Massas de Jupiter

mu = (M_1 * M_2)/(M_1 + M_2) # Massa reduzida das duas estrelas
G = constants.G.value # Constante da Gravitacao

a_S = (1 * units.AU).decompose().value # Distancia inicial das estrela (semi-eixo maior) = 1 Unidade Astronomica
a_P = (2 * units.AU).decompose().value # Distancia inicial do planeta (semi-eixo maior) = 1 Unidade Astronomica

r_1 = a_S*(mu/M_1) # Posicao relativa da primeira estrela
r_2 = -a_S*(mu/M_2) # Posicao relativa da segunda estrelas
r_p = a_P # Posicao do planeta

P = 2*np.pi*np.sqrt((a_S**3)/(G*(M_1 + M_2))) # Periodo da orbita das estrelas (entre si) pela terceira lei de kepler
P_p = 2*np.pi*np.sqrt((a_P**3)/(G*(M_1 + M_2))) # Periodo da orbita do planeta em torno das estrelas pela terceira lei de kepler

v_1 = 2*np.pi*r_1/P # Velocidade inicial da estrela 1
v_2 = 2*np.pi*r_2/P # Velocidade inicial da estrela 2
v_p = 2*np.pi*r_p/P_p # Velocidade inicial do planeta

X = np.array([r_1, 0, 0, r_2, 0, 0, r_p, 0, 0, v_1, 0, 0, v_2, 0, 0, v_p, 0, 0])
t_orbit = np.linspace(0, P_p, 5)

trajetoria = odeint(dxdt, X, t_orbit)


