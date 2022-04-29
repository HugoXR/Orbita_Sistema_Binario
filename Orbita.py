# Trajetoria orbita usando os dados do Kepler 16
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from astropy import constants
from astropy import units

import vpython

def dxdt(X, t):
    r12 = X[3:6] - X[0:3] # Posicao entre a primeira e segunda estrela
    rp1 = X[6:9] - X[0:3] # Posicao entre a primeira estrela e o planeta
    rp2 = X[6:9] - X[3:6] # Posicao entre a segunda estrela e o planeta
    mod_r12 = np.linalg.norm(r12) # Modulo da posicao r12
    mod_rp1 = np.linalg.norm(rp1) # Modulo da posicao rp1
    mod_rp2 = np.linalg.norm(rp2) # Modulo da posicao rp2
    return np.hstack((np.array(X[9:]), (G*M_2*r12)/(mod_r12**3) +
                      (G*M_P*rp1/(mod_rp1**3)),
                      -(G*M_1*r12)/(mod_r12**3) + (G*M_P*rp2)/(mod_rp2**3), ((-(G*M_1*rp1)/(mod_rp1**3)) +
                     (-(G*M_2*rp2)/(mod_rp2**3)))))


M_1 = (0.6897*units.Msun).decompose().value # Massa da primeira estrela em Massas Solares
M_2 = (0.20225*units.Msun).decompose().value # Massa da segunda estrela em Massas Solares
M_P = (0.333*units.Mjup).decompose().value # Massa do planeta em Massas de Jupiter

mu = (M_1 * M_2)/(M_1 + M_2) # Massa reduzida das duas estrelas
G = constants.G.value # Constante da Gravitacao

a_S = (0.22431 * units.AU).decompose().value # Distancia inicial das estrela (semi-eixo maior) em Unidade Astronomica
a_P = (0.7048 * units.AU).decompose().value # Distancia inicial do planeta (semi-eixo maior) em Unidade Astronomica

r_1 = -a_S*(mu/M_1) # Posicao relativa da primeira estrela
r_2 = a_S*(mu/M_2) # Posicao relativa da segunda estrelas
r_p = a_P # Posicao do planeta

P = 2*np.pi*np.sqrt((a_S**3)/(G*(M_1 + M_2))) # Periodo da orbita das estrelas (entre si) pela terceira lei de kepler
P_p = 2*np.pi*np.sqrt((a_P**3)/(G*(M_1 + M_2))) # Periodo da orbita do planeta em torno das estrelas pela terceira lei de kepler

#P = (41.079220*units.d).decompose().value # Periodo da orbita das estrelas(entre si) em dias
#P_p = (228.776*units.d).decompose().value # Periodo da orbita do planeta em torno das estrelas em dias

v_1 = 2*np.pi*r_1/P # Velocidade inicial da estrela 1
v_2 = 2*np.pi*r_2/P # Velocidade inicial da estrela 2
v_p = 2*np.pi*r_p/P_p # Velocidade inicial do planeta

#v_1 = np.sqrt(2*G*M_1/abs(r_1))
#v_2 = -np.sqrt(2*G*M_2/abs(r_2))
#v_p = np.sqrt(2*G*M_1/abs(r_1)) + np.sqrt(2*G*M_2/abs(r_2))

X = np.array([r_1, 0, 0, r_2, 0, 0, r_p, 0, 0, 0, v_1, 0, 0, v_2, 0, 0,
              v_p, 0])
t_orbit = np.linspace(0, P_p, 1000)

trajetoria = odeint(dxdt, X, t_orbit)

#Simulacao com vpython
star1 = vpython.sphere(pos=vpython.vector(r_1, 0, 0), radius=(0.6897*units.R_sun).decompose().value, color=vpython.color.yellow, make_trail=True)
star2 = vpython.sphere(pos=vpython.vector(r_2, 0, 0), radius=(0.22623*units.R_sun).decompose().value, color=vpython.color.cyan, make_trail=True)
planet = vpython.sphere(pos=vpython.vector(r_p, 0, 0),
                        radius=(0.7538*units.R_jupiter).decompose().value,
                        color=vpython.color.red, make_trail=True)

for ponto in trajetoria:
    vpython.rate(100)
    star1.pos = vpython.vector(ponto[0], ponto[1], ponto[2])
    star2.pos = vpython.vector(ponto[3], ponto[4], ponto[5])
    planet.pos = vpython.vector(ponto[6], ponto[7], ponto[8])

# Plotando grafico
plt.plot(trajetoria[:, 0], trajetoria[:, 1], label="Estrela Kepler-16A")
plt.plot(trajetoria[:, 3], trajetoria[:, 4], label="Estrela Kepler-16B")
plt.plot(trajetoria[:, 6], trajetoria[:, 7], label="Planeta Kepler-16b")
plt.plot(trajetoria[-1, 0], trajetoria[-1, 1], 'o')
plt.plot(trajetoria[-1, 3], trajetoria[-1, 4], 'o')
plt.plot(trajetoria[-1, 6], trajetoria[-1, 7], 'o')
plt.title("Sistema Circumbinário Kepler-16")
plt.legend()
plt.show()
plt.savefig("Sistema Circumbinário Kepler-16.pdf")
