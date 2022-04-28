# Trajetoria orbita usando os dados do Kepler 16
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from astropy import constants
from astropy import units


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

#P = 2*np.pi*np.sqrt((a_S**3)/(G*(M_1 + M_2))) # Periodo da orbita das estrelas (entre si) pela terceira lei de kepler
#P_p = 2*np.pi*np.sqrt((a_P**3)/(G*(M_1 + M_2))) # Periodo da orbita do planeta em torno das estrelas pela terceira lei de kepler

P = (41.079220*units.d).decompose().value
P_p = (228.776*units.d).decompose().value

v_1 = 2*np.pi*r_1/P # Velocidade inicial da estrela 1
v_2 = 2*np.pi*r_2/P # Velocidade inicial da estrela 2
v_p = 2*np.pi*r_p/P_p # Velocidade inicial do planeta

#v_1 = -13E3
#v_2 = 45E3
#v_p = 33E3

X = np.array([r_1, 0, 0, r_2, 0, 0, r_p, 0, 0, 0, v_1, 0, 0, v_2, 0, 0,
              v_p, 0])
t_orbit = np.linspace(0, P_p, 1000)

trajetoria = odeint(dxdt, X, t_orbit)

# Plotando grafico
plt.plot(trajetoria[:,0], trajetoria[:,1], label="Estrela Kepler-16A")
plt.plot(trajetoria[:,3], trajetoria[:,4], label="Estrela Kepler-16B")
plt.plot(trajetoria[:,6], trajetoria[:,7], label="Planeta Kepler-16b")
plt.plot(trajetoria[-1,0], trajetoria[-1, 1], 'o')
plt.plot(trajetoria[-1,3], trajetoria[-1, 4], 'o')
plt.plot(trajetoria[-1,6], trajetoria[-1, 7], 'o')
plt.title("Sistema Circumbinário Kepler-16")
plt.legend()
plt.show()
plt.savefig("Sistema Circumbinário Kepler-16.pdf")
