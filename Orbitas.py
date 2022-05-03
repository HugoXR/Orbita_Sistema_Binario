import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from astropy import constants
from astropy import units
import vpython
#Constantes iniciais
G = constants.G.value # Constante da Gravitacao
n=3 #Numero de Astros
d=3 #dimensoes do problema

M=np.zeros(n)
M[0] = (0.6897*units.Msun).decompose().value # Massa da primeira estrela em Massas Solares
M[1] = (0.20225*units.Msun).decompose().value # Massa da segunda estrela em Massas Solares
M[2] = (0.333*units.Mjup).decompose().value # Massa do planeta em Massas de Jupiter
mu = (M[0] * M[1])/(M[0] + M[1])

X=np.zeros(2*d*n)

a_S = (0.22431 * units.AU).decompose().value # Distancia inicial das estrela (semi-eixo maior) em Unidade Astronomica
a_P = (0.7048 * units.AU).decompose().value # Distancia inicial do planeta (semi-eixo maior) em Unidade Astronomica

P = 2*np.pi*np.sqrt((a_S**3)/(G*(M[0] + M[1]))) # Periodo da orbita das estrelas (entre si) pela terceira lei de kepler
P_p = 2*np.pi*np.sqrt((a_P**3)/(G*(M[0] + M[1]))) # Periodo da orbita do planeta em torno das estrelas pela terceira lei de kepler

r_1 = -a_S*(mu/M[0]) # Posicao relativa da primeira estrela
r_2 = a_S*(mu/M[1]) # Posicao relativa da segunda estrelas
r_p = a_P # Posicao do planeta

v_1 = 2*np.pi*r_1/P # Velocidade inicial da estrela 1
v_2 = 2*np.pi*r_2/P # Velocidade inicial da estrela 2
v_p = 2*np.pi*r_p/P_p # Velocidade inicial do planeta

X = np.array([r_1, 0, 0, r_2, 0, 0, r_p, 0, 0, 0, v_1, 0, 0, v_2, 0, 0,
              v_p, 0])


#Formula Geral (para o problema escrito em 3 dimensões)
def dxdt(X, t):
    n=len(X)//2
    a=np.zeros(n)
    for i in range(n//3):
        for j in range(n//3):
            if i!=j:
                a[3*i:3*(1+i)]+=-G*M[j]*(X[i*3:3*(1+i)]-X[j*3:3*(1+j)])/np.linalg.norm(X[i*3:3*(1+i)]-X[j*3:3*(1+j)])**3
    return np.hstack((X[n:],a))


#Grafico
t_orbit = np.linspace(0, P_p, 1000)

trajetoria = odeint(dxdt, X, t_orbit)

# plt.plot(trajetoria[:, 0], trajetoria[:, 1], label="Estrela Kepler-16A")
# plt.plot(trajetoria[:, 3], trajetoria[:, 4], label="Estrela Kepler-16B")
# plt.plot(trajetoria[:, 6], trajetoria[:, 7], label="Planeta Kepler-16b")
# plt.plot(trajetoria[-1, 0], trajetoria[-1, 1], 'o')
# plt.plot(trajetoria[-1, 3], trajetoria[-1, 4], 'o')
# plt.plot(trajetoria[-1, 6], trajetoria[-1, 7], 'o')
# plt.title("Sistema Circumbinário Kepler-16")
# plt.legend()
# plt.show()
# plt.savefig("Sistema Circumbinário Kepler-16.pdf")

