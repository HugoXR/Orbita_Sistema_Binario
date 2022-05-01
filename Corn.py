import numpy as np



G=1
X = np.array([1, 3, 0, 2, 0, 2, 3, 0, 9, 0, 4, 0, 0, 5, 0, 0,
              6, 0])
M=np.zeros(len(X)//6)
M[0]=M_1 = 1
M[1]=M_2 = 1
M[2]=M_P = 1


def dxdt(X, t):
    n=len(X)//2
    a=np.zeros(n)
    for i in range(n//3):
        for j in range(n//3):
            if i!=j:
                a[3*i:3*(1+i)]+=-G*M[j]*(X[i*3:3*(1+i)]-X[j*3:3*(1+j)])/np.linalg.norm(X[i*3:3*(1+i)]-X[j*3:3*(1+j)])**3
    return np.hstack((X[n:],a))

def Dxdt(X, t):
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

a=dxdt(X,1)
b=Dxdt(X,1)
print(a-b)