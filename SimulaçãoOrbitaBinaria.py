import vpython as vp 

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 12:37:05 2022

@author: biell
"""


#Aqui vai ficar as variaveis que vamos poder mudar
# Distancia inicial entre o planeta e sua estrela
distancia_inicial_planeta   = 20


# Velocidade inicial do planeta para com a estrela em AU (unidade astronômica) por ano.
# Unidades métricas (SI) 1 AU = 1.495978707×10e+11 m
vel_inicial_planeta  =   2 

#Massa do planeta, como um multiplo da massa do sol 
massa_planeta    =  0.001  

#Massa da Estrela 1 
massa_estrela1        = 1.3


#Massa da Estrela 2
massa_estrela2        = .7

#Quanto as estrelas estão distantes
separacao_estrelas    = 2

# Velocidade inicial das estrelas 
vel_inicialestrelas  =  8


#Quantos segundos, 1 ano simulado deve levar 
segundos_por_anos = 0.001



#inclinação do plano 
inclinacao = 30


#vetores 
#vetores posição estrelas 
pos_estrela1 = vp.vector(separacao_estrelas * massa_estrela1 / (massa_estrela1 + massa_estrela2),0,0)
pos_estrela2 = vp.vector(-separacao_estrelas * massa_estrela2 / (massa_estrela1 + massa_estrela2),0,0)


# vetores velocidade estrelas
vel_estrela1 = vp.vector(0,vel_inicialestrelas * massa_estrela1 / (massa_estrela1 + massa_estrela2),0)
vel_estrela2 = vp.vector(0,-vel_inicialestrelas * massa_estrela2 / (massa_estrela1 + massa_estrela2),0)
vel_estrela1 = vp.vector(0,0,0)
vel_estrela2 = vp.vector(0,vel_inicialestrelas,0)

#centro de massa entre as estrelas
centro_massa = (pos_estrela1*massa_estrela1 + pos_estrela2*massa_estrela2) / (massa_estrela1 + massa_estrela2)
vel_centro_massa = (vel_estrela1*massa_estrela1 + vel_estrela2*massa_estrela2) / (massa_estrela1 + massa_estrela2)

 
pos_estrela1 -= centro_massa
pos_estrela2 -= centro_massa
vel_estrela1 -= vel_centro_massa
vel_estrela2 -= vel_centro_massa


#tamanho das estrelas
tam_estrela1 = 1.0
tam_estrela2 = 0.5

inclinacao *= vp.pi/180


# vetores iniciais do planeta
pos_planeta = vp.vector(distancia_inicial_planeta*vp.cos(inclinacao), 0, distancia_inicial_planeta*vp.sin(inclinacao))
vel_planeta = vp.vector(0,vel_inicial_planeta,0)

#tamanho da esfera planeta
tam_planeta = 0.01

#Como estamos usando as unidades tranformadas para o sistema de unidades astronomicas e procurando uma orbita eliptica 4*pi^2*M sendo M=(m1+m2)
G= 4*vp.pi**2


#objetos esferas 
bola_planeta = vp.sphere(pos=pos_planeta, radius=tam_planeta, color=vp.vector(0.3,0.5,1), make_trail=True, retain = 2000)
bola_estrela1 = vp.sphere(pos=pos_estrela1, radius=tam_estrela1, color=vp.vector(1,1,0.3), make_trail=True, retain = 2000)
bola_estrela2= vp.sphere(pos=pos_estrela2, radius=tam_estrela2, color=vp.vector(1,0.6,0.3), make_trail=True, retain = 2000)
  
# centro de massa considerando a velocidade do planeta 
vel_centro_massa = (massa_planeta * vel_planeta  + vel_estrela1*massa_estrela1 + vel_estrela2*massa_estrela2) / (massa_planeta+ massa_estrela1 + massa_estrela2 )

vel_planeta -= vel_centro_massa


recheck = True
time=0
interval=0.0001


pos_planeta += vel_planeta * interval/2
ultima_distancia_estrela = distancia_estrela = (pos_planeta).mag  


steps = 0


while True:
#Simulação
       #Distancia entre cada corpo
    planeta_estrela1 = (pos_planeta - pos_estrela1).mag
    estrela1_estrela2 = (pos_estrela1 - pos_estrela2).mag
    planeta_estrela2 = (pos_planeta - pos_estrela2).mag

        #Direção do vetor unitario 
    direcao_planeta_estrela1 = (pos_estrela1 - pos_planeta).hat  #A.hat = A/|A|, a unit vector in the direction of the vector 
    direcao_estrela1_estrela2  = (pos_estrela1 - pos_estrela2).hat
    direcao_planeta_estrela2 = ( pos_estrela2 - pos_planeta).hat


        #Força entre cada corpo
    forca_estrela1_planeta = G * massa_estrela1 * massa_planeta / (planeta_estrela1**2) * direcao_planeta_estrela1   
    forca_estrela1_estrela2 = G * massa_estrela1 * massa_estrela2 / (estrela1_estrela2**2) *direcao_estrela1_estrela2  
    forca_estrela2_planeta = G * massa_estrela2 * massa_planeta / (planeta_estrela2 ** 2) * direcao_planeta_estrela2   
        
        #Reação
    forca_planeta_estrela1 = -forca_estrela1_planeta
    forca_estrela2_estrela1  = -forca_estrela1_estrela2
    forca_planeta_estrela2 = -forca_estrela2_planeta                                          
  
        #Acereleração de cada corpo
    acel_estrela1 = (forca_planeta_estrela1 + forca_estrela2_estrela1) / massa_estrela1
    acel_planeta = (forca_estrela1_planeta + forca_estrela2_planeta) / massa_planeta              
    acel_estrela2 = (forca_planeta_estrela2 + forca_estrela1_estrela2)  / massa_estrela2             
        
        #A velocidade em cada instante
    vel_estrela1 += acel_estrela1 * interval
    vel_estrela2 += acel_estrela2 * interval
    vel_planeta  += acel_planeta  * interval
    
        #Posição em cada instante 
    pos_estrela1 += vel_estrela1 * interval
    pos_estrela2 += vel_estrela2 * interval
    pos_planeta  += vel_planeta * interval
  
    time += interval
    
        #localização real do objeto bola
    bola_planeta.pos = pos_planeta
    bola_estrela1.pos = pos_estrela1
    bola_estrela2.pos = pos_estrela2

    vp.rate(1 / interval / segundos_por_anos)
    steps += 1

    ultima_ultima_distancia_estrela = ultima_distancia_estrela
    ultima_distancia_estrela = (pos_planeta).mag

#referencia do site projeto https://walterfreeman.github.io/ast101/binarysim.html