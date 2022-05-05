# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:55:36 2022

@author: Rafael Kamimura
"""
from vpython import *
import numpy as np


G=6.67e-11
g = 9.81
m1=2e30
m2=0.7e30
rdist=2e10
M=m1+m2
Mp=6e24
r1=m2*rdist/(M)
r2=rdist*(1-m2/M)
star1=sphere(pos=vector(-r1,0,0), radius=1e9, color=color.yellow, make_trail=True)
star2=sphere(pos=vector(r2,0,0), radius=0.6e9, color=color.cyan, make_trail=True)


star1.m=m1
star2.m=m2

Rcm=(star1.m*star1.pos+star2.m*star2.pos)/M
Cm=sphere(pos=Rcm, radius=0.1e9, make_trail=True)
plan1=sphere(pos=Rcm, radius=0.3e9, color=color.blue, make_trail=True)
plan1.m=Mp


v1init=np.sqrt(G*star2.m**2/(rdist*M))
star1.p=star1.m*vector(0,v1init,0)
star2.p=-star1.p
plan1.v = np.sqrt(G*star1.m/mag(star1.pos - plan1.pos))
plan1.p = plan1.m*vector(0,plan1.v,0)



t=0
dt=1000

# Vetor aceleracao que aponta do planeta a estrela Kepler-16A
star1_planet_component = (plan1.pos - star1.pos) # Componentes do vetor posicao
mod_star1_planet_component = star1_planet_component.mag # Modulo do vetor posicao
star1_acc = G*m1*(star1_planet_component)/mod_star1_planet_component**3 # Vetor aceleracao
plan1.vecstar1 = vec(star1_acc*1E+10) # Criando vetor
attach_arrow(plan1, "vecstar1", color=color.yellow, shaftwidth=0.8*1E+8) # Anexando vetor


# Vetor acaeleracao que aponta do planeta a estrela Kepler-16B
star2_planet_component = plan1.pos - star2.pos # Componentes do vetor
mod_star2_planet_component = star2_planet_component.mag # Modulo do vetor posicao
star2_acc = G*m2*(star2_planet_component)/mod_star2_planet_component**3 # Vetor aceleracao
plan1.vecstar2 = vec(star2_acc*1E+10) # Criando vetor
attach_arrow(plan1, "vecstar2", color=color.cyan, shaftwidth=0.8*1E+8) # Anexando vetor


# Vetor que aponta do planeta a direcao do vetor da forca resultante
radial_component = star1_acc + star2_acc # Componentes do vetor
plan1.rad = vec(radial_component*1E+10) # Criando vetor
attach_arrow(plan1, "rad", color=color.blue, shaftwidth=0.8*1E+8) # Anexando vetor

while t<1e7:
    rate(100)
    r=star2.pos-star1.pos
    r1p = plan1.pos - star1.pos
    r2p = plan1.pos - star2.pos
    F12=-G*star1.m*star2.m*norm(r)/mag(r)**2
    F1p=-G*star1.m*plan1.m*norm(r1p)/mag(r1p)**2
    F2p=-G*star2.m*plan1.m*norm(r2p)/mag(r2p)**2
    star1.p = star1.p - (F12+F1p)*dt
    star2.p = star2.p + (F12-F2p)*dt
    plan1.p = plan1.p + (F1p+F2p)*dt 
    star1.pos = star1.pos + star1.p*dt/star1.m
    star2.pos = star2.pos + star2.p*dt/star2.m
    plan1.pos = plan1.pos + plan1.p*dt/plan1.m
    Rcm=(star1.m*star1.pos+star2.m*star2.pos)/M

    star1_planet_component = (plan1.pos - star1.pos)
    mod_star1_planet_component = star1_planet_component.mag
    star1_acc = -G*m1*(star1_planet_component)/mod_star1_planet_component**3 # Vetor aceleracao
    plan1.vecstar1 = vec(star1_acc*1E+9)

    star2_planet_component = (plan1.pos - star2.pos)
    mod_star2_planet_component = star2_planet_component.mag
    star2_acc = -G*m2*(star2_planet_component)/mod_star2_planet_component**3 # Vetor aceleracao
    plan1.vecstar2 = vec(star2_acc*1E+11)

    radial = (star1_acc + star2_acc)
    plan1.rad = vec(radial*1E+9)

    Cm.pos=Rcm
    t=t+dt
  
  


