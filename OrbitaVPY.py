# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:55:36 2022

@author: Rafael Kamimura
"""
from vpython import *
import numpy as np
import matplotlib.pyplot as plt


G=6.67e-11
m1=2e30
m2=0.7e30
rdist=2e10
M=m1+m2

r1=m2*rdist/(M)
r2=rdist*(1-m2/M)
star1=sphere(pos=vector(-r1,0,0), radius=1e9, color=color.yellow, make_trail=True)
star2=sphere(pos=vector(r2,0,0), radius=0.6e9, color=color.cyan, make_trail=True)

star1.m=m1
star2.m=m2

Rcm=(star1.m*star1.pos+star2.m*star2.pos)/M
Cm=sphere(pos=Rcm, radius=0.3e9, make_trail=True)

v1init=np.sqrt(G*star2.m**2/(rdist*M))
star1.p=star1.m*vector(0,v1init,0)
star2.p=-star1.p

t=0
dt=1000

while t<1e7:
  rate(100)
  r=star2.pos-star1.pos
  F12=-G*star1.m*star2.m*norm(r)/mag(r)**2
  star1.p=star1.p-F12*dt
  star2.p=star2.p+F12*dt
  star1.pos=star1.pos+star1.p*dt/star1.m
  star2.pos=star2.pos+star2.p*dt/star2.m
  Rcm=(star1.m*star1.pos+star2.m*star2.pos)/M
  Cm.pos=Rcm
  
  t=t+dt
