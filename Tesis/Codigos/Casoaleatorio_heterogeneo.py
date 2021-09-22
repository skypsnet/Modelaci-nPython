#!/usr/bin/env python
#-*- coding: 850 -*-

from __future__ import print_function
from fenics import *
import numpy as np
import random as ra
import matplotlib.pyplot as plt

# Se crea la malla donde se define el dominio

mesh= RectangleMesh(Point(0,0),Point(200,100),20,10)
V = FunctionSpace(mesh, 'P', 1)
b=200/20
h=100/10
plt.figure()
plot(mesh)

# Se definen las condiciones de fronteras

def frontera_D(x,frontera):
 tol=1E-14 
 if frontera:
  if abs(x[0]-200)<=tol: 
   return True
  else:
   return False 
 else:
  return False  

def frontera_I(x,frontera):
 tol=1E-14 
 if frontera:
  if abs(x[0])<=tol: 
   return True
  else:
   return False 
 else:
  return False

F_I = DirichletBC(V, Constant(100), frontera_I)
F_D = DirichletBC(V, Constant(10), frontera_D)
bc = [F_I,F_D]

# Se define la heterogeneidad en el sistema

# Se realiza el marcado de la malla
k=range(0,400,1)
subdomains=CellFunction('size_t',mesh,0)
cont=0
for cell in cells(mesh):
 subdomains[cell]=k[cont]
 cont=cont+1

plt.figure()
im1=plot(subdomains)
print(subdomains)
plt.colorbar(im1)

# Se crea un espacio de funciones que simbolizan la conductividad hidráulica   
V0= FunctionSpace(mesh,"DG",0)
k=Function(V0)

# Se asigna para cada valor del marcado un valor del conjunto de conductividades
ka = np.loadtxt("Pruebayy.txt",delimiter=',',skiprows=1,usecols=[6])
k_values=np.exp(ka)
cont=0
for cell_no in cells(mesh):
 subdomain_no=subdomains[cell_no]
 k.vector()[cont]=k_values[subdomain_no]
 cont=cont+1

#Se define el problema variacional

u=TrialFunction(V)
v=TestFunction(V)
f=Constant(0)
a=k*dot(grad(u),grad(v))*dx
g=Constant(0)
L=f*v*dx-g*v*ds

# Se realiza el cálculo de la solución

u= Function(V)
solve(a==L,u,bc)

# Se organizan los valores para la exportación de los datos
tau=project(grad(u))

xx=np.loadtxt("Coordenadascentroide.csv",delimiter=',',skiprows=1,usecols=[1])
yy=np.loadtxt("Coordenadascentroide.csv",delimiter=',',skiprows=1,usecols=[2])
datos=np.zeros((400,5))
datos[:,0]=xx
datos[:,1]=yy
datos[:,4]=k_values
for i in range(0,399):
 datos[i,2]=tau(xx[i],yy[i])[0]
 datos[i,3]=tau(xx[i],yy[i])[1]


# Ploteo de la solución
plt.figure()
ax= plt.subplot(111)
im=plot(u)
#plot(mesh)
plt.colorbar(im)
#plot(u)
#plot(mesh)
plot(-grad(u))
plt.title('Flujo de agua subterranea')
plt.ylabel('Elevacion [m]')
plt.xlabel('Distancia [m]')
plt.show()


