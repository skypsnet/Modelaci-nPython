#!/usr/bin/env python
#-*- coding: 850 -*-

from __future__ import print_function
from fenics import *
import numpy as np
import random as ra
import matplotlib.pyplot as plt
import sys
reload(sys)
sys.setdefaultencoding('utf-8')

# Se crea la malla donde se define el dominio

mesh= RectangleMesh(Point(0,0),Point(200,100),20,10)
V = FunctionSpace(mesh, 'P', 1)
# Se definen las condiciones de fronteras

def frontera_I(x,frontera):
 tol=1E-14    
 if frontera:
  if x[0]<=tol: 
   return True
  else:
   return False 
 else:
  return False  

def frontera_D(x,frontera):
 tol=1E-14    
 if frontera:
  if abs(x[0]-200)<=tol: 
   return True
  else:
   return False 
 else:
  return False  


F_D = DirichletBC(V, Constant(10), frontera_D)
F_I = DirichletBC(V, Constant(100), frontera_I)
bc = [F_D,F_I]

#Se define el problema variacional

u=TrialFunction(V)
v=TestFunction(V)
f=Constant(0)
a=dot(grad(u),grad(v))*dx
g=Constant(0)
L=f*v*dx-g*v*ds

# Se realiza el cálculo de la solución

u= Function(V)
solve(a==L,u,bc)
tau=project(-grad(u))
flujo=tau.vector()[:] 


# Se organizan los valores para la exportación de los datos
tau=project(grad(u))

xx=np.loadtxt("Coordenadascentroide.csv",delimiter=',',skiprows=1,usecols=[1])
yy=np.loadtxt("Coordenadascentroide.csv",delimiter=',',skiprows=1,usecols=[2])
datos=np.zeros((400,7))
datos[:,0]=xx
datos[:,1]=yy
datos[:,4]=1
for i in range(0,399):
 datos[i,2]=-100*tau(xx[i],yy[i])[0]
 datos[i,3]=-100*tau(xx[i],yy[i])[1]

np.savetxt("SoluciónFlujo1.1",datos)

# Se gráfica la descarga específica

figure()
plt.plot(datos[360:,0],datos[360:,3])


# Ploteo de la solución
  
#im=plot(u)
#plt.colorbar(im) 
#plot(u)
#plot(mesh)
#lo=plot(-grad(u))
figure()
carga=plot(u)
flujo=plt.quiver(xx,yy,datos[:,2],datos[:,3])
plt.colorbar(carga)
plt.title('Flujo de agua subterránea')
plt.ylabel('Elevación [m]')
plt.xlabel('Distancia [m]')
plt.show()



 
