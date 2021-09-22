#!/usr/bin/env python
#-*- coding: 850 -*-

"""
Script que resuelve la ecuacion de flujo estacionario de Toth para condiciones complejas de valle intermontano
y un medio heterogeneo (Modelo de 2) a partir del metodo de elemento (Fenics)
Por: Ricardo Balam Chagoya Morales

Condiciones de la funcion que describe la geometria del nivel freático

Amplitud de la ecuación senoidal de pie de montaña = 1 m
Amplitud de la ecuación senoidal de intermontana = 10 m
longitud de onda de le ecuación senoidal = 0.01 m
Angulo de pendiente = 1.4711 rad 
"""


from __future__ import print_function
from fenics import *
import numpy as np
import random as ra

# Refinación de la malla  
numcelx=20
numcely=10
nele2=numcelx*numcely
# Conductividad hidraulica central del primer y el segundo estrato
k1= 10
k2= 10
# Datos del modelo
lonx=200
lony=100
pasox=lonx/numcelx
pasoy=lony/numcely

# Crea el mallado rectangular con intervalos de 20 metros
mesh = RectangleMesh(Point(0,0), Point(200,100), numcelx, numcely)
V = FunctionSpace(mesh, 'P', 1)

# Se definen las condiciones de frontera superior
# Condicion de pie de montana
u_I = Expression('500-x[0]*tan(1.4711)+1*(sin((1000*x[0])/cos(1.4711))/cos(1.4711))', degree=1)

def frontera_I(x, dentro_frontera):
 tol = 1E-14
 if dentro_frontera:
   if abs(x[1]-100)<=tol and abs(x[0])<=40:
    return True
   else:
    return False
 else:
  return False

# Condicion de valle intermontano 
u_C = Expression('100+10*(sin(1000*x[0]))', degree=1)

def frontera_C(x, dentro_frontera):
 tol = 1E-14
 if dentro_frontera:
   if abs(x[1]-100)<=tol and 40<=abs(x[0])<=160:
    return True
   else:
    return False
 else:
  return False

# COndicion de pie de montana derecho
u_D = Expression('-1500+x[0]*tan(1.4711)+1*(sin((1000*x[0])/cos(1.4711))/cos(1.4711))', degree=1)

def frontera_D(x, dentro_frontera):
 tol = 1E-14
 if dentro_frontera:
   if abs(x[1]-100)<=tol and abs(x[0])>=160:
    return True
   else:
    return False
 else:
  return False


CondIzquierda = DirichletBC(V, u_I, frontera_I)
CondCentral = DirichletBC(V, u_C, frontera_C)
CondDerecha = DirichletBC(V, u_D, frontera_D)

bc = [CondIzquierda,CondCentral,CondDerecha]

# Se define la heterogeneidad del sistema, aplicando para cada nodo un valor central

#class K(Expression):
# def set_k_values(self,k,num_ele):
#   self.k=k
#   self.k[100]=0
# def eval(self, value, x):
#   tol = 1E-14
#   cont=0
#   for i in range(0+pasox,lonx+1,pasox):
#    for j in range(0+pasoy,lony+1,pasoy):
#     if (x[0]-i)<tol and (x[1]-j)<tol:
#      value[0] = self.k[cont]
#      cont=cont+1

class K(Expression):
 def set_k_values(self,k,num_ele):
   self.k=k
   self.k[100]=0
 def eval(self, value, x):
   tol = 1E-14
   cont=0
   if abs(x[1]-10)<=tol:
    value[0] = 1
   else:
    value[0] = 1000
       
      

# and (x[0]-(i-10))>tol and (x[1]-(j-10)>tol): 
# Inicializando subdominio del medio
# Se lee el archivo de texto que contiene los resultados de la simulación no condicional 
num_ele=200
k = np.loadtxt("simnocond1.txt",delimiter=',',skiprows=1,usecols=[6])

#k=(1,10)
kappa = K(degree=0)
kappa.set_k_values(k,num_ele)

# Se define el problema variacional
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)
a = kappa*dot(grad(u), grad(v))*dx
g = Constant(0)
L = f*v*dx - g*v*ds

# Se realiza el calculo de la solucion

u = Function(V)
solve(a == L, u, bc)

# Se plotea la solucion
plot(u)
plot(mesh)

# Se salva la solucion en un archivo vtkfile 
vtkfile = File('poisson/solution.pvd')
vtkfile << u

# Se calcula el error L2
error_L2 = errornorm(u_I, u, 'L2')

# Se calcula el error maximo en los vertices
vertex_values_u_I = u_I.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
error_max = np.max(np.abs(vertex_values_u_I - vertex_values_u))

# Impresion de los errores
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Se realiza el ploteo de la solucion
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax= plt.subplot(111)
plot(u)
im=plot(-grad(u))

plt.title('Flujo de agua subterranea en valle intermontano')
plt.ylabel('Elevacion [m]')
plt.xlabel('Distancia [m]')
plt.colorbar(im)
plt.show()
