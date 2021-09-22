"""
Script que resuelve la ecuacion de flujo estacionario de Toth
a partir de elemento finito
Por: Ricardo Balam Chagoya Morales
"""


from __future__ import print_function
from fenics import *
import random as ra 

# Crea el mallado rectangular con intervalos de 20 metros
mesh = RectangleMesh(Point(0,0), Point(200,100), 10, 5)
V = FunctionSpace(mesh, 'P', 1)

# Se definen las condiciones de frontera superior
u_D = Expression('0.02*x[0]+100', degree=1)

tol = 1E-14
def boundary_D(x, on_boundary):
 if on_boundary:
   if near(x[1], 100, tol):
      return True
   else:
      return False
 else:
   return False

bc = DirichletBC(V, u_D, boundary_D)

#Se define la heterogeneidad del medio

class K(Expression):
 def set_k_values(self, k_0, k_1):
  self.k_0, self.k_1 = k_0, k_1
 def eval(self, value, x):
  tol = 1E-14
  if (x[1] <= 50 + tol): 
   for i in range(0,201,20):
    for j in range(0,50,20):
     if (x[0]-i)<tol and (x[1]-j)<tol: 
      value[0] = self.k_0+ra.random()
  else:
   for m in range(0,201,20):
    for n in range(60,101,20):
     if (x[0]-m)<tol and (x[1]-n)<tol:
      value[0] = self.k_1+ra.random()

# Inicializando subdominio del medio
kappa = K(degree=0)
kappa.set_k_values(1,10)

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
error_L2 = errornorm(u_D, u, 'L2')

# Se calcula el error maximo en los vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Impresion de los errores
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Se realiza el ploteo de la solucion
import matplotlib.pyplot as plt
plot(u)
plt.show()

