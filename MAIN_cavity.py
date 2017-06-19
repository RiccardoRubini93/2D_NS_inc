from __future__ import division
import numpy as np
from matplotlib import pyplot, cm
import math
from Solver import *

#setting up domain
nx = 201
ny = 101
nt = 10000
nit = 500
L_x = 2
L_y = 1
tau = 4
dx = L_x / (nx - 1)
dy = L_y / (ny - 1)
x = np.linspace(0, L_x, nx)
y = np.linspace(0, L_y, ny)
X, Y = np.meshgrid(x, y)

#some physical parameters

beta = 1
V_top = 20
T_top = 100
rho = 1
nu = 0.2
alpha = 0.2
dt = 0.0001

#initializa matrices

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
T = np.zeros((ny, nx))
p = np.zeros((ny, nx)) 
b = np.zeros((ny, nx))

#computing solution

u, v, p, T = solver(nt, u, v,T, dt, dx, dy, p, rho, nu, alpha,V_top,T_top,nx,ny,nit)

magU = (u**2+v**2)**0.5

#plotting

f1= pyplot.figure(1)
pyplot.rcParams['axes.facecolor'] = 'white'
f1.patch.set_facecolor('white')

pyplot.subplot(121)
pyplot.contourf(X, Y, T, 8, alpha=.75, cmap='jet')
C = pyplot.contour(X, Y, T, 8, colors='black', linewidth=.5) 
pyplot.clabel(C, inline=1, fontsize=10)
#pyplot.colorbar()
pyplot.xlabel('X')
pyplot.ylabel('Y')
pyplot.xlim([0,L_x])
pyplot.ylim([0,L_y])
pyplot.title('Temperature')

pyplot.subplot(122)
#f2 = pyplot.figure(2)
#pyplot.rcParams['axes.facecolor'] = 'white'
#f2.patch.set_facecolor('white')
pyplot.streamplot(x, y, u, v, density=5, color=magU,linewidth=1)
pyplot.xlabel('X')
pyplot.ylabel('Y')
pyplot.xlim([0,L_x])
pyplot.ylim([0,L_y])
pyplot.title('Streamlines')

pyplot.show()
