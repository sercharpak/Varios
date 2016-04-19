# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as Gamma

# <codecell>

#Hecho por Sergio Daniel Hernandez Charpak
#Codigo para el punto 3 de la tarea 3 de Biologia Sintetica en el semestre 2015-1 en Uniandes dictado por Juan Manuel Pedraza
#Graficas distribuciones de probabilidad

# <codecell>

def dist_poisson(lam,x):
    return (lam**x)*(np.exp(-lam))/(Gamma(x))
def dist_exponential(lam,x):
    return (lam)*(np.exp(-lam*x))
def dist_gaussian(mu,sigma,x):
    return (1.0/(2.0*sigma*np.pi)**(0.5))*np.exp(-0.5*((x-mu)**2.0)/(sigma**2.0))

# <codecell>

#Formemos los linspaces
x_p=np.linspace(0,50,1000)
x_exp=np.linspace(0,10,1000)
x_g=np.linspace(-10,10,2000)
lam_p=3.0
y_p=dist_poisson(lam_p,x_p)
lam_exp=1.0
y_exp=dist_exponential(lam_exp,x_exp)
mu=0.0 #normal
sigma=1.0 #normal
y_g=dist_gaussian(mu,sigma,x_g)

# <codecell>

plt.figure(1)
plt.title('Distribucion de Poisson con lambda = '+str(lam_p))
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x_p, y_p)
plt.show()

# <codecell>

plt.figure(2)
plt.title('Distribucion Exponencial con lambda = '+str(lam_exp))
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x_exp, y_exp)
plt.show()

# <codecell>

plt.figure(3)
plt.title('Distribucion Gaussiana con mu='+str(mu) +' y sigma='+str(sigma))
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x_g, y_g)
plt.show()

# <codecell>


