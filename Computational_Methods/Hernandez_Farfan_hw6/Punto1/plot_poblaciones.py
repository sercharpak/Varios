# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages

# <codecell>

nombreArchivo=(sys.argv[1])
#nombreArchivo='poblaciones_20.000000_30.000000.dat'
arregloStrings = nombreArchivo.split('_')
x0=float(arregloStrings[1])
y0=float(arregloStrings[2].strip('.dat'))
data=np.loadtxt(nombreArchivo)
t=data[:,0]
x=data[:,1]
xmax=max(x)
xmin=min(x)
y=data[:,2]
ymax=max(y)
ymin=min(y)

# <codecell>

fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
plt.xlabel('x - Presas')
plt.ylabel('y - Predadores')
plt.title('Poblacion Predadores vs Presas para condiciones iniciales x0='+str(x0)+' y y0='+str(y0))
plt.grid()
plt.plot(x,y)
pp = PdfPages(nombreArchivo.strip('.dat')+'.pdf')
pp.savefig(fig)
pp.close()
#plt.show()

# <codecell>

"""
#En este caso no me interesa la animacion
import matplotlib.animation as animation
line, = ax.plot([], [], lw=2)
def init():
    line.set_data([], [])
    return line,
# Funcion de Animacion
tam=len(t)
def animate(i):
    line.set_data(x, y)
    plt.grid()
    return line,
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=tam-1, interval=20, blit=True,repeat=False)
"""

# <codecell>


