import numpy as np
import sys,string,os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

datos = sys.argv[1]
data = np.loadtxt(datos)
nombreSinDat=datos.strip('.dat')
arregloStrings = nombreSinDat.split('_')
rho=float(arregloStrings[1].strip('.dat'))

x=np.linspace(0,100,101)
t=np.linspace(0,120,121)
xy,ty=np.meshgrid(x,t)
fig3D = plt.figure()
fig3D.set_size_inches(20.0,11.0)
ax = Axes3D(plt.gcf())
ax.plot_wireframe(xy,ty,data)
ax.set_title("Evolucion temporal cuerda para rho="+str(rho)+"kg/m^3")
ax.set_xlabel('x (m)')
ax.set_ylabel('t (s)')
ax.set_zlabel('U (m)')
pp = PdfPages('cuerda_'+str(rho)+'.pdf')
pp.savefig(fig3D)
pp.close()

plt.show();
