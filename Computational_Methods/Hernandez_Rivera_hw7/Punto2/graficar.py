import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from matplotlib.backends.backend_pdf import PdfPages

datos = sys.argv[1]
nombreSinDat=datos.strip('.dat')
#estado_0.010000.dat
arregloStrings = nombreSinDat.split('_')
tf=float(arregloStrings[1].strip('.dat'))
data = np.loadtxt(datos)

figRho = plt.figure()
figRho.set_size_inches(18.5,10.5)
plt.xlabel('x (m)')
plt.ylabel('Densidad Rho (kg/m^3)')
plt.title('Densidad del sistema en el tiempo t='+str(tf)+' s')
plt.plot(data[:,0],data[:,3])

figVel = plt.figure()
figVel.set_size_inches(18.5,10.5)
plt.xlabel('x (m)')
plt.ylabel('Velocidad(m/s)')
plt.title('Velocidad del sistema en el tiempo t='+str(tf)+' s')
plt.plot(data[:,0],data[:,1])

figP = plt.figure()
figP.set_size_inches(18.5,10.5)
plt.xlabel('x (m)')
plt.ylabel('P (Pa)')
plt.title('Presion del sistema en el tiempo t='+str(tf)+' s')
plt.plot(data[:,0],data[:,2])

pp = PdfPages(nombreSinDat+'.pdf')
pp.savefig(figRho)
pp.savefig(figVel)
pp.savefig(figP)
pp.close()

