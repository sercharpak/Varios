# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#Codigo realizado para completar el punto 3 de la tarea 2 de Biologia Sintetica en el semestre 2015-I en Uniandes
#Codigo realizado por Sergio Hernandez

# <codecell>

import numpy as np
import matplotlib.pyplot as plt

# <markdowncell>

# Numeral a: Para la ecuación de diferencia $x_{i+1} = -3 x_i (x_i - 1 ) $ dibuje las trayectorias para dos puntos iniciales cualquiera en el espacio de fase $x_{i+1}$ vs. $x_i$ y en el plano $x_i$  vs.$i$. 

# <markdowncell>

# Podemos observar que tenemos la libertad de escoger los puntos iniciales. Vamos a escoger dos puntos aleatorios enteros entre 0 y 20 y vamos a graficar la secuencia para 100 puntos.

# <codecell>

puntos_iniciales_int=np.random.randint(20, size=2)

# <codecell>

print "Los puntos iniciales son: "+str(puntos_iniciales_int)

# <codecell>

#Inicializamos
x_i_1_punto_1=np.zeros(100)
x_i_punto_1=np.zeros(100)
x_i_1_punto_2=np.zeros(100)
x_i_punto_2=np.zeros(100)
arreglo_i=np.zeros(100)
x_i_punto_1[0]=puntos_iniciales_int[0]
x_i_punto_2[0]=puntos_iniciales_int[1]
#Arranquemos
for i in range (1,100):
    x_i_1_punto_1[i]= - 3 * x_i_punto_1[i-1] * (x_i_punto_1[i-1] -1)
    x_i_punto_1[i]=x_i_1_punto_1[i]
    x_i_1_punto_2[i]= - 3 * x_i_punto_2[i-1] * (x_i_punto_2[i-1] -1)
    x_i_punto_2[i]=x_i_1_punto_2[i]
    arreglo_i[i]=i

# <codecell>

#Ahora grafiquemos
#Primer Punto
plt.figure(1)
plt.plot(x_i_punto_1,x_i_1_punto_1)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[0]))
plt.savefig('Punto_3_a_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales_int[0])+'.png')
plt.figure(2)
plt.plot(arreglo_i,x_i_1_punto_1)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[0]))
plt.savefig('Punto_3_a_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales_int[0])+'.png')
#Segundo Punto
plt.figure(3)
plt.plot(x_i_punto_2,x_i_1_punto_2)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[1]))
plt.savefig('Punto_3_a_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales_int[1])+'.png')
plt.figure(4)
plt.plot(arreglo_i,x_i_1_punto_2)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[1]))
plt.savefig('Punto_3_a_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales_int[1])+'.png')
plt.show()

# <markdowncell>

# Se puede observar que no se puede sacar muchas conclusiones al respecto. Salvo que si el punto inicial es 1 no mucho ocurre. Se sospecha que se debe a la divergencia de la ecuacion de diferencia. Vamos a escoger dos numeros aleatorios entre 0 y 1 para que converja.

# <codecell>

puntos_iniciales=np.random.random(2)
print "Ahora los puntos iniciales son: "+str(puntos_iniciales)

# <codecell>

#Inicializamos
x_i_1_punto_1=np.zeros(100)
x_i_punto_1=np.zeros(100)
x_i_1_punto_2=np.zeros(100)
x_i_punto_2=np.zeros(100)
arreglo_i=np.zeros(100)
x_i_punto_1[0]=puntos_iniciales[0]
x_i_punto_2[0]=puntos_iniciales[1]
#Arranquemos
for i in range (1,100):
    x_i_1_punto_1[i]= - 3 * x_i_punto_1[i-1] * (x_i_punto_1[i-1] -1)
    x_i_punto_1[i]=x_i_1_punto_1[i]
    x_i_1_punto_2[i]= - 3 * x_i_punto_2[i-1] * (x_i_punto_2[i-1] -1)
    x_i_punto_2[i]=x_i_1_punto_2[i]
    arreglo_i[i]=i

# <codecell>

#Ahora grafiquemos
#Primer Punto
plt.figure(5)
plt.plot(x_i_punto_1,x_i_1_punto_1)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[0]))
plt.savefig('Punto_3_a_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales[0])+'.png')
plt.figure(6)
plt.plot(arreglo_i,x_i_1_punto_1)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[0]))
plt.savefig('Punto_3_a_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales[0])+'.png')
#Segundo Punto
plt.figure(7)
plt.plot(x_i_punto_2,x_i_1_punto_2)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[1]))
plt.savefig('Punto_3_a_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales[1])+'.png')
plt.figure(8)
plt.plot(arreglo_i,x_i_1_punto_2)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -3 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[1]))
plt.savefig('Punto_3_a_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales[1])+'.png')
plt.show()

# <markdowncell>

# Podemos observar que el comportamiento si el punto inicial esta entre 0 y 1 es oscilatorio luego de cierto numero de iteraciones.

# <markdowncell>

# Numeral b: Para la ecuación de diferencia $x_{i+1}=4 x_i (x_i - 1)$, dibuje las trayectorias para dos puntos iniciales cualquiera en el espacio de fase $x_{i+1}$ vs. $x_i $y en el plano $x_i$  vs.$i$.

# <markdowncell>

# Vamos a repetir exactamente el mismo procedimiento que en el numeral anterior.

# <codecell>

puntos_iniciales_int=np.random.randint(20, size=2)
print "Los puntos iniciales son: "+str(puntos_iniciales_int)

# <codecell>

#Inicializamos
x_i_1_punto_1=np.zeros(100)
x_i_punto_1=np.zeros(100)
x_i_1_punto_2=np.zeros(100)
x_i_punto_2=np.zeros(100)
arreglo_i=np.zeros(100)
x_i_punto_1[0]=puntos_iniciales_int[0]
x_i_punto_2[0]=puntos_iniciales_int[1]
#Arranquemos
for i in range (1,100):
    x_i_1_punto_1[i]= -4 * x_i_punto_1[i-1] * (x_i_punto_1[i-1] -1)
    x_i_punto_1[i]=x_i_1_punto_1[i]
    x_i_1_punto_2[i]= -4 * x_i_punto_2[i-1] * (x_i_punto_2[i-1] -1)
    x_i_punto_2[i]=x_i_1_punto_2[i]
    arreglo_i[i]=i

# <codecell>

#Ahora grafiquemos
#Primer Punto
plt.figure(9)
plt.plot(x_i_punto_1,x_i_1_punto_1)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[0]))
plt.savefig('Punto_3_b_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales_int[0])+'.png')
plt.figure(10)
plt.plot(arreglo_i,x_i_1_punto_1)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[0]))
plt.savefig('Punto_3_b_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales_int[0])+'.png')
#Segundo Punto
plt.figure(11)
plt.plot(x_i_punto_2,x_i_1_punto_2)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[1]))
plt.savefig('Punto_3_b_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales_int[1])+'.png')
plt.figure(12)
plt.plot(arreglo_i,x_i_1_punto_2)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales_int[1]))
plt.savefig('Punto_3_b_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales_int[1])+'.png')
plt.show()

# <codecell>

puntos_iniciales=np.random.random(2)
print "Ahora los puntos iniciales son: "+str(puntos_iniciales)

# <codecell>

#Inicializamos
x_i_1_punto_1=np.zeros(100)
x_i_punto_1=np.zeros(100)
x_i_1_punto_2=np.zeros(100)
x_i_punto_2=np.zeros(100)
arreglo_i=np.zeros(100)
x_i_punto_1[0]=puntos_iniciales[0]
x_i_punto_2[0]=puntos_iniciales[1]
#Arranquemos
for i in range (1,100):
    x_i_1_punto_1[i]= -4.0 * x_i_punto_1[i-1] * (x_i_punto_1[i-1] -1)
    x_i_punto_1[i]=x_i_1_punto_1[i]
    x_i_1_punto_2[i]= -4.0 * x_i_punto_2[i-1] * (x_i_punto_2[i-1] -1)
    x_i_punto_2[i]=x_i_1_punto_2[i]
    arreglo_i[i]=i

# <codecell>

#Ahora grafiquemos
#Primer Punto
plt.figure(13)
plt.plot(x_i_punto_1,x_i_1_punto_1)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[0]))
plt.savefig('Punto_3_b_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales[0])+'.png')
plt.figure(14)
plt.plot(arreglo_i,x_i_1_punto_1)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[0]))
plt.savefig('Punto_3_b_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales[0])+'.png')
#Segundo Punto
plt.figure(15)
plt.plot(x_i_punto_2,x_i_1_punto_2)
plt.xlabel("$x_{i}$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $x_{i}$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[1]))
plt.savefig('Punto_3_b_x_i_1_vs_x_i_Iniciando_'+str(puntos_iniciales[1])+'.png')
plt.figure(16)
plt.plot(arreglo_i,x_i_1_punto_2)
plt.xlabel("$i$")
plt.ylabel("$x_{i+1}$")
plt.title("$x_{i+1}$ vs. $i$ para $x_{i+1} = -4 x_i (x_i - 1 ) $ con punto inicial: "+str(puntos_iniciales[1]))
plt.savefig('Punto_3_b_x_i_1_vs_i_Iniciando_'+str(puntos_iniciales[1])+'.png')
plt.show()

# <markdowncell>

# Podemos observar que obtenemos un comportamiento oscilatorio, al igual que el numeral anterior. Sin embargo las oscilaciones son mucho mas violentas en este caso. Modificar este factor multiplicatorio (en este caso por una unidad) tiene un gran impacto en el sistema.

