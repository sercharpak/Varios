# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#Codigo realizado para completar el razonamiento del punto 2 de la tarea 2 de Biologia Sintetica en el semestre 2015-I en Uniandes
#Codigo realizado por Sergio Hernandez

# <codecell>

import numpy as np
import networkx as nx

# <markdowncell>

# Importamos la matriz de adjacencia que se fabrico en Excel con base a la imagen adjuntada en la tarea.

# <codecell>

matriz_A=np.loadtxt('Tabla_2.tsv')

# <markdowncell>

# Ahora obtenemos una matriz de NetworkX con todas las longitudes de todos los menores caminos gracias al comando all_pairs_dijkstra_path_length luego de haber pasado la matriz al formato de NetworkX para su correcto uso.

# <codecell>

G = nx.from_numpy_matrix(matriz_A, create_using=nx.DiGraph())
matriz_C_Todos=nx.all_pairs_dijkstra_path_length(G)

# <markdowncell>

# Los cuatro nodos que nos interesan son PABN2 (Aca, 5), PABN3(Aca, 6), FIPS3(Aca, 9) y PCSF1(Aca, 16).

# <markdowncell>

# Ahora borremos a PABN2 (Nodo 5) y miremos cuantos caminos cambiaron.

# <codecell>

matriz_A_sin_5=np.loadtxt('Tabla_2.tsv')
matriz_A_sin_5[:,5]=0.0
matriz_A_sin_5[5,:]=0.0
G_sin_5 = nx.from_numpy_matrix(matriz_A_sin_5, create_using=nx.DiGraph())
matriz_C_sin_5=nx.all_pairs_dijkstra_path_length(G_sin_5)

# <codecell>

matriz_A_sin_6=np.loadtxt('Tabla_2.tsv')
matriz_A_sin_6[:,6]=0.0
matriz_A_sin_6[6,:]=0.0
G_sin_6 = nx.from_numpy_matrix(matriz_A_sin_6, create_using=nx.DiGraph())
matriz_C_sin_6=nx.all_pairs_dijkstra_path_length(G_sin_6)

# <codecell>

matriz_A_sin_9=np.loadtxt('Tabla_2.tsv')
matriz_A_sin_9[:,9]=0.0
matriz_A_sin_9[9,:]=0.0
G_sin_9 = nx.from_numpy_matrix(matriz_A_sin_9, create_using=nx.DiGraph())
matriz_C_sin_9=nx.all_pairs_dijkstra_path_length(G_sin_9)

# <codecell>

matriz_A_sin_16=np.loadtxt('Tabla_2.tsv')
matriz_A_sin_16[:,16]=0.0
matriz_A_sin_16[16,:]=0.0
G_sin_16 = nx.from_numpy_matrix(matriz_A_sin_16, create_using=nx.DiGraph())
matriz_C_sin_16=nx.all_pairs_dijkstra_path_length(G_sin_16)

# <codecell>

contador_diferencias_5=0
contador_diferencias_6=0
contador_diferencias_9=0
contador_diferencias_16=0
for i in range (0,24):
    for j in range (0,24):
        if (i!=5 and j!=5):
            if(matriz_C_Todos[i][j] != matriz_C_sin_5[i][j]):
                contador_diferencias_5=contador_diferencias_5+1
        if(i!=6 and j!=6):
            if(matriz_C_Todos[i][j]!= matriz_C_sin_6[i][j]):
                contador_diferencias_6=contador_diferencias_6+1
        if(i!=9 and j!=9):
            if(matriz_C_Todos[i][j]!= matriz_C_sin_9[i][j]):
                contador_diferencias_9=contador_diferencias_9+1
        if(i!=16 and j!=16):
            if(matriz_C_Todos[i][j]!= matriz_C_sin_16[i][j]):
                contador_diferencias_16=contador_diferencias_16+1
print "Numero de caminos minimos cuya longitud fue modificada al eliminar PABN2 fueron: "+str(contador_diferencias_5)
print "Numero de caminos minimos cuya longitud fue modificada al eliminar PABN3 fueron: "+str(contador_diferencias_6)
print "Numero de caminos minimos cuya longitud fue modificada al eliminar FIPS3 fueron: "+str(contador_diferencias_9)
print "Numero de caminos minimos cuya longitud fue modificada al eliminar PCSF1 fueron: "+str(contador_diferencias_16)

# <markdowncell>

# Podemos observar que eliminar PABN2 (0) y PABN3 (0) no modifican ninguna longitud de camino minima mientras que eliminar FIPS3 (12) y PCSF1 (8) si. 
# Por lo tanto se puede reorganizar estos nodos en orden de importancia: FIPS3, PCSF1 y PABN2 y PABN3.

