# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline
import numpy as np
import matplotlib.pyplot as plt
import time
#Hecho por Sergio Daniel Hernandez Charpak
#Codigo para el punto 4 de la tarea 3 de Biologia Sintetica en el semestre 2015-1 en Uniandes dictado por Juan Manuel Pedraza
#Simulacion Estocastica Primitiva

# <codecell>

def evento_rna(num_rna,delta_t):
    e1=np.random.random()
    e2=np.random.random()
    if((e1<(k_r*delta_t)) and (e2>=(gamma_r*num_rna*delta_t))):
        num_rna=num_rna+1.0
    elif((e1>=(k_r*delta_t)) and (e2<(gamma_r*num_rna*delta_t))):
        if(num_rna>0.0):
            num_rna=num_rna-1.0
    return num_rna
def evento_prot(num_rna,num_prot,delta_t):
    e1=np.random.random()
    e2=np.random.random()
    if((e1<(k_p*num_rna*delta_t)) and (e2>=(gamma_p*num_prot*delta_t))):
        num_prot=num_prot+1.0
    elif((e1>=(k_p*num_rna*delta_t)) and (e2<(gamma_p*num_prot*delta_t))):
        if(num_prot>0.0):
            num_prot=num_prot-1.0
    return num_prot
def evento_rna_retro(num_rna,num_prot,delta_t):
    e1=np.random.random()
    e2=np.random.random()
    if((e1<((k_r/(1.0+((num_prot/K)**2.0)))*delta_t)) and (e2>=(gamma_r*num_rna*delta_t))):
        num_rna=num_rna+1.0
    elif((e1>=((k_r/(1.0+((num_prot/K)**2.0)))*delta_t)) and (e2<(gamma_r*num_rna*delta_t))):
        if(num_rna>0.0):
            num_rna=num_rna-1.0
    return num_rna
def evento_prot_retro(num_rna,num_prot,delta_t):
    e1=np.random.random()
    e2=np.random.random()
    if((e1<(k_p*num_rna*delta_t)) and (e2>=(gamma_p*num_prot*delta_t))):
        num_prot=num_prot+1.0
    elif((e1>=(k_p*num_rna*delta_t)) and (e2<(gamma_p*num_prot*delta_t))):
        if(num_prot>0.0):
            num_prot=num_prot-1.0
    return num_prot

# <codecell>

def simulacion_primitiva(num_cells,delta_t,t_limite):
    num_puntos=int(t_limite/delta_t)
    matriz_rna=np.zeros((num_cells,num_puntos))
    matriz_prot=np.zeros((num_cells,num_puntos))
    
    for j in range (num_cells):
        
        t_walk=np.zeros(num_puntos)
        r_walk=np.zeros(num_puntos)
        prot_walk=np.zeros(num_puntos)
        for i in range (1,num_puntos):
            t_walk[i]=t_walk[i-1]+delta_t
            r_walk[i]=evento_rna(r_walk[i-1],delta_t)
            prot_walk[i]=evento_prot(r_walk[i-1],prot_walk[i-1],delta_t)
        matriz_rna[j]=r_walk
        matriz_prot[j]=prot_walk
        
    r_mean=np.zeros(num_puntos)
    r_std=np.zeros(num_puntos)
    prot_mean=np.zeros(num_puntos)
    prot_std=np.zeros(num_puntos)
    
    for i in range(num_puntos):
        r_mean[i]=np.mean(matriz_rna[:,i])
        r_std[i]=np.std(matriz_rna[:,i])
        prot_mean[i]=np.mean(matriz_prot[:,i])
        prot_std[i]=np.std(matriz_prot[:,i])
        
    return t_walk, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std
def simulacion_primitiva_retro(num_cells,delta_t,t_limite):
    num_puntos=int(t_limite/delta_t)
    matriz_rna=np.zeros((num_cells,num_puntos))
    matriz_prot=np.zeros((num_cells,num_puntos))
    
    for j in range (num_cells):
        
        t_walk=np.zeros(num_puntos)
        r_walk=np.zeros(num_puntos)
        prot_walk=np.zeros(num_puntos)
        for i in range (1,num_puntos):
            t_walk[i]=t_walk[i-1]+delta_t
            r_walk[i]=evento_rna_retro(r_walk[i-1],prot_walk[i-1],delta_t)
            prot_walk[i]=evento_prot_retro(r_walk[i-1],prot_walk[i-1],delta_t)
        matriz_rna[j]=r_walk
        matriz_prot[j]=prot_walk
        
    r_mean=np.zeros(num_puntos)
    r_std=np.zeros(num_puntos)
    prot_mean=np.zeros(num_puntos)
    prot_std=np.zeros(num_puntos)
    
    for i in range(num_puntos):
        r_mean[i]=np.mean(matriz_rna[:,i])
        r_std[i]=np.std(matriz_rna[:,i])
        prot_mean[i]=np.mean(matriz_prot[:,i])
        prot_std[i]=np.std(matriz_prot[:,i])
        
    return t_walk, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std

# <codecell>

#Constantes biologicamente relevantes
k_r=1.0 #rna/min
k_p=50.0 # 1 prot/(min*rna)
K=0.5 #cte de retroalimentacion
gamma_p=1.0/30.0 #  prot/min
gamma_r=1.0/5.0 #  rna/min

# <codecell>

delta_t=0.005
t_limite=300
num_cells=150

# <codecell>

#Correr el no Retroalimentado
t_ini=time.time()
t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std = simulacion_primitiva(num_cells,delta_t,t_limite)
t_fin=time.time()
print "Tiempo transcurrido: "+str((t_fin-t_ini))
#Correr el Retroalimentado
t_ini=time.time()
t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std = simulacion_primitiva_retro(num_cells,delta_t,t_limite)
t_fin=time.time()
print "Tiempo transcurrido: "+str((t_fin-t_ini))

# <codecell>

#Graficar
plt.figure(4)
plt.plot(t_final,r_mean)
plt.title("Numero de Rnas promedio vs tiempo")
xlabel("tiempo (min)")
ylabel("Numero de Rnas promedio")

plt.figure(5)
plt.plot(t_final,prot_mean)
plt.title("Numero de Proteinas promedio vs tiempo")
xlabel("tiempo (min)")
ylabel("Numero de Proteinas promedio")

plt.figure(6)
plt.plot(t_final,r_std/r_mean)
plt.title("Ruido para Rna vs tiempo")
xlabel("tiempo (min)")
ylabel("Ruido para Rna")

plt.figure(7)
plt.plot(t_final,prot_std/prot_mean)
plt.title("Ruido para proteina vs tiempo")
xlabel("tiempo (min)")
ylabel("Ruido para proteina")

num_puntos=int(t_limite/delta_t)

plt.figure(8)
hist(matriz_rna[:, num_puntos-1], bins = 20) #cogemos el ultimo
plt.title("Distribucion de RNA en estado estable")
xlabel("numero de Rna")
ylabel("numero de celulas")

plt.figure(9)
hist(matriz_prot[:, num_puntos-1], bins = 60) #cogemos el ultimo
plt.title("Distribucion de proteinas en estado estable")
xlabel("numero de proteinas")
ylabel("numero de celulas")

