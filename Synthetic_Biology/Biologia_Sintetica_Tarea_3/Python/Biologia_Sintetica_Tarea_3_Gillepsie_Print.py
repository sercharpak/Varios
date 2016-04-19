# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline
import numpy as np
import matplotlib.pyplot as plt
import time
#Hecho por Sergio Daniel Hernandez Charpak
#Codigo para el punto 5 de la tarea 3 de Biologia Sintetica en el semestre 2015-1 en Uniandes dictado por Juan Manuel Pedraza
#Algoritmo Gillepsie

# <codecell>

def evento_gillepsie(num_rna,num_prot,t,k_r,k_p,gamma_p,gamma_r):
    e1=np.random.random() #se crea un numero aleatorio entre 0 y 1 mediante una distribucion uniforme 
    #Normalicemos las constantes
    k_total=k_r+k_p*num_rna+gamma_p*num_prot+gamma_r*num_rna
    paso_t=np.random.exponential()/(k_total)
    t_nuevo=t+paso_t
    #Ahora mimremos que tipo de evento sucedio
    if(e1<(k_r/k_total)):#Habemus rna
        num_rna_nuevo=num_rna+1.0
        num_prot_nuevo=num_prot
    elif((e1<((k_r+gamma_r*num_rna)/k_total)) and (e1>=(k_r/k_total))):#Destruimos Rna
        num_rna_nuevo=num_rna-1.0
        num_prot_nuevo=num_prot
    elif((e1<((k_r+k_p*num_rna+gamma_r*num_rna)/k_total)) and (e1>=((k_r+gamma_r*num_rna)/k_total))):#HAbemus proteina
        num_rna_nuevo=num_rna
        num_prot_nuevo=num_prot+1.0
    else:#Destruimos proteina
        num_rna_nuevo=num_rna
        num_prot_nuevo=num_prot-1.0
    return num_rna_nuevo, num_prot_nuevo,t_nuevo

def evento_gillepsie_retro(num_rna,num_prot,t,k_r,K,k_p,gamma_p,gamma_r):
    e1=np.random.random() #se crea un numero aleatorio entre 0 y 1 mediante una distribucion uniforme 
    #Normalicemos las constantes
    k_total=(k_r/(1.0+((num_prot/K)**2.0)))+k_p*num_rna+gamma_p*num_prot+gamma_r*num_rna
    paso_t=np.random.exponential()/(k_total)
    t_nuevo=t+paso_t
    #Ahora mimremos que tipo de evento sucedio
    if(e1<((k_r/(1.0+((num_prot/K)**2.0)))/k_total)):#Habemus rna
        num_rna_nuevo=num_rna+1.0
        num_prot_nuevo=num_prot
    elif((e1<(((k_r/(1.0+((num_prot/K)**2.0)))+gamma_r*num_rna)/k_total)) and (e1>=((k_r/(1.0+((num_prot/K)**2.0)))/k_total))):#Destruimos Rna
        num_rna_nuevo=num_rna-1.0
        num_prot_nuevo=num_prot
    elif((e1<(((k_r/(1.0+((num_prot/K)**2.0)))+k_p*num_rna+gamma_r*num_rna)/k_total)) and (e1>=(((k_r/(1.0+((num_prot/K)**2.0)))+gamma_r*num_rna)/k_total))):#HAbemus proteina
        num_rna_nuevo=num_rna
        num_prot_nuevo=num_prot+1.0
    else:#Destruimos proteina
        num_rna_nuevo=num_rna
        num_prot_nuevo=num_prot-1.0
    return num_rna_nuevo, num_prot_nuevo,t_nuevo

# <codecell>

def simulacion_gillepsie(num_cells,t_limite,delta_t,k_r,k_p,gamma_r,gamma_p):
    contador_cells=0.0
    num_puntos=int(t_limite/delta_t)
    matriz_rna=np.zeros((num_cells,num_puntos))
    matriz_prot=np.zeros((num_cells,num_puntos))
    
    while (contador_cells<num_cells):
        num_rna=0.0
        num_prot=0.0
        t_total=0.0
        t_walk=[]
        r_walk=[]
        prot_walk=[]
        while (t_total<t_limite):
            num_rna, num_prot,t_total=evento_gillepsie(num_rna,num_prot,t_total,k_r,k_p,gamma_p,gamma_r)
            t_walk.append(t_total)
            r_walk.append(num_rna)
            prot_walk.append(num_prot)
        #Ahora se deben reformar los arreglos por los tiempos
        t_final = np.linspace(0, t_limite-delta_t, t_limite/delta_t)
        r_final = np.zeros(len(t_final))
        prot_final = np.zeros(len(t_final))
        pos=0
        for i in range (len(t_final)):
            while((t_walk[pos] < t_final[i])&((pos+1)<len(prot_walk))):
                pos=pos+1
            r_final[i] = r_walk[pos]
            prot_final[i] = prot_walk[pos]
        matriz_rna[contador_cells]=r_final
        matriz_prot[contador_cells]=prot_final
        contador_cells=contador_cells+1.0
    
    #Ahora calculamos la media y la desviacion
    r_mean=np.zeros(num_puntos)
    r_std=np.zeros(num_puntos)
    prot_mean=np.zeros(num_puntos)
    prot_std=np.zeros(num_puntos)
    for i in range(num_puntos):
        r_mean[i]=np.mean(matriz_rna[:,i])
        r_std[i]=np.std(matriz_rna[:,i])
        prot_mean[i]=np.mean(matriz_prot[:,i])
        prot_std[i]=np.std(matriz_prot[:,i])
        
        
    return t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std

def simulacion_gillepsie_retro(num_cells,t_limite,delta_t,k_r,K,k_p,gamma_r,gamma_p):
    contador_cells=0.0
    num_puntos=int(t_limite/delta_t)
    matriz_rna=np.zeros((num_cells,num_puntos))
    matriz_prot=np.zeros((num_cells,num_puntos))
    
    while (contador_cells<num_cells):
        num_rna=0.0
        num_prot=0.0
        t_total=0.0
        t_walk=[]
        r_walk=[]
        prot_walk=[]
        while (t_total<t_limite):
            num_rna, num_prot,t_total=evento_gillepsie_retro(num_rna,num_prot,t_total,k_r,K,k_p,gamma_p,gamma_r)
            t_walk.append(t_total)
            r_walk.append(num_rna)
            prot_walk.append(num_prot)
        #Ahora se deben reformar los arreglos por los tiempos
        t_final = np.linspace(0, t_limite-delta_t, t_limite/delta_t)
        r_final = np.zeros(len(t_final))
        prot_final = np.zeros(len(t_final))
        pos=0
        for i in range (len(t_final)):
            while((t_walk[pos] < t_final[i])&((pos+1)<len(prot_walk))):
                pos=pos+1
            r_final[i] = r_walk[pos]
            prot_final[i] = prot_walk[pos]
        matriz_rna[contador_cells]=r_final
        matriz_prot[contador_cells]=prot_final
        contador_cells=contador_cells+1.0
    
    #Ahora calculamos la media y la desviacion
    r_mean=np.zeros(num_puntos)
    r_std=np.zeros(num_puntos)
    prot_mean=np.zeros(num_puntos)
    prot_std=np.zeros(num_puntos)
    for i in range(num_puntos):
        r_mean[i]=np.mean(matriz_rna[:,i])
        r_std[i]=np.std(matriz_rna[:,i])
        prot_mean[i]=np.mean(matriz_prot[:,i])
        prot_std[i]=np.std(matriz_prot[:,i])
        
        
    return t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std

# <codecell>

#Constantes biologicamente relevantes
k_r=1.0 #rna/min
k_p=50.0 # 1 prot/(min*rna)
K=0.5 #cte de retroalimentacion
gamma_p=1.0/30.0 #  prot/min
gamma_r=1.0/5.0 #  rna/min

# <codecell>

#Correr el no Retroalimentado
t_ini=time.time()
t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std = simulacion_gillepsie(num_cells,t_limite,delta_t,k_r,k_p,gamma_r,gamma_p)
t_fin=time.time()
print "Tiempo transcurrido: "+str((t_fin-t_ini))
#Correr el Retroalimentado: 
t_ini=time.time()
t_final, matriz_rna, matriz_prot, r_mean, r_std,prot_mean,prot_std = simulacion_gillepsie_retro(num_cells,t_limite,delta_t,k_r,K,k_p,gamma_r,gamma_p)
t_fin=time.time()
print "Tiempo transcurrido: "+str((t_fin-t_ini))

# <codecell>

#Graficar
plt.figure(4)
plt.plot(t_final,r_mean)
plt.title("Numero de Rnas promedio vs tiempo con retroalimentacion")
xlabel("tiempo (min)")
ylabel("Numero de Rnas promedio")

plt.figure(5)
plt.plot(t_final,prot_mean)
plt.title("Numero de Proteinas promedio vs tiempo con retroalimentacion")
xlabel("tiempo (min)")
ylabel("Numero de Proteinas promedio")

plt.figure(6)
plt.plot(t_final,r_std/r_mean)
plt.title("Ruido para Rna vs tiempo con retroalimentacion")
xlabel("tiempo (min)")
ylabel("Ruido para Rna")

plt.figure(7)
plt.plot(t_final,prot_std/prot_mean)
plt.title("Ruido para proteina vs tiempo con retroalimentacion")
xlabel("tiempo (min)")
ylabel("Ruido para proteina")

num_puntos=int(t_limite/delta_t)

plt.figure(8)
hist(matriz_rna[:, num_puntos-1], bins = 20) #cogemos el ultimo
plt.title("Distribucion de RNA en estado estable con retroalimentacion")
xlabel("numero de Rna")
ylabel("numero de celulas")

plt.figure(9)
hist(matriz_prot[:, num_puntos-1], bins = 60) #cogemos el ultimo
plt.title("Distribucion de proteinas en estado estable con retroalimentacion")
xlabel("numero de proteinas")
ylabel("numero de celulas")

