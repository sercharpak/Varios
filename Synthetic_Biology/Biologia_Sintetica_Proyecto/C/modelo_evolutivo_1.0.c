#include <stdio.h>            // Needed for printf()
#include <stdlib.h>           // Needed for exit() and ato*()
#include <math.h>             // Needed for log()
#include <string.h>

double aleatorio_exponencial(int semilla);
double *evento_gillepsie(double num_vec_Ecol,double num_vec_Salm,double num_vec_Rhod,double num_vec_Cont, double t, double fit_Ecol, double fit_Salm,double fit_Rhod, double fit_Cont, int semilla);
float *crear_matriz(int n, int m);
double *simulacion_gillepsie(int num_corridas,double t_limite, double delta_t, double fit_Ecol, double fit_Salm, double fit_Rhod, double fit_Cont);
float*crearArregloCero(int n_points);
float*crearArregloEquiEspaciado(float x_ini,float x_fin,int n_points);
float*copiarArreglo(float* original, int n_points);


int main(int argc, char **argv){
double x;
int i;
for(i=0;i<10;i++){
x=aleatorio_exponencial(i);
printf("%f\n",x);
}
return 0;
}

double aleatorio_exponencial(int semilla){
double random;
double random_exponencial;
    srand48(semilla);
    random = drand48();
random_exponencial=exp(-random);
return random_exponencial;
}

double *evento_gillepsie(double num_vec_Ecol,double num_vec_Salm,double num_vec_Rhod,double num_vec_Cont, double t, double fit_Ecol, double fit_Salm,double fit_Rhod, double fit_Cont, int semilla){
double rand_exp, e1;
double k_total;
double paso_t;
double *arreglo_salida;
arreglo_salida = malloc(5* sizeof(double));
rand_exp=aleatorio_exponencial(semilla);
e1=drand48();
k_total=fit_Ecol+fit_Salm+fit_Rhod+fit_Cont;
paso_t=rand_exp/(k_total);
arreglo_salida[4]=t+paso_t;

//Ahora mimremos que tipo de evento sucedio
 if(e1<(fit_Ecol/k_total)){//Habemus Vec_E_Coli
        arreglo_salida[0]=num_vec_Ecol+1.0;
        arreglo_salida[1]=num_vec_Salm;
        arreglo_salida[2]=num_vec_Rhod;
        arreglo_salida[3]=num_vec_Cont;
}
else if((e1<((fit_Ecol+fit_Salm)/k_total)) && (e1>=(fit_Ecol/k_total))){//Habemus Vec_Salm
        arreglo_salida[0]=num_vec_Ecol;
	arreglo_salida[1]=num_vec_Salm+1.0;      
        arreglo_salida[2]=num_vec_Rhod;
        arreglo_salida[3]=num_vec_Cont;
}
else if((e1<((fit_Ecol+fit_Salm+fit_Rhod)/k_total)) && (e1>=((fit_Ecol+fit_Salm)/k_total))){//Habemus Vec_Rhodn
        arreglo_salida[0]=num_vec_Ecol;
        arreglo_salida[1]=num_vec_Salm;
        arreglo_salida[2]=num_vec_Rhod+1.0;
        arreglo_salida[3]=num_vec_Cont;
}
    else{//Habemus Vec_Cont
        arreglo_salida[0]=num_vec_Ecol;
        arreglo_salida[1]=num_vec_Salm;
        arreglo_salida[2]=num_vec_Rhod;
        arreglo_salida[3]=num_vec_Cont+1.0;
}
    return arreglo_salida;
}


double* simulacion_gillepsie(num_corridas,t_limite,delta_t,fit_Ecol,fit_Salm,fit_Rhod,fit_Cont){
int contador_corridas, contador_gillepsie, pos, i, j,k;
double num_puntos, t_total,tf;
double *matriz_vec_Ecol;
double *matriz_vec_Salm;
double *matriz_vec_Rhod;
double *matriz_vec_Cont;
double *vec_Ecol_walk;
double *vec_Salm_walk;
double *vec_Rhod_walk;
double *vec_Cont_walk;
double *t_walk;
double *t_final;
double *vec_Ecol_final;
double *vec_Salm_final;
double *vec_Rhod_final;
double *vec_Cont_final;
double *arreglo_gillepsie;
double num_vec_Ecol, num_vec_Salm, num_vec_Rhod, num_vec_Cont;
    contador_corridas=0;
    num_puntos=(int) (t_limite/delta_t);
    matriz_vec_Ecol=crear_matriz(num_corridas,num_puntos);
    matriz_vec_Salm=crear_matriz(num_corridas,num_puntos);
    matriz_vec_Rhod=crear_matriz(num_corridas,num_puntos);
    matriz_vec_Cont=crear_matriz(num_corridas,num_puntos);
    tf=t_limite-delta_t;

    while (contador_corridas<num_corridas){
        num_vec_Ecol=45;
        num_vec_Salm=45;
        num_vec_Rhod=45;
        num_vec_Cont=45;
        t_total=0.0;
	t_final=crearArregloEquiEspaciado(0,tf,num_puntos);
        t_walk=crearArregloCero(num_puntos);
        vec_Ecol_walk=crearArregloCero(num_puntos+100);
        vec_Salm_walk=crearArregloCero(num_puntos+100);
        vec_Rhod_walk=crearArregloCero(num_puntos+100);
        vec_Cont_walk=crearArregloCero(num_puntos+100);
	contador_gillepsie=0;
        while (t_total<t_limite){
            arreglo_gillepsie=evento_gillepsie(num_vec_Ecol,num_vec_Salm,num_vec_Rhod,num_vec_Cont,t_total,fit_Ecol,fit_Salm,fit_Rhod,fit_Cont,contador_gillepsie);
	    num_vec_Ecol=arreglo_gillepsie[0];
	    num_vec_Salm=arreglo_gillepsie[1];
	    num_vec_Rhod=arreglo_gillepsie[2];
	    num_vec_Cont=arreglo_gillepsie[3];
	    t_total=arreglo_gillepsie[4];
            t_walk[contador_gillepsie]=t_total;
            vec_Ecol_walk[contador_gillepsie]=num_vec_Ecol;
            vec_Salm_walk[contador_gillepsie]=num_vec_Salm;
            vec_Rhod_walk[contador_gillepsie]=num_vec_Rhod;
            vec_Cont_walk[contador_gillepsie]=num_vec_Cont;
	    contador_gillepsie++;
	}
        //Ahora se deben reformar los arreglos por los tiempos
        vec_Ecol_final =crearArregloCero(num_puntos);
        vec_Salm_final =crearArregloCero(num_puntos);
        vec_Rhod_final =crearArregloCero(num_puntos);
        vec_Cont_final =crearArregloCero(num_puntos);
        pos=0;
        for (i=0;i<num_puntos+100;i++){
            while((t_walk[pos] < t_final[i])&((pos+1)<(num_puntos+100))){
                pos=pos+1;
		}
            vec_Ecol_final[i] = vec_Ecol_walk[pos];
            vec_Salm_final[i] = vec_Salm_walk[pos];
            vec_Rhod_final[i] = vec_Rhod_walk[pos];
            vec_Cont_final[i] = vec_Cont_walk[pos];
	}
	for(k=0;k<num_puntos;k++){
	matriz_vec_Ecol[contador_corridas*num_puntos+k]=vec_Ecol_final[k];
	matriz_vec_Salm[contador_corridas*num_puntos+k]=vec_Salm_final[k];
        matriz_vec_Rhod[contador_corridas*num_puntos+k]=vec_Rhod_final[k];
        matriz_vec_Cont[contador_corridas*num_puntos+k]=vec_Cont_final[k];
	}        
        contador_corridas=contador_corridas+1;
    }
/**
    //Ahora calculamos la media y la desviacion
    vec_Ecol_mean=crearArregloCero(num_puntos);
    vec_Ecol_std=crearArregloCero(num_puntos);
    vec_Salm_mean=crearArregloCero(num_puntos);
    vec_Salm_std=crearArregloCero(num_puntos);
    vec_Rhod_mean=crearArregloCero(num_puntos);
    vec_Rhod_std=crearArregloCero(num_puntos);
    vec_Cont_mean=crearArregloCero(num_puntos);
    vec_Cont_std=crearArregloCero(num_puntos);
    for(j=0;j<num_puntos;j++){
        vec_Ecol_mean[j]=np.mean(matriz_vec_Ecol[:,i])
        vec_Ecol_std[j]=np.std(matriz_vec_Ecol[:,i])
        vec_Salm_mean[j]=np.mean(matriz_vec_Salm[:,i])
        vec_Salm_std[i]=np.std(matriz_vec_Salm[:,i])
        vec_Rhod_mean[i]=np.mean(matriz_vec_Rhod[:,i])
        vec_Rhod_std[i]=np.std(matriz_vec_Rhod[:,i])
        vec_Cont_mean[i]=np.mean(matriz_vec_Cont[:,i])
        vec_Cont_std[i]=np.std(matriz_vec_Cont[:,i])
        }
    return t_final, matriz_vec_Ecol, matriz_vec_Salm, matriz_vec_Rhod, matriz_vec_Cont, vec_Ecol_mean, vec_Ecol_std, vec_Salm_mean,vec_Salm_std,vec_Rhod_mean,vec_Rhod_std,vec_Cont_mean,vec_Cont_std
*/
}

/*
 *Funcion que dado un arreglo y su tamano, crea y retorna una copia de este.
 */
float*copiarArreglo(float* original, int n_points)
{
  int i;
  float* arregloRespuesta;
  arregloRespuesta=malloc(n_points * sizeof(float));
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=original[i];
    }
  return arregloRespuesta;
}

/*
 * Metodo que crea una matriz n*m con 0 en todas las entradas.
 * @params n numero de filas
 * @params m numero de columnas
 * @returns *float matriz de ceros n*m
 */
float *crear_matriz( int n, int m){
  float *matrix;
  int i;
  int j;
  matrix = malloc(n * m * sizeof(float));

  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
       matrix[i*m + j]=0.0;
    }
  }    
  return matrix;
}


/*
 *Funcion que crea y retorna un vector de n_points puntos equiespaciados entre x_ini y x_fin
 */
float*crearArregloEquiEspaciado(float x_ini,float x_fin,int n_points)
{
  int i;
  float* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(float));
  float paso;
  paso=(x_fin-x_ini)/n_points;
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=0.0+paso*i;
    }
  return arregloRespuesta;
}


/*
 *Funcion que crea y retorna un vector de n_points posiciones con 0.0 en ellas.
 */
float*crearArregloCero(int n_points)
{
  int i;
  float* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(float));
  //Inicialicemoslo en 0
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=0.0;
    }
  return arregloRespuesta;
}
