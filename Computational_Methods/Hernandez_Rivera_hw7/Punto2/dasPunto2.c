//Creado para trabajar al punto 2 de la tarea 7 por Sergio
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define gamma 1.4
/*
 *EJECUTAR USANDO MAKEFILE
 *Los parametros se ingresan en el siguiente orden: tf
 */
double*crearArregloCero(int n_points);
double*crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points);
void escribirArreglos(double*x,double*rho, double*vel, double*p,int n_x, double tf);
double*copiarArreglo(double* original, int n_points);
double*copiarMatriz(double* original, int n, int m);
double *crear_matriz(int n, int m);
double *suma(double *matrix_A, double *matrix_B, int n, int m);
double *multiplicacionPorEscalar(double *matrix_A, double escalar, int n, int m);
double calcularEnergia(double p, double rho, double vel);
double*calcularF(double*u, int m);
double*formarMatrizUInicial(double*x,double et_i_l, double rho_i_l, double vel_i_l, double et_i_r, double rho_i_r, double vel_i_r, int m);
double*calcularUEstadoTf(double*u_init,double*f_init,int n_t,int n_x, double dt, double dx);
double*calcularUMedios(double*u_init,double*f_init,int n_x,double dt, double dx);
double*obtenerRho(double*u_futura,int n_x);
double*obtenerVel(double*u_futura,int n_x);
double*obtenerP(double*u_futura,int n_x);
int  main(int argc, char **argv)
{
  double tf,ti,dt,xi,xf,dx,et_i_l, rho_i_l, vel_i_l,et_i_r,rho_i_r,vel_i_r,p_i_l,p_i_r;
  int n_x,n_t,n,m;
  double*t,*x,*rho,*p,*vel;
  double*u_init,*f_init,*u_futura;
  tf=atof(argv[1]);
  ti=0.0;
  dt=0.00001;
  n_t=(int)round(tf/dt);
  xi=-10.0;
  xf=10.0;
  n_x=1000;
  x=crearArregloEquiEspaciado(xi,xf,n_x);
  dx=x[1]-x[0];
  rho_i_l=1.0;
  rho_i_r=0.125;
  vel_i_l=0.0;
  vel_i_r=0.0;
  p_i_l=100.0*pow(10.0,3.0);
  p_i_r=10.0*pow(10.0,3.0);
  et_i_l=calcularEnergia(p_i_l,rho_i_l,vel_i_l);
  et_i_r=calcularEnergia(p_i_r,rho_i_r,vel_i_r);
  u_init=formarMatrizUInicial(x,et_i_l,rho_i_l,vel_i_l,et_i_r,rho_i_r,vel_i_r, n_x);
  f_init=calcularF(u_init,n_x);
  /* u_futura=calcularUMedios(u_init,f_init,n_x,dt,dx); */
  u_futura=calcularUEstadoTf(u_init,f_init,n_t,n_x,dt,dx);
  rho=obtenerRho(u_futura,n_x);
  p=obtenerP(u_futura,n_x);
  vel=obtenerVel(u_futura,n_x);
  escribirArreglos(x,rho,vel,p,n_x,tf);
  return 0;
}
/**
 * Metodo que calcular el arreglo rho dada la matriz de estado u_futura
 **/
double*obtenerRho(double*u_futura,int n_x)
{
  double*rho;
  int j;
  rho=crearArregloCero(n_x);
  for(j=0;j<n_x;j++){
    rho[j]=u_futura[j];
  }
  return rho;
}
/**
 * Metodo que calcular el arreglo vel dada la matriz de estado u_futura
 **/
double*obtenerVel(double*u_futura,int n_x)
{
  double*vel;
  int j;
  vel=crearArregloCero(n_x);
  for(j=0;j<n_x;j++){
    vel[j]=u_futura[1*n_x+j]/u_futura[j];
  }
  return vel;
}
/**
 * Metodo que calcular el arreglo p dada la matriz de estado u_futura
 **/
double*obtenerP(double*u_futura,int n_x)
{
  double*p;
  int j;
  p=crearArregloCero(n_x);
  for(j=0;j<n_x;j++){
    p[j]=(gamma-1)*(u_futura[2*n_x+j]-(0.5)*(u_futura[1*n_x+j]*u_futura[1*n_x+j]/u_futura[j]));
  }
  return p;
}
/**
 * Metodo que calcular UMedios dados los parametros necesarios
 **/
double*calcularUMedios(double*u_init,double*f_init,int n_x,double dt, double dx)
{
  double*u_medios;
  int i,j,n,m;
  n=3;
  m=n_x;
  u_medios=copiarMatriz(u_init,n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m-1;j++){
      u_medios[i*m+j]=0.5*(u_init[i*m+(j+1)]+u_init[i*m+j])-(dt/(2.0*dx))*(f_init[i*m+(j+1)]-f_init[i*m+j]);
    }
  }
  return u_medios;
}
/**
 * Metodo que calcula la matriz u hasta el estado tf que correspondo a ejecutar n_t-1 veces 
 * el loop sobre el espacio n_x
 **/
double*calcularUEstadoTf(double*u_init,double*f_init,int n_t,int n_x,double dt, double dx)
{
  double*u_futura,*f_futura,*u_medios,*f_medios;
  int i,j,k,h;
  u_futura=crear_matriz(3,n_x);
  for(i=0;i<3;i++){
    u_futura[i*n_x+0]=u_init[i*n_x+0];
    u_futura[i*n_x+n_x-1]=u_init[i*n_x+n_x-1];
  }
  for(k=1;k<n_t;k++){
    u_medios=calcularUMedios(u_init,f_init,n_x,dt,dx);
    f_medios=calcularF(u_medios,n_x);
    for(i=0;i<3;i++){
      for(j=1;j<n_x-1;j++){
	u_futura[i*n_x+j]=u_init[i*n_x+j]-(dt/dx)*(f_medios[i*n_x+j]-f_medios[i*n_x+j-1]);
      }
    }
    f_futura=calcularF(u_futura,n_x);
    u_init=copiarMatriz(u_futura,3,n_x);
    f_init=copiarMatriz(f_futura,3,n_x);
  }
  
  return u_futura;
}

/**
 * Metodo que forma la matriz de los valores iniciales de la matriz u
 **/
double*formarMatrizUInicial(double*x,double et_i_l, double rho_i_l, double vel_i_l, double et_i_r, double rho_i_r, double vel_i_r, int m)
{
  int i,j;
  double* matriz_u_i;
  matriz_u_i=crear_matriz(3,m);
  for(j=0;j<m;j++)
    {
      if(x[j]<0.0)
	{
	  matriz_u_i[0*m+j]=rho_i_l;
	  matriz_u_i[1*m+j]=rho_i_l*vel_i_l;
	  matriz_u_i[2*m+j]=rho_i_l*et_i_l;
	}
      else
	{
	  matriz_u_i[0*m+j]=rho_i_r;
	  matriz_u_i[1*m+j]=rho_i_r*vel_i_r;
	  matriz_u_i[2*m+j]=rho_i_r*et_i_r;
	}
    }
  return matriz_u_i;
}

/**
 *Metodo que calcula la matriz F dada la matriz u respectiva y su tamano (3*m)
 **/
double*calcularF(double*u,int m)
{
  /*
    f[0]=u[1]
    f[1]=((u[1]*u[1])/u[0])+(gamma-1)*(u[2]-((1.0/2.0)*((u[1]*u[1])/u[0])))
    f[2]=(u[2]+(gamma-1)*(u[2]-((1.0/2.0)*((u[1]*u[1])/u[0]))))*((u[1])/u[0])
  */
  int i,j;
  double*matrizF;
  matrizF=crear_matriz(3,m);
  for(j=0;j<m;j++){
    matrizF[0*m+j]=u[1*m+j];
    matrizF[m+j]=((u[m+j]*u[m+j])/u[j])+(gamma-1.0)*(u[2*m+j]-((0.5)*((u[m+j]*u[m+j])/u[j])));
    matrizF[2*m+j]=(u[2*m+j]+(gamma-1.0)*(u[2*m+j]-((1.0/2.0)*((u[1*m+j]*u[1*m+j])/u[0*m+j]))))*((u[1*m+j])/u[0*m+j]);
  }
  return matrizF;
}

/**
 *Metodo que calcula la energia dada la presion, la densidad y la velocidad
 **/
double calcularEnergia(double p, double rho, double vel)
{
  return (1.0/rho)*(p/(gamma-1.0))+(0.5*vel*vel);
}
/*
 * Metodo que suma dos matriz y pone el resultado en la tercera.
 * @params matricita1, matricita2 las dos matrices a sumar
 * @params matricitaResultado en donde se imprime el resultado
 * @params n numero de filas, m numero de columnas
 */
double *multiplicacionPorEscalar(double *matrix_A, double escalar, int n, int m)
{
  double *matricitaResultado=crear_matriz( n, m);
  int i,j,k;
  for (i=0;i<n;i++)
    {
    for(j=0;j<m;j++)
      {
	matricitaResultado[i*m +j ]=0.0; // Ponemos toda la matriz en ceros para evitar problemas
	matricitaResultado[i*m +j ]=matrix_A[i*m + j]*escalar;
      }
    }
  return matricitaResultado;
}
/*
 * Metodo que suma dos matriz y pone el resultado en la tercera.
 * @params matricita1, matricita2 las dos matrices a sumar
 * @params matricitaResultado en donde se imprime el resultado
 * @params n numero de filas, m numero de columnas
 */
double *suma(double *matricita1, double *matricita2,  int n, int m)
{
  double *matricitaResultado=crear_matriz( n, m);
  int i,j,k;
  for (i=0;i<n;i++)
    {
    for(j=0;j<m;j++)
      {
	matricitaResultado[i*m +j ]=0.0; // Ponemos toda la matriz en ceros para evitar problemas
	matricitaResultado[i*m +j ]=matricita1[i*m+j]+matricita2[i*m+j];
      }
    }
  return matricitaResultado;
}

/*
 * Metodo que crea una matriz n*m con 0 en todas las entradas.
 * @params n numero de filas
 * @params m numero de columnas
 * @returns *double matriz de ceros n*m
 */
double *crear_matriz( int n, int m){
  double *matrix;
  int i;
  int j;
  matrix = malloc(n * m * sizeof(double));

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
double*crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(double));
  double paso;
  paso=(x_fin-x_ini)/n_points;
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=x_ini+paso*i;
    }
  return arregloRespuesta;
}

/*
 *Funcion que crea y retorna un vector de n_points posiciones con 0.0 en ellas.
 */
double*crearArregloCero(int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta= malloc(n_points * sizeof(double));
  //Inicialicemoslo en 0
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=0.0;
    }
  return arregloRespuesta;
}
/*
 *Funcion que dado un arreglo y su tamano, crea y retorna una copia de este.
 */
double*copiarArreglo(double*original, int n_points)
{
  int i;
  double* arregloRespuesta;
  arregloRespuesta=malloc(n_points * sizeof(double));
  for(i=0;i<n_points;i++)
    {
      arregloRespuesta[i]=original[i];
    }
  return arregloRespuesta;
}
/*
 *Funcion que dada una matriz y su tamano, crea y retorna una copia de esta.
 */
double*copiarMatriz(double*original, int n, int m)
{
double *matricitaResultado=crear_matriz( n, m);
 int i,j;
  for (i=0;i<n;i++)
    {
    for(j=0;j<m;j++)
      {
	matricitaResultado[i*m+j]=0.0;//Ponemos la matriz en ceros para evitar problemas
	matricitaResultado[i*m+j]=original[i*m+j];
      }
    }
  return matricitaResultado;
}
/*
 *Funcion que escribe 4 arreglos TODOS de tamano n_x  en el tiempo t en un archivo .dat (dado su tamano).
 */
void escribirArreglos(double*x,double*rho, double*vel, double*p,int n_x, double tf)
{
  FILE* archivo;
  //Sabemos que solo hay 4 arreglos

  char bufX[20];
  int i;
  sprintf(bufX, "%f", tf);
  char n1[50], n3[50];
   strcpy(n1,  "estado_");
   strcpy(n3, ".dat");
   strcat(n1, bufX);
   strcat(n1, n3);
   archivo = fopen(n1, "w");
  for(i=0;i<n_x;i++){
    fprintf(archivo, "%f \t %f \t %f \t %f \n", x[i], vel[i], p[i], rho[i]);
  }
  fclose(archivo);
}
