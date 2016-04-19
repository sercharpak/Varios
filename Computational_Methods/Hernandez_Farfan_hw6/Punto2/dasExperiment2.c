#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define masa 1.6726e-27
#define c 3.0e8
#define b0 3.0e-5
#define rT 6.3781e6
#define e 1.6e-19
#define pi 3.14159265359

double *crear_matriz(int n, int m);
double *velocidadesIniciales(double energia, double alpha);
double *posicionesIniciales();
double magnitudVelocidadInicial(double energia);
double magnitudVector(double*vector,int n);
void escribirArreglos(double*t,double*x,double*y, double*z,int n_x,double energio, double alpha);
void rungerKuttaMagnetico(double* x, double* y, double*z, double*velI, double* t, int n_points, double energia, double alpha,double gamma);
double* pasoRungerKutta4(double paso,double* posiciones,double* velocidades,double gamma);
double*calcularAceleraciones(double*posiciones,double*velocidades,double gamma);
double*crearArregloCero(int n_points);
double*crearArregloEquiEspaciado(double x_ini,double x_fin,int n_points);


int  main(int argc, char **argv)
{
  double energia_calc,energia,alpha,alpha_calc,magnitudVelocidadInicia,magnitudPosicion,ti,tf,v,gamma;
  double *velI,*posI,*x,*y,*z,*t;
  double *matrizAResolver,*matricitaL,*matricitaLT,*matricitaY,*matricitaZ,*matricitaX;
  int orden,n_points,k,h;
  int n_filas;
  n_points=1000000;
  ti=0.0;
  tf=100.0;
  n_filas=3;
  orden =2;
  energia=atof(argv[1]);
  energia_calc=e*pow(10.0,6)*atof(argv[1]);// LA ENERGIO ESTA EN MeV!!!!
  alpha=atof(argv[2]);
  alpha_calc = atof(argv[2])*pi/180.0; 
  velI=velocidadesIniciales(energia_calc,alpha_calc);
  v=magnitudVector(velI,3);
  //Freaking Gamma
  // gamma=1.0/(sqrt(1.0-(v*v)/(c*c)));
  gamma=energia_calc/(masa*c*c)+1.0;
  posI=posicionesIniciales();
  t=crearArregloEquiEspaciado(ti,tf,n_points);
  x=crearArregloCero(n_points);
  x[0]=posI[0];
  y=crearArregloCero(n_points);
  y[0]=posI[1];
  z=crearArregloCero(n_points);
  z[0]=posI[2];
  //RungerKutta Time!
  rungerKuttaMagnetico(x,y,z,velI,t,n_points,energia,alpha,gamma);
  return 0;
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
      arregloRespuesta[i]=0.0+paso*i;
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
 *Funcion que calcula las aceleraciones del sistema en el tiempo t dadas las posiciones y las velocidades en el mismo tiempo t
 *Esta es la funcion que realiza el analisis matricial (despomposicion de cholesky)
 */
double*calcularAceleraciones(double*pos,double*vel,double gamma)
{
  double dasConstante;
  double r;
  double ax,ay,az;
  double*acc=crear_matriz(3,1);
  r=magnitudVector(pos,3);
  dasConstante=-(b0*rT*rT*rT*e)/(gamma*masa);
  ax=dasConstante*((2.0*pos[2]*pos[2]-pos[0]*pos[0]-pos[1]*pos[1])*vel[1]-(3.0*pos[1]*pos[2])*vel[2])/(pow(r,5));
ay=dasConstante*((3.0*pos[0]*pos[2])*vel[2]-(2.0*pos[2]*pos[2]-pos[0]*pos[0]-pos[1]*pos[1])*vel[0])/(pow(r,5));
az=dasConstante*((3.0*pos[1]*pos[2])*vel[0]-(3.0*pos[0]*pos[2])*vel[1])/(pow(r,5));
acc[0]=ax;
acc[1]=ay;
acc[2]=az;
return acc;

}

    /*
     *Funcion que dado un paso (en el tiempo), un x y un y (al tiempo t) calcula el x, y y z al tiempo t+paso usando Runger Kutta 4to orden
     *Devuelve una matriz 3x2 con las posiciones en la primera columna y las velocidades en la segunda
     */
double* pasoRungerKutta4(double paso,double* posiciones,double* velocidades,double gamma)
    {
      //x(t+paso) = x+paso*(k1+2*k2+2*k3+k4)/6;
      //vx(t+paso) = vx+paso*(kv1+2*kv2+2*kv3+kv4)/6;
      //y(t+paso) = y+paso*(m1+2*m2+2*m3+m4)/6;
      //vy(t+paso) = vy+paso*(mv1+2*mv2+2*mv3+mv4)/6;
      //z(t+paso) = z+paso*(n1+2*n2+2*n3+n4)/6;
      //vz(t+paso) =v z+paso*(nv1+2*nv2+2*nv3+nv4)/6;
      double k1,k2,k3,k4,kF; //Para x
      double kv1,kv2,kv3,kv4,kvF; //Para vx   
      double ka1,ka2,ka3,ka4; //Para ax 
      double m1,m2,m3,m4,mF; //Para y
      double mv1,mv2,mv3,mv4,mvF; //Para vy
      double ma1,ma2,ma3,ma4; //Para ay
      double n1,n2,n3,n4,nF; //Para z
      double nv1,nv2,nv3,nv4,nvF; //Para vz
      double na1,na2,na3,na4; //Para az
      double*pos2,*pos3,*pos4,*posF;
      double*vel2,*vel3,*vel4,*velF;
      double*acc1,*acc2,*acc3,*acc4,*retorno;
      int i,j;
      //Primero - Comienzo Intervalo
      k1=posiciones[0];
      m1=posiciones[1];
      n1=posiciones[2];
      kv1=velocidades[0];
      mv1=velocidades[1];
      nv1=velocidades[2];
      acc1=calcularAceleraciones(posiciones,velocidades,gamma);
      ka1=acc1[0];
      ma1=acc1[1];
      na1=acc1[2];
      //Segundo - Mitad Intervalo
      k2=posiciones[0]+paso*0.5*kv1;
      m2=posiciones[1]+paso*0.5*mv1;
      n2=posiciones[2]+paso*0.5*nv1;
      kv2=velocidades[0]+paso*0.5*ka1;
      mv2=velocidades[1]+paso*0.5*ma1;
      nv2=velocidades[2]+paso*0.5*na1;
      pos2=crear_matriz(3,1);
      vel2=crear_matriz(3,1);
      pos2[0]=k2;
      pos2[1]=m2;
      pos2[2]=n2;
      vel2[0]=kv2;
      vel2[1]=mv2;
      vel2[2]=nv2;
      acc2=calcularAceleraciones(pos2,vel2,gamma);
      ka2=acc2[0];
      ma2=acc2[1];
      na2=acc2[2];
      //Tercero - Mitad Intervalo
      k3=posiciones[0]+paso*0.5*kv2;
      m3=posiciones[1]+paso*0.5*mv2;
      n3=posiciones[2]+paso*0.5*nv2;
      kv3=velocidades[0]+paso*0.5*ka2;
      mv3=velocidades[1]+paso*0.5*ma2;
      nv3=velocidades[2]+paso*0.5*na2;
      pos3=crear_matriz(3,1);
      vel3=crear_matriz(3,1);
      pos3[0]=k3;
      pos3[1]=m3;
      pos3[2]=n3;
      vel3[0]=kv3;
      vel3[1]=mv3;
      vel3[2]=nv3;
      acc3=calcularAceleraciones(pos3,vel3,gamma);
      ka3=acc3[0];
      ma3=acc3[1];
      na3=acc3[2];
      //Cuarto - Final Intervalo
      k4=posiciones[0]+paso*kv3;
      m4=posiciones[1]+paso*mv3;
      n4=posiciones[2]+paso*nv3;
      kv4=velocidades[0]+paso*ka3;
      mv4=velocidades[1]+paso*ma3;
      nv4=velocidades[2]+paso*na3;
      pos4=crear_matriz(3,1);
      vel4=crear_matriz(3,1);
      pos4[0]=k4;
      pos4[1]=m4;
      pos4[2]=n4;
      vel4[0]=kv4;
      vel4[1]=mv4;
      vel4[2]=nv4;
      acc4=calcularAceleraciones(pos4,vel4,gamma);
      ka4=acc4[0];
      ma4=acc4[1];
      na4=acc4[2];
      //Final - Promediemos esta @#%
      kF=posiciones[0]+paso*(kv1+2.0*kv2+2.0*kv3+kv4)/6.0;
      mF=posiciones[1]+paso*(mv1+2.0*mv2+2.0*mv3+nv4)/6.0;
      nF=posiciones[2]+paso*(nv1+2.0*nv2+2.0*nv3+mv4)/6.0;
      kvF=velocidades[0]+paso*(ka1+2.0*ka2+2.0*ka3+ka4)/6.0;
      mvF=velocidades[1]+paso*(ma1+2.0*ma2+2.0*ma3+na4)/6.0;
      nvF=velocidades[2]+paso*(na1+2.0*na2+2.0*na3+ma4)/6.0;
      posF=crear_matriz(3,1);
      velF=crear_matriz(3,1); 
      posF[0]=kF;
      posF[1]=mF;
      posF[2]=nF;
      velF[0]=kvF;
      velF[1]=mvF;
      velF[2]=nvF;
      retorno = crear_matriz(6,1);
      for(i=0;i<3;i++)
	{
	  retorno[i]=posF[i];
	  retorno[i+3]=velF[i];
	}
      return retorno;
    }

  /*
   *Funcion que Ejecuta RungerKutta para el sistema magnetico
   */
void rungerKuttaMagnetico(double* x, double* y, double*z, double*velI, double* t, int n_points, double energia, double alpha,double gamma)
  {
    int i,j;
    double paso;
    //Que el paso es el periodo de ciclotron
    paso=2*pi*gamma*masa/(50.0*e*b0);
    double* velActual;
    velActual=crear_matriz(3,1);
    velActual[0]=velI[0];
    velActual[1]=velI[1];
    velActual[2]=velI[2];
    for(i=0;i<n_points;i++)
      {
	int k;
	k=i+1;
	double* posActual;
	posActual=crear_matriz(3,1);
	posActual[0]=x[i];
	posActual[1]=y[i];
	posActual[2]=z[i];
	double* temp;
	temp=pasoRungerKutta4(paso, posActual,velActual,gamma);
	x[k]=temp[0];
	y[k]=temp[1];
	z[k]=temp[2];
	for(j=0;j<3;j++)
	  {
	    velActual[j]=temp[j+3];
	  }
      }
    escribirArreglos(t,x,y,z,n_points,energia,alpha);
  }

/*
 *Funcion que escribe 4 arreglos TODOS de tamano n_x  en el tiempo t en un archivo .dat (dado su tamano).
 */
void escribirArreglos(double*t,double*x,double*y, double*z,int n_x,double energia, double alpha)
{
  FILE* archivo;
  //Sabemos que solo hay 4 arreglos
  double x0;
  double y0;
  double z0;
  x0=x[0];
  y0=y[0];
  char bufX[20];
  char bufY[20];
  char nmx= x0 -'0';
  char nmy= y0 -'0';
  int i;
  sprintf(bufX, "%f", energia);
  sprintf(bufY, "%f", alpha);
  char n1[50], n3[50], n2[50];
   strcpy(n1,  "trayectoria_");
   strcpy(n2, "_");
   strcpy(n3, ".dat");

   strcat(n1, bufX);
   strcat(n1, n2);
   strcat(n1, bufY);
   strcat(n1, n3);
   archivo = fopen(n1, "w");
  for(i=0;i<n_x;i++){
    fprintf(archivo, "%f \t %f \t %f \t %f \n", t[i], x[i], y[i], z[i]);
  }
  fclose(archivo);
}

/*
 * Dado un vector de n posiciones me calcula su magnitud
 */
double magnitudVector(double*vector,int n)
{
  int i;
  double suma;
  for(i=0;i<n;i++)
    {
      suma=suma+vector[i]*vector[i];
    }
  return (sqrt(suma));
}

/*
 *Metodo que me da el arreglo de las posiciones iniciales
 */
double*posicionesIniciales()
{
  double*posiciones;
  posiciones=crear_matriz(3,1);
  posiciones[0]=2.0*rT;
  posiciones[1]=0.0;
  posiciones[2]=0.0;
  return posiciones;
}

//Devuelve el vector de velocidadesInicialese la velocidad de acuerdo a la energia cinetica entrada por parametro
double*velocidadesIniciales(double energia, double alpha)
{
  double *velocidades;
  double magnitud;
  magnitud=magnitudVelocidadInicial(energia);
  velocidades=crear_matriz(3,1);
  velocidades[0]=0.0;
  velocidades[1]=magnitud*sin(alpha);
  velocidades[2]=magnitud*cos(alpha);
  return velocidades;
}

/*
 * Funcion que dada la energia cinetica me devuelve la magnitud de la velocidad
 */
double magnitudVelocidadInicial(double energia)
{
  return sqrt((c*c*(1.0-(1.0/((energia/(masa*c*c))+1.0)))));
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

