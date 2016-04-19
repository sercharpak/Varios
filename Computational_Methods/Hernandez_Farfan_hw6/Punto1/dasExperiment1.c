//Para el punto 1 Sergio
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define pi 3.14159265359

float*crearArregloCero(int n_points);
float*crearArregloEquiEspaciado(float x_ini,float x_fin,int n_points);
float*copiarArreglo(float* original, int n_points);
float derivadaX(float x, float y);
float derivadaY(float x, float y);
float* pasoRungerKutta4(float paso, float x, float y);
void rungerKuttaPredatorPrey(float* x, float* y, float* t, int n_points);
void escribirArreglosPoblaciones(float*t,float*x,float*y, int n_x);

int  main(int argc, char **argv)
{
  float x0 = atof(argv[1]); 
  float y0 = atof(argv[2]); 
  /* 
  float xeq,yeq;
  float A=20.0;
  float B=1.0;
  float C=30.0;
  float D=1.0;
  xeq=C/D;
  yeq=A/B;
  */
  float ti=0.0;
  float tf=1.0;
  int n_points=100;
  float* t;
  float* x;
  float* y;
  t=crearArregloEquiEspaciado(ti,tf,n_points);
  x=crearArregloCero(n_points);
  x[0]=x0;
  y=crearArregloCero(n_points);
  y[0]=y0;
  //RungerKutta Time!
  rungerKuttaPredatorPrey(x,y,t,n_points);
  return 0;
}

    /*
     *Funcion que dado un paso (en el tiempo), un x y un y (al tiempo t) calcula el x y y al tiempo t+paso usando Runger Kutta 4to orden
     */
    float* pasoRungerKutta4(float paso,float x, float y)
    {
      //x(t+paso) = x+(k1+2*k2+2*k3+k4)/6;
      //y(t+paso) = y+(m1+2*m2+2*m3+m4)/6;
      float k1,k2,k3,k4;
      float m1,m2,m3,m4;
      k1=paso*derivadaX(x,y);
      m1=paso*derivadaY(x,y);
      k2=paso*derivadaX(x+0.5*k1,y+0.5*m1);
      m2=paso*derivadaY(x+0.5*k1,y+0.5*m1);
      k3=paso*derivadaX(x+0.5*k2,y+0.5*m2);
      m3=paso*derivadaY(x+0.5*k2,y+0.5*m2);
      k4=paso*derivadaX(x+k3,y+m3);
      m4=paso*derivadaY(x+k3,y+m3);
      float* retorno;
      retorno = malloc(sizeof(float) * 2);
      retorno[0]=x+(k1+2.0*k2+2.0*k3+k4)/6.0;
      retorno[1]=y+(m1+2.0*m2+2.0*m3+m4)/6.0;
      return retorno;
    }

  /*
   *Funcion que Ejecuta RungerKutta para el sistema Predador Presa
   */
  void rungerKuttaPredatorPrey(float* x, float* y, float* t, int n_points)
  {
    int i;
    float x0, y0, paso;
     paso=t[1]-t[0];
    for(i=0;i<n_points;i++)
      {
	float* temp;
	temp=pasoRungerKutta4(paso, x[i], y[i]);
	x[i+1]=temp[0];
	y[i+1]=temp[1];
      }
    escribirArreglosPoblaciones(t,x,y, n_points);
  }

  /*
   *Funcion que especifica la derivada de X en funcion del tiempo con A y B iguales a 20 y 1 respectivamente
   */
float derivadaX(float x, float y)
{
    int A = 20;
    int B = 1;
    return A*x-B*y*x;
}

  /*
   *Funcion que especifica la derivada de Y en funcion del tiempo con C y D iguales a 30  y 1 respectivamente
   */
float derivadaY(float x, float y)
{
  int C = 30;
  int D = 1;
  return -C*y+D*x*y;
}

/*
 *Funcion que escribe 3 arreglos TODOS de tamano n_x  en el tiempo t en un archivo .dat (dado su tamano).
 */
void escribirArreglosPoblaciones(float*t,float*x,float*y, int n_x)
{
  FILE* archivo;
  //Sabemos que solo hay 3 arreglos
  float x0;
  float y0;
  x0=x[0];
  y0=y[0];
  char bufX[20];
  char bufY[20];
  char nmx= x0 -'0';
  char nmy= y0 -'0';
  int i;
  sprintf(bufX, "%f", x0);
  sprintf(bufY, "%f", y0);
  char n1[50], n3[50], n2[50];
   strcpy(n1,  "poblaciones_");
   strcpy(n2, "_");
   strcpy(n3, ".dat");

   strcat(n1, bufX);
   strcat(n1, n2);
   strcat(n1, bufY);
   strcat(n1, n3);
   archivo = fopen(n1, "w");
  for(i=0;i<n_x;i++){
    fprintf(archivo, "%f \t %f \t %f \n", t[i], x[i], y[i]);
  }
  fclose(archivo);
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
