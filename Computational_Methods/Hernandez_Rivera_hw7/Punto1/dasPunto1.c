#include <stdio.h>
#include <stdlib.h>
#include <math.h>
float *CondicionInicial(int lengx);
int main(int argc, char **argv)
{
  //Inicializacion y Constantes
  int i;
  int j;
  int lx=101;
  float *U_initial;
  FILE *in;
  float tiempo=120;
  
  U_initial=CondicionInicial(lx);
  float delta_x=1;
  float delta_t=1;
  float nt=tiempo/delta_t;
  float T=40;
  float rho=atof(argv[1]);
  char filename[100];
  snprintf(filename, sizeof(char)*100,"string_%f.dat",(float) rho);
  float c=sqrt(T/rho);
  float r=c*delta_t/delta_x;
  float *U_future;
  float *U_present;
  float *U_past;
  U_future=malloc(sizeof(float)*lx);
  U_present=malloc(sizeof(float)*lx);
  U_past=malloc(sizeof(float)*lx);
  printf("%f\n",r);
//Abrir archivo para escribir
  in =fopen(filename,"w");
  if(!in)
    {
      printf("Problems opening the file \n");
    }
  // Escribir inicial
  for(i=0;i<lx;i++)
    {
        fprintf(in,"%f\t",U_initial[i]);
    }
   fprintf(in,"\n");
  //Primera Iteracion
  U_future[0]=0;
  U_future[lx-1]=0;
  for(i=1;i<lx-1;i++)
    {
      U_future[i]=U_initial[i]+(pow(r,2.0)/2.0)*(U_initial[i+1]-2.0*U_initial[i]+U_initial[i-1]);
    }
  //Copiar vectores
  for(i=0;i<lx;i++)
    {
      U_past[i]=U_initial[i];
      U_present[i]=U_future[i];
        fprintf(in,"%f\t",U_present[i]);
    }
   fprintf(in,"\n");

  // Resto de iteraciones
   for(j=1;j<nt;j++)
    {
      for(i=1;i<lx-1;i++)
	{
	  U_future[i]=2*(1-pow(r,2))*U_present[i]-U_past[i]+pow(r,2)*(U_present[i+1]+U_present[i-1]);
	}
      //Copiar vectores de nuevo
      for(i=0;i<lx;i++)
	{
	 
	  U_past[i]=U_present[i];
	  U_present[i]=U_future[i];
	  // if(i%10==0)
	  // {
	            fprintf(in,"%f\t",U_present[i]);
		    //}
	}
        fprintf(in,"\n");
    }
  


  
  //Escribir
  for(i=0;i<lx;i++)
    {
      // fprintf(in,"%f\n",U_present[i]);
    }
  //float test=U_initial[80]+(pow(r,2.0)/2.0)*(U_initial[80+1]-2.0*U_initial[80]+U_initial[80-1]);
  //printf("%f\n",tiempo/0.01);
  return 0;
}
float *CondicionInicial(int lengx)
{
  int i;
  float *ui;
  ui=malloc(sizeof(float)*lengx);
  for(i=0;i<lengx;i++)
      {
      float ii=i+0.0;
      //Aqui va la condicion inicial
      
      if(ii<=0.8*100)
	{
	  ui[i]=1.25*ii/100;
	}
      else
	{
	  ui[i]=5-5*ii/100;
	  }
      // ui[i]=exp(-pow((ii-50)/10,2)); 
	}
      
  return ui;
}
  
