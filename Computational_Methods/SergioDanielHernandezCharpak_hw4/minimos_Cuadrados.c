#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
float *descompocision_cholesky(float *matricita1, int n) ;
float *transpuesta(float *matricita1, int n, int m);
float *load_matrix(char *filename, int *n, int *m);
float *crear_matriz(int n, int m);
float *multiplicacion(float *matrix_A, float *matrix_B,  int n, int m, int p);
float *formarMatrizB(float *matrizDatos, float *matrizAT,int n, int m);
float *suma(float *matrix_A, float *matrix_B, int n, int m);
float *forwardSubstitution(float *matricitaLT,float *matricitaB, int n, int m, int n1, int m1);
float *backSubstitution(float *matricitaL, float*matricitaZ, int n, int m, int n1, int m1);
float *prepararRetorno(float*matrizX,int orden,float r);
float chiCuadrado(float *matrizDatos, float *matrizA,float*matricitaX, int n,int k, int m);
float coeficienteCorrelacion(float *matrizDatos, float *matrizA,float*matricitaX, int n,int k, int m);
float *formarMatrizEntradaA(float *matricitaDatosint ,int n_filas,int n_cols,int orden); 
float *formarMatrizY(float *matrizDatos, int n, int m);
/*
 *La funcion main retorna una matriz con la solucion (matrizX + coeficiente de determinacion r)
 */
float * main(int argc, char **argv){
  FILE *archivo;
  float *matricitaResultado;
  float *matricitaL;
  float* matrizModeloAB;
  float *matricitaB;
  float *matricitaLT;
  float *matricitaDatos;
  float *matricitaAT;
  float *matricitaA;
  float *matricitaBeforeCholesky;
  float *matricitaX;
  float *matricitaY;
  float *matricitaZ;
  float *matricitaARetornar;
  float r;
  int n_filas,n_cols ;
  int i, j, k,h;
  int orden = atoi(argv[2]); 
  //La matriz de Datos sera ahora o cos(2Theta) (modelo A) o Theta (modelo B) (columna1) y F (columna2) y tendra 1000 filas
  matricitaDatos = load_matrix(argv[1], &n_filas, &n_cols); //Cargamos la matriz con los datos
  //Ahora debemos formar la matriz pra Minimos Cuadrados para el modelo A
  matricitaA=formarMatrizEntradaA(matricitaDatos,n_filas,n_cols,orden); //Los Thetas estan en la primera columna de matricitaDatos
  matricitaAT = transpuesta(matricitaA,n_filas,orden+1); //Le sacamos la transpuesta
  matricitaBeforeCholesky=multiplicacion(matricitaAT,matricitaA,orden+1,n_filas,orden+1); //Las multiplicamos y obtenemos una matriz simetrica positiva
  matricitaL = descompocision_cholesky(matricitaBeforeCholesky,orden+1); //Realizamos la descompocision de cholesky y obtenemos la matriz triangular inferior L
  matricitaLT=transpuesta(matricitaL,orden+1,orden+1); //Obtenemos la matriz triangular superior LT, transpuesta de L
  matricitaY=formarMatrizY(matricitaDatos, n_filas, n_cols);
  matricitaB=multiplicacion(matricitaAT,matricitaY,orden+1,n_filas,1);//Formamos la Matriz B
  //Ya tenemos todas las matrices, ahora resolvamos esta @%(#$)%
  //Hacemos forwardSubstituttion en LZ=B para encontrar a Z
  //Hacemos back substitution para encontrar x con LTX=Y
  matricitaZ=forwardSubstitution(matricitaL,matricitaB,orden+1,orden+1,orden+1,1);
  matricitaX=backSubstitution(matricitaLT,matricitaZ,orden+1,orden+1,orden+1,1);
  //  r=chiCuadrado(matricitaY, matricitaA,matricitaX, n_filas, orden,n_cols);
  r=coeficienteCorrelacion(matricitaY, matricitaA,matricitaX, n_filas, orden,n_cols);
  //Codigo para imprimir los resultados
  matricitaARetornar=prepararRetorno(matricitaX,orden,r);
  // printf("Matriz A Retornar\n");
  for(k=0;k<orden+2;k++){
    for(h=0;h<1;h++){
      printf(" %f ", matricitaARetornar[k*1 + h]);
    }
    printf("\n");
  }
  return matricitaARetornar;
}

/*
 *Metodo que prepara la matriz a retornar, corrigiendo el tercer dato de x (que corresponde a 1/2g0) y agregando a la matriz una posicion con el coeficiente de determinacion r.
 */
float *prepararRetorno(float*matrizX,int orden, float r)
{
  int i;
  int n=orden+2;
  float* matrizR;
  matrizR = malloc(n* sizeof(float));
  for(i=0;i<n;i++)
  {
    if(i!=n-1)
      {
	matrizR[i]=matrizX[i];
      }
      else
      {
	matrizR[i]=r;
      }
  }
  return matrizR;
}

/*
 *Metodo que me da el coeficiente de correlacion del conjunto de datos con el modelo
 */
float coeficienteCorrelacion(float *matrizY, float * matrizA, float*matricitaX, int n, int orden, int m)
{
  int i,j,l,k,h;
  l=1;
  float r,rNominador,rDenominador,sumax,sumay,sumax2,sumay2,sumaxy;
  r=0;
  rNominador=0.0;
  rDenominador=0.0;
  //Ahora calculamos la matriz de Y teoricos
  float *matrizYteo=crear_matriz(n,l);
  matrizYteo=multiplicacion(matrizA,matricitaX,n,orden+1,l);
  sumax=0.0;
  sumax2=0.0;
  sumay=0.0;
  sumay2=0.0;
  sumaxy=0.0;
  for(i=0;i<n;i++)
  {
    for(j=0;j<l;j++)
    {
      sumax=sumax+matrizYteo[i*l+j];
      sumax2=sumax2+((matrizYteo[i*l+j])*(matrizYteo[i*l+j]));
      sumay=sumay+matrizY[i*l+j];
      sumay2=sumay2+((matrizY[i*l+j])*(matrizY[i*l+j]));
      sumaxy=sumaxy+((matrizYteo[i*l+j])*(matrizY[i*l+j]));
    }
  }
  rNominador=((n*sumaxy)-(sumax*sumay));
  rDenominador=(sqrt((n*sumax2)-(sumax*sumax)))*(sqrt((n*sumay2)-(sumay*sumay)));
  r=rNominador/rDenominador;
  return r;
}

/*
 *Metodo que me da el chicuadrado  del conjunto de datos con el modelo
 */
float chiCuadrado(float *matrizY, float * matrizA, float*matricitaX, int n, int orden, int m)
{
  int i,j,l,k,h;
  l=1;
  float r,sumaNominador,sumaDenominador,suma,sumaYTeo,sumaYTeo2,varYTeo;
  r=0;
  //Ahora calculamos la matriz de Y teoricos
  float *matrizYteo=crear_matriz(n,l);
  matrizYteo=multiplicacion(matrizA,matricitaX,n,orden+1,l);
  sumaYTeo=0.0;
  sumaYTeo2=0.0;
  sumaNominador=0.0;
  suma=0.0;
  sumaDenominador=0.0;
  for(i=0;i<n;i++)
  {
    for(j=0;j<l;j++)
    {
      sumaYTeo=sumaYTeo+(matrizYteo[i*l+j]);
      sumaYTeo2=sumaYTeo2+((matrizYteo[i*l+j])*(matrizYteo[i*l+j]));
      sumaNominador=(((matrizYteo[i*l+j])-(matrizY[i*l+j]))*((matrizYteo[i*l+j])-(matrizY[i*l+j])));
    }
  }
  varYTeo=(sumaYTeo2/n)-((sumaYTeo/n)*(sumaYTeo/n));
  r=sumaNominador/varYTeo;
  return r;
}


/*
 * Metodo que implementa la forward substitution dadas dos matrices y sus dimenciones
 */
float *forwardSubstitution(float *matricitaLT,float *matricitaB, int n, int m, int n1, int m1)
{
  int i,k,l;
  float suma;
  float* matrizZ=crear_matriz(n,m1);
  l=m1-1;
  for ( i = 0; i < n; i++)
  {
    suma=0.0;
    for ( k = 0; k < i; k++) 
      {
	suma=suma+(matrizZ[k*m1+l]*((matricitaLT[i*m+k])/(matricitaLT[i*m+i])));
      }
    matrizZ[i*m1+l]=((matricitaB[i*m1+l])/(matricitaLT[i*m+i]))-suma;
  }
  return matrizZ;
}
/*
 * Metodo que implementa backsubstitution dadas dos matrices y sus dimenciones
 */
float *backSubstitution(float *matricitaL, float*matricitaZ, int n, int m, int n1, int m1)
{
  float* matrizX=crear_matriz( n, m1);
  int i,k,l;
  float suma;
  l=m1-1;
  for ( i = n-1; i >=0; i--)
  {
    suma=0.0;
    for ( k = i+1; k < n; k++) 
      {
	suma=suma+(matrizX[k*m1+l]*((matricitaL[i*m+k])/(matricitaL[i*m+i])));
      }
    matrizX[i*m1+l]=((matricitaZ[i*m1+l])/(matricitaL[i*m+i]))-suma;
  }
  return matrizX;
}



/*
 *Metodo que forma la matriz para los minimos cuadrados dada una matriz, sus dimensiones y sus ordenes
 */
float *formarMatrizEntradaA(float *matrizDatos, int n, int m, int orden)
{
  int i,j,l;
  l=orden+1;
  float* matrizA=crear_matriz( n, l);
  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < l; j++) 
      {
	 matrizA[i*l + j]= pow((matrizDatos[i*m + 0]),j);
      }
  }
  return matrizA;
}

/*
 *Metodo que forma la matriz y
 */
float *formarMatrizY(float *matrizDatos, int n, int m)
{
  int i,j,l,k,h;
  l=1;
  float *matrizY=crear_matriz(n,l);
  //Debemos coger la segunda columna de MatrizDatos
  //Metemos en la matriz Y la segunda columna de matrizDatos
  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < l; j++) 
      {
	matrizY[i*l + j]= matrizDatos[i*m + l];
      }
  }
  return matrizY;
}


/*
 *Metodo que forma la matriz At*y
 */
float *formarMatrizB(float *matrizDatos, float *matrizAT,int n, int m)
{
  int i,j,l,k,h;
  l=m-1;
  k=m+1;
  float *matrizB=crear_matriz(k,l);
  float *matrizY=crear_matriz(n,l);
  //Debemos coger la segunda columna de MatrizDatos y Multiplicarla por matrizAT
  //Metemos en la matriz Y la segunda columna de matrizDatos
  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < l; j++) 
      {
	matrizY[i*l + j]= matrizDatos[i*m + l];
      }
  }
  //Ahora multiplicamos la matrizAT con la matrizY
  matrizB=multiplicacion(matrizAT,matrizY,k,n,l);
  return matrizB;
}

/*
 * Metodo que me da la transupuesta de la matriz  pasada como parametro
*/
float *transpuesta(float *matricita1, int n, int m)
{
  float* matrizT=  crear_matriz( m, n);
   int i,j;
    for ( i = 0; i < m; i++)
      {
        for ( j = 0; j < n; j++) 
	  {
	    matrizT[i * n + j]= (matricita1[j * m + i]) ;
	  }
       }
    return matrizT;
}


/*
 * Metodo que descompone una matriz definida positiva y real en la multiplicacion de una matriz triangular baja L
 * Usando la descompocision de Cholsky
 * http://en.wikipedia.org/wiki/Cholesky_decomposition
*/

float *descompocision_cholesky(float *matricita1, int n) 
{
   float* matrizL=  crear_matriz( n, n);
   double suma = 0;
   int i,j,k;
    for (i=0;i<n;i++)
      {
        for (j= 0;j<(i+1);j++) 
	  {
            for (k=0;k<j;k++)
	      {
                suma=suma+((matrizL[i*n+k])*(matrizL[j*n+k]));
	      }
	    if(i==j)
	      {
		matrizL[i*n+j]= sqrt(matricita1[i*n+i]-suma) ;
	      }
	    else
	      {
		matrizL[i*n+j]=((1.0/matrizL[j*n+j])*(matricita1[i*n+j]-suma));
	      }    
	    suma=0;
	  }
       }
    return matrizL;
}

/*
 * Metodo que multiplica dos matriz y pone el resultado en la tercera.
 * @params matricita1, matricita2 las dos matrices a multiplicar
 * @params matricitaResultado en donde se imprime el resultado
 * @params n numero de filas primera matriz, m numero de columnas primera matriz y numero filas segunda matriz, p numero de columnas segunda matriz
 */
float *multiplicacion(float *matricita1, float *matricita2,  int n, int m,int p)
{
  float *matricitaResultado=crear_matriz(n,p);
  int i,j,k;
  float suma;
  for (i=0;i<n;i++)
    {
    for(j=0;j<p;j++)
      {
	suma = 0.0;
      for(k=0;k<m;k++)
	{
	suma=suma+(matricita1[i*m+k]*matricita2[k*p+j]);
	}
      matricitaResultado[i*p+j]=suma;
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
float *suma(float *matricita1, float *matricita2,  int n, int m)
{
  float *matricitaResultado=crear_matriz( n, m);
  int i,j,k;
  for (i=0;i<n;i++)
    {
    for(j=0;j<m;j++)
      {
	matricitaResultado[i*m +j ]=0.0; // Ponemos toda la matriz en ceros para evitar problemas
	matricitaResultado[i*m +j ]=matricita1[i*m + j]+matricita2[i*m +j ];
      }
    }
  return matricitaResultado;
}


float *load_matrix(char *filename, int *n, int *m){
  float *matrix;
  FILE *in;
  int n_row, n_cols;
  int i;
  int j;
  if(!(in=fopen(filename, "r"))){
    printf("Problem opening file %s\n", filename);
    exit(1);
  }
  fscanf(in, "%d %d\n", &n_row, &n_cols);
  matrix = malloc(n_row * n_cols * sizeof(float));
  for(i=0;i<n_row;i++){
    for(j=0;j<n_cols;j++){
      fscanf(in, "%f", &matrix[i*n_cols + j]);
    }
  }    
  *n = n_row;
  *m = n_cols;
  return matrix;
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

