/*The Mandelbrot set is a fractal that is defined as the set of points c
in the complex plane for which the sequence z_{n+1} = z_n^2 + c
with z_0 = 0 does not tend to infinity.*/

/*This code computes an image of the Mandelbrot set.*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define DEBUG 1

#define          X_RESN  1024  /* x resolution */
#define          Y_RESN  1024  /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

typedef struct complextype
{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end)
{
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 +
         (t_end.tv_sec - t_ini.tv_sec);
}

int main (int argc, char** argv)
{
  
  int numprocs, rango, flops=0, f_totales=0;

  MPI_Init(&argc, &argv);
        
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rango);
  
  /* Mandelbrot variables */
  Compl   z, c;
  float   lengthsq, temp;
  int *vres, *res[Y_RESN];
  int i, j, k, f_procesos[numprocs];   //Operaciones en punto flotante por proceso
  
  /* Timestamp variables */
  struct timeval  ti, tf;

  /* Allocate result matrix of Y_RESN x X_RESN */
  //Lo reserva si es root
  if(rango == 0){
    vres = (int *) malloc(Y_RESN * X_RESN * sizeof(int));
    if (!vres)
    {
      fprintf(stderr, "Error allocating memory\n");
      return 1;
    }
    for (i=0; i<Y_RESN; i++)
      res[i] = vres + i*X_RESN;
  }
  
  //Resultado de cada proceso
  int resultados_proceso[Y_RESN/numprocs+1][X_RESN];
  
  //División del trabajo
  int filas_totales_proceso[numprocs];
  int fila_comienzo[numprocs];
  int elementos_totales_proceso[numprocs];
  int desplazamiento[numprocs];
  
  for(i=0; i<numprocs; i++){
    if( i < Y_RESN%numprocs ){
      filas_totales_proceso[i] = Y_RESN/numprocs + 1;
      fila_comienzo[i] = i*filas_totales_proceso[i];
    }else{
      filas_totales_proceso[i] = Y_RESN/numprocs;
      fila_comienzo[i] = fila_comienzo[i-1] + filas_totales_proceso[i-1];
    }
    elementos_totales_proceso[i] = filas_totales_proceso[i] * X_RESN;
    desplazamiento[i] = fila_comienzo[i] * X_RESN;		
  }
  
  /* Start measuring time */
  gettimeofday(&ti, NULL);
  
  /* Calculate and draw points */
  //Cada proceso calcula su parte  
  for(i=fila_comienzo[rango]; i<fila_comienzo[rango]+filas_totales_proceso[rango]; i++)
  {
    for(j=0; j < X_RESN; j++)
    {
      z.real = z.imag = 0.0;
      c.real = X_MIN + j * (X_MAX - X_MIN)/X_RESN;
      c.imag = Y_MAX - i * (Y_MAX - Y_MIN)/Y_RESN;
      flops += 4;
      k = 0;

      do
      {    /* iterate for pixel color */
        temp = z.real*z.real - z.imag*z.imag + c.real;
        z.imag = 2.0*z.real*z.imag + c.imag;
        z.real = temp;
        lengthsq = z.real*z.real+z.imag*z.imag;
        flops += 10;
        k++;
      } while (lengthsq < 4.0 && k < maxIterations);

      if (k >= maxIterations) resultados_proceso[i-fila_comienzo[rango]][j] = 0;
      else resultados_proceso[i-fila_comienzo[rango]][j] = k;
    }
  }

  /* End measuring time */
  gettimeofday(&tf, NULL);
  fprintf (stderr, "Procesador: %d (PERF) Computation time (seconds) = %lf\n",rango, get_seconds(ti,tf));
  
  //Esperamos por todos los procesos para medir el tiempo
  MPI_Barrier(MPI_COMM_WORLD);
  
  gettimeofday(&ti, NULL);
  
  //Metemos en vres los resultados de cada proceso
  MPI_Gatherv(resultados_proceso, elementos_totales_proceso[rango], MPI_INT, vres, elementos_totales_proceso, desplazamiento, MPI_INT, 0, MPI_COMM_WORLD);
  
  gettimeofday(&tf, NULL);
  fprintf (stderr, "Procesador: %d (PERF) Communication time (seconds) = %lf\n",rango, get_seconds(ti,tf));
  
  //Metemos en flops el número de operaciones de punto flotante de cada proceso
  MPI_Gather(&flops, 1, MPI_INT, f_procesos, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  //Root suma todos los flops y calcula el trabajo
  if(rango == 0){
    for(i = 0; i < numprocs; i++){
      f_totales += f_procesos[i];
    }
  //Si es menor que 1 es que trabajó de más y si es mayor trabajó de menos
    for(i = 0; i < numprocs; i++){
      fprintf(stderr, "Procesador %d: %.6f  flops = %d\n", i, f_totales*1.0 / (f_procesos[i]*1.0*numprocs), f_procesos[i]);
    }
  }


  /* Print result out */
  if(DEBUG && (rango == 0)) {
    for(i=0;i<Y_RESN;i++) {
      for(j=0;j<X_RESN;j++)
              printf("%3d ", res[i][j]);
      printf("\n");
    }
  }

  MPI_Finalize();
}
