//-------------------------SE USA GATHERV----------------------------------------
//FUNCIONA AUNQUE Y_RESN NO SEA MULTIPLO DE NUMPROCS
//
//Distingue entre tiempos de cómputo y tiempos de comunicaciones
//
//Se mide el balanceo de carga mediante el nº de flops consumidos
//----------------------------------------------------------------------

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
#define          Y_RESN  1024 /* y resolution */

/* Boundaries of the mandelbrot set */
#define           X_MIN  -2.0
#define           X_MAX   2.0
#define           Y_MIN  -2.0
#define           Y_MAX   2.0

/* More iterations -> more detailed image & higher computational cost */
#define   maxIterations  1000

typedef struct complextype{
  float real, imag;
} Compl;

static inline double get_seconds(struct timeval t_ini, struct timeval t_end){
  return (t_end.tv_usec - t_ini.tv_usec) / 1E6 +
         (t_end.tv_sec - t_ini.tv_sec);
}

int get_k(int i, int j, int *flops, float pasoX, float pasoY){
	
	Compl z, c;
	float lengthsq, temp;
	int k;
	
	z.real = z.imag = 0.0;
	c.real = X_MIN + j * pasoX;
	c.imag = Y_MAX - i * pasoY;
	//Lo anterior ejecuta 4 flops
	*flops += 4;
	k = 0;

	//Iteramos la ecuanción de Mandelbrot mientras el cuadrado
	//del módulo no supere 4 y k no supere el nº máximo de
	//iteraciones permitido
	do { /* iterate for pixel color (falso color)*/
		temp = z.real*z.real - z.imag*z.imag + c.real;
		z.imag = 2.0*z.real*z.imag + c.imag;
		z.real = temp;
		lengthsq = z.real*z.real+z.imag*z.imag;
		k++;
		//Lo anterior ejecuta 10 flops
		*flops += 10;
	} while (lengthsq < 4.0 && k < maxIterations);
	
	return k;
}


int main (int argc, char** argv){
	
	int numprocs, rank, flops=0, flops_total=0;
	
	//Calculamos (antes de entrar en los procesos paralelos) los
	//pasos (incrementos) de X e Y en cada cambio de fila y columna
	float pasoX = (X_MAX - X_MIN)/X_RESN;
	float pasoY = (Y_MAX - Y_MIN)/Y_RESN;

	//Se inician los procesos paralelos
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			
	/* Mandelbrot variables */
	//Estos cálculos no se pueden hacer antes de iniciar el código
	//paralelizado, ya que no sabemos el número de procesos que se arrancan
	int i, j, k;
	int filas_por_proceso = Y_RESN/numprocs;
	
	//Matriz donde cada proceso guardará su resultados parciales
	//de las filas que le correspondan
	int res_parcial[filas_por_proceso+1][X_RESN];
	
	//Matriz que almacenará el nº de iteraciones realizado en cada punto
	int *vres, *res[Y_RESN];
	
	//Vector que almacena el nº de flops realizadas por cada proceso
	int flops_por_proceso[numprocs];
	
	//Vectores que usará GatherV
	//Los procesos de rank < Y_RESN % numprocs harán Y_RESN/numprocs + 1,
	//y los restantes harán Y_RESN/numprocs
	int filas[numprocs];
	int filaInicial[numprocs];
	int count[numprocs];
	int displ[numprocs];
	
	//Llenado de los arrays
	for(i=0; i<numprocs; i++){
		if( i < Y_RESN%numprocs ){
			filas[i] = Y_RESN/numprocs + 1;
			filaInicial[i] = i*filas[i];
		}else{
			filas[i] = Y_RESN/numprocs;
			filaInicial[i] = filaInicial[i-1] + filas[i-1];
		}
		count[i] = filas[i] * X_RESN;
		displ[i] = filaInicial[i] * X_RESN;		
	}
	//Variables para la medición de tiempos.
	struct timeval ti, tf;
	int microseconds;
	
	//Si soy rank==0
	if(rank==0){
		//Reservo espacio para la matriz resultado de tamaño Y_RESN x X_RESN
		vres = (int *) malloc( Y_RESN * X_RESN * sizeof(int));
		if (!vres){
			fprintf(stderr, "Error allocating memory\n");
			return 1;
		}
		for(i=0; i<Y_RESN; i++){
			res[i] = vres + i*X_RESN;
		}
	}

	/* Start measuring time */
	gettimeofday(&ti, NULL);

	/* Calculate and draw points */
	//Cada proceso calcula la fracción de filas que le corresponde
	for(i=filaInicial[rank]; i<filaInicial[rank]+filas[rank]; i++){
		for(j=0; j < X_RESN; j++){
			k = get_k(i, j, &flops, pasoX, pasoY);
			//Si k alcanzó (sin superar el módulo es valor 2) el máximo
			//nº de iteraciones permitidas colocamos 0, si el módulo superó 2
			//antes de alcanzar maxIterations colocamos k
			if (k >= maxIterations) res_parcial[i-filaInicial[rank]][j] = 0;
			else res_parcial[i-filaInicial[rank]][j] = k;
		}
	}

	
	//Fin de la computación. Se calcula el tiempo de cálculo consumido y se muestra
	gettimeofday(&tf, NULL);
	fprintf (stderr, "Rank:%d (PERF) Computation Time (seconds) = %lf\n", rank, get_seconds(ti,tf));
	
	//Esperamos a que todos los procesos terminen sus cálculos para que las medidas de tiempos de 
	//comunicación no se vean afectadas por los tiempos de cálculo
	MPI_Barrier(MPI_COMM_WORLD);
	
	//Inicio de las comunicaciones
	gettimeofday(&ti, NULL);

	//Se integra en res (matriz con la info de todos los puntos analizados) los resultados
	//parciales de cada uno de los procesos
	MPI_Gatherv(res_parcial, count[rank], MPI_INT, vres, count, displ, MPI_INT, 0, MPI_COMM_WORLD);
	
	//Fin de las comunicaciones. Se calcula el tiempo de comunicaciones consumido y se muestra
	gettimeofday(&tf, NULL);
	fprintf (stderr, "Rank:%d (PERF) Communication Time (seconds) = %lf\n", rank, get_seconds(ti,tf));
	
	//El proceso root recupera el nº de floats de cada proceso.
	MPI_Gather(&flops, 1, MPI_INT, flops_por_proceso, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//El proceso 0 totaliza el nº de flops consumidos por todos los procesos
	if(rank == 0){
		for(i = 0; i < numprocs; i++){
			flops_total += flops_por_proceso[i];
		}
		//Se muestra una métrica del reparto del trabajo realizado del tipo:
		//"¿He trabajado más(<1) o menos(>1) de lo que debería haber trabajado
		//si la carga de trabajo (en flops) estuviera bien repartida?"
		for(i = 0; i < numprocs; i++){
			fprintf(stderr, "b%d = %5f  flops = %d\n", i, (flops_total*1.0) / (flops_por_proceso[i]*1.0*numprocs), flops_por_proceso[i]);
		}
	}

	//Se saca por la salida estándar (stdout) los valores de k almacenados en la matriz
	//tal cual (posición) se corresponden a un barrido de pantalla.
	//Una separación (espacio en blanco) entre elementos de la misma fila y un cambio
	//de línea entre líneas. Este es un formato que aceptará la librería de Python
	//para mostrar la matriz de valor es de k (nº de iteraciones) como un mapa de color
	if( DEBUG && (rank == 0)) {
		for(i=0;i<Y_RESN;i++) {
		  for(j=0;j<X_RESN;j++)
				  printf("%3d ", res[i][j]);
		  printf("\n");
		}
	}
	
	MPI_Finalize();
}
