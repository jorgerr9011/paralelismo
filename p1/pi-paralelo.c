#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

/*int calculos (int rank, int numprocs, int subCount, int n, double x, double y, double z){
	int i;
	//double pi,x,y,z;

	for (i = rank+1; i <= n; i=i+numprocs) {
				// Get the random numbers between 0 and 1
				x = ((double) rand()) / ((double) RAND_MAX);
				y = ((double) rand()) / ((double) RAND_MAX);

				// Calculate the square root of the squares
				z = sqrt((x*x)+(y*y));

				// Check whether z is within the circle
				if(z <= 1.0)
						subCount++;
			}
}*/

int main(int argc, char *argv[]){
	
    int i, n, done = 0, count, subCount;
    double PI25DT = 3.141592653589793238462643;
    double pi,x,y,z;
   
	//Variables e inicialización del entorno MPI
	int numprocs, rank;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Si no se inicializa no genera aleatorios bien
	srand(time(NULL)+rank);

	while(!done){
		count = 0;
		subCount = 0;
		if(rank == 0){
			printf("Enter the number of points: (0 quits) \n");
			scanf("%d",&n);

			//Primero envío n a los otros procesos 
			for(i=1; i<numprocs; i++){
				MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
			
			if(n==0){
				break;
			}	
			
			for (i = rank+1; i <= n; i=i+numprocs) {
				// Get the random numbers between 0 and 1
				x = ((double) rand()) / ((double) RAND_MAX);
				y = ((double) rand()) / ((double) RAND_MAX);

				// Calculate the square root of the squares
				z = sqrt((x*x)+(y*y));

				// Check whether z is within the circle
				if(z <= 1.0)
						subCount++;
			}
			//subCount = calculos(rank, numprocs, subCount, n, x, y, z);
			
			//Esto sería el Send de root
			count = subCount;
			
			//Recibo los resultados de los otros procesos
			for(i=1; i<numprocs; i++){
				MPI_Recv(&subCount, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				count = count + subCount;
			}
			pi = ((double) count/(double) n)*4.0;
			
			printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
			
		}else{
			//Espero la recepción de n (el número de puntos a generar)
			MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			
			if(n==0){
				break;
			}
					
			for (i = rank+1; i <= n; i=i+numprocs) {
				// Get the random numbers between 0 and 1
				x = ((double) rand()) / ((double) RAND_MAX);
				y = ((double) rand()) / ((double) RAND_MAX);

				// Calculate the square root of the squares
				z = sqrt((x*x)+(y*y));

				// Check whether z is within the circle
				if(z <= 1.0)
						subCount++;
			}
			//subCount = calculos(rank, numprocs, subCount, n, x, y, z);
	    	
			//Devuelvo el subresultado
			MPI_Send(&subCount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);	
		}
	
	}
	
	MPI_Finalize();
	   
    return 0;
}
