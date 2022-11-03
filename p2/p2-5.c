//P2 Paralelismo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define FOREVER 1

//MPI_Bcast propia
int MPI_BinomialColectiva(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
	int i, source, dest, rank, numprocs;
	int devuelto;														
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	//i es el paso (tal como explica la diapo 4 del tema 2)
	for(i=1;;i++){
		//Tal vez me toque enviar
		if(rank < pow(2,i-1)){
			dest = rank + pow(2,i-1);
			if(dest >= numprocs){
				break;
			}
			if((devuelto=MPI_Send(buffer, count, datatype, dest, 0, comm)) != MPI_SUCCESS){
				return devuelto;
			}
		}
		// Tal vez me toque recibir
		else {
			source = rank-pow(2,i-1);
			// Efectivamente me toca recibir
			if(source < pow(2,i-1)){
				if((devuelto=MPI_Recv(buffer, count, datatype, source, 0, comm, MPI_STATUS_IGNORE)) != MPI_SUCCESS){
					return devuelto;
				}
			}
		}
	}
	
	return MPI_SUCCESS; 
}

//MPI_Reduce propia
int MPI_FlattreeColectiva(void *buff, void *recvbuff, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm){
	
	int i, myrank, rank, numprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	
	//Si soy el root recibo e integro (sumo) mis datos y los datos del resto de los procesos																									
	if(myrank == root){
		int devuelto;
		//Sumo (inicializo copiando) mis datos (buff) al buffer destino (recvbuff)
		for(i=0;i<count;i++){
			((int*)recvbuff)[i] = ((int*)buff)[i];
		}
		//Recibo del resto y los integro en recvbuff
		int tmp[count];
		for(rank=0; rank<numprocs; rank++){
			if(rank != myrank){
				if((devuelto=MPI_Recv(&tmp,count,datatype,rank,0,comm,MPI_STATUS_IGNORE)!=MPI_SUCCESS)){
					return devuelto;
				}
				for(i=0;i<count;i++){
					((int*)recvbuff)[i] = ((int*)recvbuff)[i] + tmp[i];

				}
			}
		}
		return MPI_SUCCESS;
	//Si no soy el root, envío al root los datos
	} else {
		return MPI_Send(buff, count, datatype, root, 0, comm);
	}
}

int main(int argc, char *argv[]){
	
    int i, n, generados, count, subCount;
    double PI25DT = 3.141592653589793238462643;
    double pi,x,y,z;
   
	//Variables e inicialización del entorno MPI
	int numprocs, rank;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	srand(time(NULL)+rank);

	//Bucle eterno (sin condición de salida). Finalizará (con break) si
	//se introduce 0 como número de intervalos
	while(FOREVER){
		count = 0;
		subCount = 0;
		
		if(rank == 0){
			//Leo desde teclado el numero de puntos.
			printf("Enter the number of points: (0 quits) \n");
			scanf("%d",&n);
		}
		
		//Difusión desde el proceso 0 al resto de los procesos de N
		MPI_BinomialColectiva(&n,1,MPI_INT,0,MPI_COMM_WORLD);
			
		//Si el valor introducido para n es 0 --> FIN
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
		
		//Integro mi resultado
		count = subCount;
		
		//Cada proceso (incluído el 0) entrega su subcuenta y la cuenta
		//total se recoge en la variable cuenta 
		MPI_FlattreeColectiva(&subCount, &count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		
		pi = ((double) count/(double) n)*4.0;
			
		if(rank == 0){
			//Mostramos el PI calculado y el error respecto al proporcionado
			printf("pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
		}
	
	}
	
	MPI_Finalize();
	   
    return 0;
}
