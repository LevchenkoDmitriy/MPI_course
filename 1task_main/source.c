#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]){
	int size, rank; //Размер текущего MPI_COMM_WORLD и ранг текущего процесса в MPI_COMM_WORLD
	int buf; //Буфер для приема сообщений 
	MPI_Status Status; //Информация о результате приёма сообщения

	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	if(rank == 0){
		printf("Hello from %d\n", rank);
		
		for(int i = 1; i < size; i++){
			MPI_Recv(&buf, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			printf("Hello from %d\n", buf);
		}

	}else{
		MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}