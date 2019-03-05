#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[]){
	int size, rank;
	int buf;
	MPI_Status Status; 

	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	if(rank == 0) {
        printf("Hello from %d\n", rank);
        MPI_Send(&rank, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&buf, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(rank == buf + 1) {
            printf("Hello from %d\n", rank);
            buf++;
            if(rank != size - 1) {
                MPI_Send(&buf, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            }
        }
    }

	MPI_Finalize();
	return 0;
}