#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

void myBcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
void myScatterInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm comm);
void myReduceAdd(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
void myGatherInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm comm);

int main(int argc, char *argv[]){
	int size, rank;
	MPI_Status Status; 
	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/*Тест myBcast
	int arg;
	if(rank == 0) 
		arg = 5;
	 else 
	 	arg = 0;

	printf("[RANK %d] %d\n", rank, arg);
	myBcast(&arg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

	printf("[RANK %d] %d\n", rank, arg);
	*/

	/*
	//Тест Reduce
	int buf[4] = {0, 5, 2, 8};
	int recvbuf[4] = {0};
	
		myReduceAdd(&buf, &recvbuf, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);
	
		printf("[RANK %d] %d, %d, %d, %d\n", rank, recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);
	*/

	/*
	//Тест Scatter
	int buf[8] = {1, 2, 3, 4, 5, 6, 7, 8};
	int recv[2] = {0};

	myScatterInt(&buf, 2, MPI_INTEGER, &recv, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

	printf("[RANK %d] %d, %d\n", rank, recv[0], recv[1]);
	*/
	/*
	//Тест Gather
	int buf[2] = {0};
	if(rank == 0) {buf[0] = 1; buf[1] = 2;}
	if(rank == 1) {buf[0] = 3; buf[1] = 4;}
	if(rank == 2) {buf[0] = 5; buf[1] = 6;}
	if(rank == 3) {buf[0] = 7; buf[1] = 8;}

	int recv[8] = {0};

	myGatherInt(&buf, 2, MPI_INTEGER, &recv, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

	for(int i = 0; i < 8; i++){
		printf("[Rank %d], %d = %d\n", rank, i, recv[i]);
	}
*/

	MPI_Finalize();

	return 0;		
}

void myBcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
	int size, rank;
	MPI_Comm_size (comm, &size);
	MPI_Comm_rank (comm, &rank);
	MPI_Status Status;

	MPI_Barrier(comm);

	if(rank == root){
		for(int i = 0; i < size; i++)
			MPI_Send(buffer, count, datatype, i, 0, comm);
	}else
		MPI_Recv(buffer, count, datatype, root, MPI_ANY_TAG, comm, &Status);

	MPI_Barrier(comm);
}

void myReduceAdd(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
	int size, rank;
	MPI_Comm_size (comm, &size);
	MPI_Comm_rank (comm, &rank);
	MPI_Status Status;

	MPI_Barrier(comm);

	int** a = calloc(count, sizeof(int));
	for(int i = 0; i < count; i++)
		a[i] = calloc(size, sizeof(int));
	

	if(rank == root){
		for(int i = 0; i < size; i++){
			for(int j = 0; j < count; j++){
				if(i != root)
					MPI_Recv(&a[i][j], 1, MPI_INTEGER, i, MPI_ANY_TAG, comm, &Status);
				else
					a[i][j] = *(int*)(sendbuf + j*sizeof(int));

			}
		}
	}else{
		for(int i = 0; i < count; i++)
			MPI_Send(((void*)(int*)(sendbuf + i*sizeof(int))), 1, MPI_INTEGER, root, 0, comm);
	}

	int* b = calloc(count, sizeof(int));
	if(rank == root){
		for(int i = 0; i < size; i++){
			for(int j = 0; j < count; j++){
				b[i] += a[j][i];
				*((int*)recvbuf + i) = b[i];
			}
		}
	}
	
	MPI_Barrier(comm);
}

void myScatterInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm comm){
	int size, rank;
	MPI_Comm_size (comm, &size);
	MPI_Comm_rank (comm, &rank);
	MPI_Status Status;

	MPI_Barrier(comm);

	if(rank == root){
		for(int i = 0; i < size; i++)
			MPI_Send(send_data + i * send_count * sizeof(int), send_count, send_datatype, i, 0, comm);
	}

	MPI_Recv(recv_data, recv_count, recv_datatype, root, MPI_ANY_TAG, comm, &Status);

	MPI_Barrier(comm);
}

void myGatherInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm comm){
	int size, rank;
	MPI_Comm_size (comm, &size);
	MPI_Comm_rank (comm, &rank);
	MPI_Status Status;

	MPI_Barrier(comm);
	
	MPI_Send(send_data, send_count, send_datatype, root, 0, comm);

	if(rank == root)
		for(int i = 0; i < size; i++)
			MPI_Recv(recv_data + i * recv_count * sizeof(int), recv_count, recv_datatype, i, MPI_ANY_TAG, comm, &Status);

	MPI_Barrier(comm);
}
