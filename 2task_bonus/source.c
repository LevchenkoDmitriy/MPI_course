#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <string.h>


void myBcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

void myScatterInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, 
	MPI_Datatype recv_datatype, int root, MPI_Comm comm);

void myReduceAdd(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

void myGatherInt(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, 
	MPI_Datatype recv_datatype, int root, MPI_Comm comm);


int main(int argc, char *argv[]){
	int size, rank;
	int N = 100000;
	MPI_Status Status; 
	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

//---------------------Тест myBcast----------------------//
	//Передача значения
	int arg;
	if(rank == 0) 
		arg = 5;
	 else 
	 	arg = 0;

	double* array = calloc(N, sizeof(double));
	double average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		myBcast(&arg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[31m[myBcast]\x1b[0m averageTime    = %.10lf\n", average_time);
		printf("\x1b[31m[myBcast]\x1b[0m MPI_Wtick      = %.10lf\n", MPI_Wtick());
		printf("\x1b[31m[myBcast]\x1b[0m averageError   = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
	/*
	//Тест работы myBcast
	printf("[RANK %d] %d\n", rank, arg);
	myBcast(&arg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	printf("[RANK %d] %d\n", rank, arg);
	*/

//--------------------------------------------//
//---------------------Тест Bcast----------------------//
	//Передача значения
	if(rank == 0) 
		arg = 5;
	 else 
	 	arg = 0;

	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		MPI_Bcast(&arg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[32m[Bcast]\x1b[0m averageTime      = %.10lf\n", average_time);
		printf("\x1b[32m[Bcast]\x1b[0m MPI_Wtick        = %.10lf\n", MPI_Wtick());
		printf("\x1b[32m[Bcast]\x1b[0m averageError     = %.10lf\n", error_time );
	}
	memset(array, 0, N*sizeof(array));

//--------------------------------------------//
//---------------------Тест myReduce----------------------//
	//Тест Reduce
	int buf[4] = {0, 5, 2, 8};
	int recvbuf[4] = {0};

	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		myReduceAdd(&buf, &recvbuf, 4, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[33m[myReduce]\x1b[0m averageTime   = %.10lf\n", average_time);
		printf("\x1b[33m[myReduce]\x1b[0m MPI_Wtick     = %.10lf\n", MPI_Wtick());
		printf("\x1b[33m[myReduce]\x1b[0m averageError  = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
	
	//printf("[RANK %d] %d, %d, %d, %d\n", rank, recvbuf[0], recvbuf[1], recvbuf[2], recvbuf[3]);
//--------------------------------------------//		
//---------------------Тест myReduce----------------------//
	//Тест Reduce
	memset(recvbuf, 0, 4 * sizeof(int));

	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		MPI_Reduce(&buf, &recvbuf, 4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[34m[Reduce]\x1b[0m averageTime     = %.10lf\n", average_time);
		printf("\x1b[34m[Reduce]\x1b[0m MPI_Wtick       = %.10lf\n", MPI_Wtick());
		printf("\x1b[34m[Reduce]\x1b[0m averageError    = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
//--------------------------------------------//		
//---------------------Тест myScatter----------------------//
	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	int buf1[8] = {1, 2, 3, 4, 5, 6, 7, 8};
	int recv1[2] = {0};

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		myScatterInt(&buf1, 2, MPI_INTEGER, &recv1, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[35m[myScatter]\x1b[0m averageTime  = %.10lf\n", average_time);
		printf("\x1b[35m[myScatter]\x1b[0m MPI_Wtick    = %.10lf\n", MPI_Wtick());
		printf("\x1b[35m[myScatter]\x1b[0m averageError = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
	//printf("[RANK %d] %d, %d\n", rank, recv[0], recv[1]);
//--------------------------------------------//	
//---------------------Тест Scatter----------------------//
	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;

	memset(recv1, 0, 2 * sizeof(int));

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		MPI_Scatter(&buf1, 2, MPI_INTEGER, &recv1, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[36m[Scatter]\x1b[0m averageTime    = %.10lf\n", average_time);
		printf("\x1b[36m[Scatter]\x1b[0m MPI_Wtick      = %.10lf\n", MPI_Wtick());
		printf("\x1b[36m[Scatter]\x1b[0m averageError   = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
//--------------------------------------------//	
//---------------------Тест myGather----------------------//
	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;
	int buf2[2] = {0};
	if(rank == 0) {buf2[0] = 1; buf2[1] = 2;}
	if(rank == 1) {buf2[0] = 3; buf2[1] = 4;}
	if(rank == 2) {buf2[0] = 5; buf2[1] = 6;}
	if(rank == 3) {buf2[0] = 7; buf2[1] = 8;}

	int recv2[8] = {0};

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		myGatherInt(&buf2, 2, MPI_INTEGER, &recv2, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[32m[myGather]\x1b[0m averageTime   = %.10lf\n", average_time);
		printf("\x1b[32m[myGather]\x1b[0m MPI_Wtick     = %.10lf\n", MPI_Wtick());
		printf("\x1b[32m[myGather]\x1b[0m averageError  = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
	/*
	for(int i = 0; i < 8; i++){
		printf("[Rank %d], %d = %d\n", rank, i, recv[i]);
	}*/
//--------------------------------------------//	
//---------------------Тест Gather----------------------//
	average_time = 0, sum_time = 0, start_time = 0, end_time = 0, error_time = 0, array_sum = 0, time = 0;
	memset(recv2, 0, 8*sizeof(int));

	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		MPI_Gather(&buf2, 2, MPI_INTEGER, &recv2, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

		array[i] = end_time - start_time;

		time = end_time - start_time;
		sum_time += time;
	}

	//Подсчёт среднеквадратичной ошибки
	average_time = sum_time/N;

	for(int i = 0; i < N; i++){
		array[i] = (array[i] - average_time)*(array[i] - average_time);
		array_sum += array[i];
	}

	error_time = sqrt(array_sum/(N*(N-1)));

	if(rank == 0){
		printf("\x1b[33m[Gather]\x1b[0m averageTime     = %.10lf\n", average_time);
		printf("\x1b[33m[Gather]\x1b[0m MPI_Wtick       = %.10lf\n", MPI_Wtick());
		printf("\x1b[33m[Gather]\x1b[0m averageError    = %.10lf\n", error_time );
	}

	memset(array, 0, N*sizeof(array));
//--------------------------------------------//	

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

	int** a = calloc(size, sizeof(int*));
	for(int i = 0; i < size; i++)
		a[i] = calloc(count, sizeof(int));

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