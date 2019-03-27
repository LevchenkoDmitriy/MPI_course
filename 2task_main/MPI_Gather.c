#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>

int main(int argc, char *argv[]){
	int size, rank;
	int buf = 0;//Пересылаемое сообщение
	double average_time, sum_time, start_time, end_time, error_time, array_sum, time = 0;
	MPI_Status Status; 
	int sbuf = 1;
	int recvbuf = 0;


	int N = 100000; //Количество запусков
	double* array = calloc(N, sizeof(double));//Массив для записи времени(нужно для вычисления среднего)

	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	//Измеряем время функции
	for(int i = 0; i < N; i++){
		MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

		//Собирает данные со всех и склеивает в 1 вектор.
		MPI_Gather(&sbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

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

	printf("[RANK %d] Av_Time  = %.10lf\n", rank, average_time);
	printf("[RANK %d] Wtick    = %.10lf\n", rank, MPI_Wtick());
	printf("[RANK %d] av_error = %.10lf\n", rank, error_time );

	MPI_Finalize();

	return 0;
}