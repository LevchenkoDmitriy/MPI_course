#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

void openNumbers(char num[2][512], char* str);

int main(int argc, char *argv[]){
	MPI_Status Status; 
	int size, rank;
	int block = 0;
	char num[2][512] = {0};
	
	char flag_transfer_between_blocks = 0;

	MPI_Init(&argc, &argv);

	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	block = 512/size;

	if(rank == 0){
		openNumbers(num, "numbers.txt");
	}

	char* recv_buf_1 = (char*)calloc(block, sizeof(char));
	char* recv_buf_2 = (char*)calloc(block, sizeof(char));

	char* gather_send_buf = (char*)calloc(block, sizeof(char));
	char* gather_recv_buf = (char*)calloc(512, sizeof(char));

	MPI_Scatter(num[0], block, MPI_CHAR, recv_buf_1, block, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Scatter(num[1], block, MPI_CHAR, recv_buf_2, block, MPI_CHAR, 0, MPI_COMM_WORLD);

    char bit_transfer_in_block = 0, bit_transfer_in_block_speculative = 0;
 	char* send_buf_1 = (char*)calloc(block, sizeof(char));
 	char* send_buf_2 = (char*)calloc(block, sizeof(char));

 	//Складывем числа в каждом блоке 
    for (int i = block - 1; i >= 0; i--) {
        send_buf_1[i] = recv_buf_1[i] + recv_buf_2[i] + bit_transfer_in_block;
        bit_transfer_in_block = 0;
        if (((unsigned char)send_buf_1[i]) >= 100) {
            bit_transfer_in_block = 1;
            send_buf_1[i] = send_buf_1[i] - 100;
        }
    }

 
    if(rank != size - 1) {
    	//В самый правый процесс ничего не придет, для остальных считаем дополнительно с битом переноса
    	bit_transfer_in_block_speculative = 1;

    	for (int i = block - 1; i >= 0; i--) {
        	send_buf_2[i] = recv_buf_1[i] + recv_buf_2[i] + bit_transfer_in_block_speculative;
        	bit_transfer_in_block_speculative = 0;

        	if (((unsigned char)send_buf_2[i]) >= 100) {
            	bit_transfer_in_block_speculative = 1;
            	send_buf_2[i] = send_buf_2[i] - 100;

        	}
    	}
    	//Ожидаем бит переноса от предыдущего блока. Самый правый сразу отправит
    	MPI_Recv(&flag_transfer_between_blocks, 1, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &Status);

    	//Если от предыдущего блока нет бита переноса, отправляем первый буфер
    	if(flag_transfer_between_blocks == 0){

    		if (rank != 0) 
    			MPI_Bsend(&bit_transfer_in_block, 1, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD);

    		for (int i = 0; i < block; i++){
 				gather_send_buf[i] = send_buf_1[i];
    		}

    	}else {
        	if (rank != 0) 
        		MPI_Bsend(&bit_transfer_in_block_speculative, 1, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD);
       		for (int i = 0; i < block; i++){
            	gather_send_buf[i] = send_buf_2[i];
        	}
    	}
   	}else {
        MPI_Bsend(&bit_transfer_in_block, 1, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD);
        for(int i = 0; i < block; i++){
            gather_send_buf[i] = send_buf_1[i];
        }
    }

    MPI_Gather(gather_send_buf, block, MPI_CHAR, gather_recv_buf, block, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0){
        if(flag_transfer_between_blocks == 1) 
        	printf("%d", 1);
        for (int j = 0; j < 512; j++) {
            printf("%d",(unsigned char)gather_recv_buf[j]);
        }
            printf("\n");
	}

	MPI_Finalize();
	return 0;

}



void openNumbers(char num[2][512], char* str){
	FILE *numbers;
	int counter = 0, count = 0, i = 0, j = 0;//counter считает количество цифр в числе(<1024)
									  //count считает запись в вспомогательные переменные(<2)
									  //j считает числа(<2)
								 	  //i нужна записи в массив
	char sym[2] = {0};

	numbers = fopen (str, "r");

	if(numbers == NULL){
		printf("No such file in directory\n");
	}

	//Записываем 2 числа в 2 переменные, "склеиваем" и записываем в 1 из массивов
	int c;
	while ((c = fgetc(numbers)) != EOF) {
    	if(counter != 1024){
    		if(count != 1){
    			//Считываем первый символ
    			sym[0] = c;
    			count++;
    		}else{
    			//Считываем второй символ, склеиваем, записываем
    			sym[1] = c;
    			sym[0] = sym[0] - '0';
    			sym[1] = sym[1] - '0';
    			num[j][i] = sym[0]*10 + sym[1];
    			count = 0;
    			i++;

    		}
    	}else{
    		//если дошли до конца первого числа
    		j = 1;
    		counter = 0;
    		i = 0;
    	}
    	counter++;
    	if((counter == 1024)&&(j == 1)){
    		sym[0] = c;
    		c = fgetc(numbers);
    		sym[1] = c;
    		sym[0] = sym[0] - '0';
    		sym[1] = sym[1] - '0';
    		num[j][i] = sym[0]*10 + sym[1];
    		break;
    	}
	}

	fclose(numbers);
}

