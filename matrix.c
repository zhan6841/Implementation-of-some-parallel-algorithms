#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <pthread.h>

int main(int argc, char *argv[]){
	int M, N, P;
	int myid, numprocs;
	int i,j;
	int numsend, sender;
	MPI_Status status;
	int numthreads;
	pthread_t *tids;
	float *A_row, *C_row;
	struct threadArg *targs;

	M = atoi(argv[1]);
	N = atoi(argv[2]);
	P = atoi(argv[3]);

	struct threadArg{
		int tid;
		float (*B)[P];
		float *A_row;
		float *C_row;
		int numthreads;
	};

	void* worker(void* arg){
		int i, j;
		struct threadArg* myarg = (struct threadArg*)arg;
		for(i = myarg->tid; i < P; i += myarg->numthreads){
			myarg->C_row[i+i*60] = 0.0;
			// calculate C[i][j]
			for(j = 0; j < N; j++){
				myarg->C_row[i+i*60] += myarg->A_row[j] * myarg->B[j][i];
			}
		}
		return NULL;
	}

	float A[M][N], B[N][P], C[M][P];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	// initialize matrix A and B in process 0
	if(!myid){
		printf("M = %d, N = %d, P = %d\n", M, N, P);
		printf("matrix A:\n");
		for(i = 0; i < M; i++){
			for(j = 0; j < N; j++){
				A[i][j] = i * j + 1;
				printf("%6.2f ", A[i][j]);
			}
			printf("\n");
		}
		printf("matrix B:\n");
		for(i = 0; i < N; i++){
			for(j = 0; j < P; j++){
				B[i][j] = i * j + 1;
				printf("%6.2f ", B[i][j]);
			}
			printf("\n");
		}
		// printf("matrix initialize\n");
	}

	// broadcast matrix B to all processes
	MPI_Bcast(B[0], N * P, MPI_FLOAT, 0, MPI_COMM_WORLD);

	// process 0: allocate assignments and receive results
	if(!myid){
		// allocate one line of matrix A to each processes
		j = (numprocs - 1) < M ? (numprocs - 1) : M;
		for(i = 1; i <= j; i++){
			MPI_Send(A[i-1], N, MPI_FLOAT, i, 99, MPI_COMM_WORLD);
		}
		numsend = j;

		for(i = 1; i <= M; i++){
			sender = (i - 1) % (numprocs - 1) + 1;
			MPI_Recv(C[i-1], P, MPI_FLOAT, sender, 100, MPI_COMM_WORLD, &status);
			// if there are lines left in matrix A that was not allocated
			// allocate one line to the process that send the result
			if(numsend < M){
				MPI_Send(A[numsend], N, MPI_FLOAT, sender, 99, MPI_COMM_WORLD);
				numsend++;
			}
			else{
				MPI_Send(&j, 0, MPI_INT, sender, 0, MPI_COMM_WORLD);
			}
		}

		// printf("assignments allocate\n");

		// output the result
		printf("matrix C:\n");
		for(i = 0; i < M; i++){
			for(j = 0; j < P; j++){
				printf("%6.2f ", C[i][j]);
			}
			printf("\n");
		}
	}
	// other processes: receive the assignments from process 0, calculate
	// send the result back to process 0
	else{
		// get the number of CPU in the node where the process is running
		numthreads = get_nprocs();
		// printf("in process %d numthreads = %d\n", myid, numthreads);
		// numthreads = (numthreads < P) ? numthreads : P;
		// the array tids is used to save the thread ID
		tids = (pthread_t *)malloc(numthreads * sizeof(pthread_t));
		// the line from matrix A
		A_row = (float *)malloc(N * sizeof(float));
		// one line of the result matrix C
		C_row = (float *)malloc((P + 60) * sizeof(float));
		// parameters of threads
		// since the there are many paras, just send result
		targs = (struct threadArg *)malloc(numthreads * sizeof(struct threadArg));
		// printf("malloc targs\n");
		for(i = 0; i < numthreads; i++){
			targs[i].tid = i;
			targs[i].B = B;
			targs[i].A_row = A_row;
			targs[i].C_row = C_row;
			targs[i].numthreads = numthreads;
		}
		// printf("targs\n");
		while(1){
			// receive one line of matrix A from process 0
			MPI_Recv(A_row, N, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			// printf("recv\n");
			// if tag == 0, break
			if(status.MPI_TAG == 0){
				break;
			}
			// create threads to calculate
			for(i = 0; i < numthreads; i++){
				pthread_create(&tids[i], NULL, worker, &targs[i]);
				// printf("create thread %d\n", i);
			}
			// wait until all the threads finish
			for(i = 0; i < numthreads; i++){
				pthread_join(tids[i], NULL);
				// printf("thread %d join\n", i);
			}
			for(i = 0; i < P; i++){
				C_row[i] = C_row[i+i*60];
			}
			// send the results to process 0
			MPI_Send(C_row, P, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
}
