#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>

// Preprocess the command line, read in matrix A, B from input files, allocate
// memory for buffers, i.e, fstreama, fstreamb to cache them. Suppose A, B's 
// size are n1*n2, n2*n3, then n1~n3 will be stored at dim[0~2]
// Return value 0 means no error occurred during preprocessing, otherwise a
// non-zero returns.
int setup(int argc, char** argv, char* &fstreama, char* &fstreamb, int* dim){
	int error = 0;
	if(argc < 4){
		printf("Invalid arguments!\n");
		printf("Usage: ./Cannon filea fileb filec\n");
		printf("filea, fileb and filec are file names for matrix A, B and C\n");
		return 1;
	}

	FILE *fha, *fhb;
	int fsizea, fsizeb;

	if(!(fha = fopen(argv[1], "r"))){
		printf("Can't open matrix file %s, Errno=%d\n", argv[1], errno);
		return 1;
	}

	if(!(fhb = fopen(argv[2], "r"))){
		printf("Can't open matrix file %s, Errno=%d\n", argv[2], errno);
		return 1;
	}

	struct stat fstata, fstatb;
	stat(argv[1], &fstata);
	stat(argv[2], &fstatb);
	fsizea = fstata.st_size;
	fsizeb = fstatb.st_size;

	fstreama = (char *)malloc(fsizea);
	fstreamb = (char *)malloc(fsizeb);
	fread(fstreama, sizeof(char), fsizea, fha);
	fread(fstreamb, sizeof(char), fsizeb, fhb);

	int n1, n2, n3, n4;
	n1 = ((int*)fstreama)[0];
	n2 = ((int*)fstreama)[1];
	n3 = ((int*)fstreamb)[0];
	n4 = ((int*)fstreamb)[1];

	if(n1 <= 0 || n2 <= 0 || n3 <= 0 || n4 <= 0 || n2 != n3){
		printf("Matrix size error, %d*%d with %dx%d\n", n1, n2, n3, n4);
		return 1;
	}

	if(fsizea < (sizeof(int)*2 + sizeof(double)*n1*n2)){
		printf("Actual size of A mismatches with stated size\n");
		return 1;
	}

	if(fsizeb < (sizeof(int)*2 + sizeof(double)*n3*n4)){
		printf("Actual size of B mismatches with stated size\n");
		return 1;
	}

	dim[0] = n1;
	dim[1] = n2;
	dim[2] = n4;

	fclose(fha);
	fclose(fhb);
	return 0;
}

void scatter_matrix(double *fstream, int n1, int n2, double *buf, int rootp, int tag){
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int maxrows = (n1 + rootp - 1) / rootp;
	int maxcols = (n2 + rootp - 1) / rootp;
	if(myrank == 0){
		int pad_size = sizeof(int)*128;
		int buf_size = sizeof(int)*2 + sizeof(double)*maxrows*maxcols;
		double *tempbuf;
		tempbuf=(double *)malloc(buf_size+pad_size);
		if(!tempbuf){
			// fail("Memory allocation failed\n");
			printf("Memory allocation failed\n");
		}
		int i, j;
		// memset(buf, 0, buf_size);
		memset((int *)buf, 0, sizeof(int)*2);
		memset((double *)(buf + sizeof(int)*2), 0, buf_size-sizeof(int)*2);
		for(i = 0; i < rootp; i++){
			for(j = 0; j < rootp; j++){
				int p = 0, q = 0;
				int imin = i*maxrows;
				int jmin = j*maxcols;
				// memset(tempbuf, 0, buf_size);
				memset((int *)tempbuf, 0, sizeof(int)*2);
				memset((double *)(tempbuf + sizeof(int)*2), 0, buf_size-sizeof(int)*2);
				((int *)tempbuf)[0] = (i == rootp - 1) ? (n1 - imin) : maxrows;
				((int *)tempbuf)[1] = (j == rootp - 1) ? (n2 - jmin) : maxcols;
				// printf("(%d, %d): row = %d, col = %d\n", i, j, ((int *)tempbuf)[0], ((int *)tempbuf)[1]);
				for(p = 0; p < ((int *)tempbuf)[0]; p++, imin++){
					for(q = 0; q < ((int *)tempbuf)[1]; q++){
						((double *)(tempbuf + sizeof(int)*2))[p*maxcols+q] = fstream[imin*n2+jmin+q];
						// printf("%f ", ((double *)(tempbuf + sizeof(int)*2))[p*maxcols+q]);
					}
					// printf("\n");
				}
				if(i*rootp + j == 0){
					// memcpy(buf, tempbuf, buf_size);
					memcpy((int *)buf, (int *)tempbuf, sizeof(int)*2);
					memcpy((double *)(buf + sizeof(int)*2), (double *)(tempbuf + sizeof(int)*2), buf_size-sizeof(int)*2);
				}
				else{
					MPI_Send((int *)tempbuf, 2, MPI_INT, i*rootp+j, tag, MPI_COMM_WORLD);
					MPI_Send((double *)(tempbuf+sizeof(int)*2), maxrows*maxcols, MPI_DOUBLE, i*rootp+j, tag, MPI_COMM_WORLD);
				}
			}
		}
	}
	else{
		MPI_Recv((int *)buf, 2, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv((double *)(buf+sizeof(int)*2), maxrows*maxcols, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	}
	/*
	for(int k = 0; k < rootp*rootp; k++){
		if(myrank == k){
			printf("core%d, ((int*)buf)[0] = %d, ((int*)buf)[1] = %d\n", k, ((int*)buf)[0], ((int*)buf)[1]);
			for(int i = 0; i < ((int*)buf)[0]; i++){
				for(int j = 0; j < ((int*)buf)[1]; j++){
					printf("%f ", ((double *)(buf+ sizeof(int)*2))[i*maxcols+j]);
				}
				printf("\n");
			}
		}
	}
	*/
}

void gather_matrix(double *fstream, int n1, int n3, double *buf, int rootp){
	MPI_Status status;
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int maxrows = (n1 + rootp - 1) / rootp;
	int maxcols = (n3 + rootp - 1) / rootp;
	if(myrank == 0){
		int buf_size = sizeof(int)*2 + sizeof(double)*maxrows*maxcols;
		int i, j;
		for(i = 0; i < rootp; i++){
			for(j = 0; j < rootp; j++){
				int p = 0, q = 0;
				int imin = i*maxrows;
				int jmin = j*maxcols;
				if(i*rootp + j != 0){
					memset(buf, 0, buf_size);
					// MPI_Recv((int *)buf, 2, MPI_INT, i*rootp+j, 102, MPI_COMM_WORLD, &status);
					MPI_Recv((double *)(buf+sizeof(int)*2), maxrows*maxcols, MPI_DOUBLE, i*rootp+j, 102, MPI_COMM_WORLD, &status);
				}
				((int *)buf)[0] = (i == rootp - 1) ? (n1 - imin) : maxrows;
				((int *)buf)[1] = (j == rootp - 1) ? (n3 - jmin) : maxcols;
				// printf("(%d, %d): row = %d, col = %d\n", i, j, ((int *)buf)[0], ((int *)buf)[1]);
				for(p = 0; p < ((int *)buf)[0]; p++, imin++){
					for(q = 0; q < ((int *)buf)[1]; q++){
						fstream[imin*n3+jmin+q] = ((double *)(buf + sizeof(int)*2))[p*maxcols+q];
						// printf("%f ", ((double *)(buf + sizeof(int)*2))[p*maxcols+q]);
					}
					// printf("\n");
				}
			}
		}
	}
	else{
		// MPI_Send((int *)buf, 2, MPI_INT, 0, 102, MPI_COMM_WORLD);
		MPI_Send((double *)(buf+sizeof(int)*2), maxrows*maxcols, MPI_DOUBLE, 0, 102, MPI_COMM_WORLD);
	}
}

// Compute C = A * B. A is a n1*n2 matrix. B is a n2*n3 matrix.
void matmul(double* A, double* B, double* C, int n1, int n2, int n3, int maxcols_a, int maxcols_b){
	#define A(i,j)  *(A + i*maxcols_a + j)
	#define B(i,j)  *(B + i*maxcols_b + j)
	#define C(i,j)  *(C + i*maxcols_b + j)

	for(int i = 0; i < n1; i++){
		for(int j = 0; j < n3; j++){
			// C(i,j) = 0.0;
			for(int k = 0; k < n2; k++){
				C(i,j) += A(i,k)*B(k,j);
			}
		}
	} 
}

// get the index of the next processor
int next_proc(int row, int col, int rootp){
	int temp = ((row + rootp) % rootp) * rootp + (col + rootp) % rootp;
	return temp;
}

void shuffle(double *A, double *bufA, int bufA_size, double *B, double *bufB, int bufB_size, 
	int n1, int n2, int n3, int rootp, int myrank){
	int i, j;
	MPI_Status status;
	int maxrows_a = (n1 + rootp - 1) / rootp;
	int maxcols_a = (n2 + rootp - 1) / rootp;
	int maxcols_b = (n3 + rootp - 1) / rootp;
	int cur_row = myrank / rootp;
	int cur_col = myrank % rootp;
	// for matrix block in row i, shift left for i times
	for(i = 0; i < cur_row; i++){
		// send matrix in A to the next proc and receive matrix from the last proc to bufA
		MPI_Sendrecv((int *)A, 2, MPI_INT, next_proc(cur_row, cur_col-1, rootp), 103, 
			(int *)bufA, 2, MPI_INT, next_proc(cur_row, cur_col+1, rootp), 103, MPI_COMM_WORLD, &status);
		MPI_Sendrecv((double *)(A+sizeof(int)*2), maxrows_a*maxcols_a, MPI_DOUBLE, next_proc(cur_row, cur_col-1, rootp), 106, 
			(double *)(bufA+sizeof(int)*2), maxrows_a*maxcols_a, MPI_DOUBLE, next_proc(cur_row, cur_col+1, rootp), 106, MPI_COMM_WORLD, &status);
		// move matrix in bufA to A
		// memcpy(A, bufA, bufA_size);
		memcpy((int *)A, (int *)bufA, sizeof(int)*2);
		memcpy((double *)(A + sizeof(int)*2), (double *)(bufA + sizeof(int)*2), bufA_size-sizeof(int)*2);
		// clear bufA
		memset(bufA, 0, bufA_size);
	}
	/*
	for(int k = 0; k < rootp*rootp; k++){
		if(myrank == k){
			printf("shuffle A core %d, ((int*)A)[0] = %d, ((int*)A)[1] = %d\n", k, ((int*)A)[0], ((int*)A)[1]);
			for(i = 0; i < ((int*)A)[0]; i++){
				for(j = 0; j < ((int*)A)[1]; j++){
					printf("%f ", ((double *)(A + sizeof(int)*2))[i*maxcols_a+j]);
				}
				printf("\n");
			}
		}
	}
	*/
	// for matrix block in col i, shift up for j times
	for(j = 0; j < cur_col; j++){
		// send matrix in B to the next proc and receive matrix from the last proc to bufB
		MPI_Sendrecv((int *)B, 2, MPI_INT, next_proc(cur_row-1, cur_col, rootp), 104, 
			(int *)bufB, 2, MPI_INT, next_proc(cur_row+1, cur_col, rootp), 104, MPI_COMM_WORLD, &status);
		MPI_Sendrecv((double *)(B+sizeof(int)*2), maxcols_a*maxcols_b, MPI_DOUBLE, next_proc(cur_row-1, cur_col, rootp), 107, 
			(double *)(bufB+sizeof(int)*2), maxcols_a*maxcols_b, MPI_DOUBLE, next_proc(cur_row+1, cur_col, rootp), 107, MPI_COMM_WORLD, &status);
		// move matrix in bufB to B
		// memcpy(B, bufB, bufB_size);
		memcpy((int *)B, (int *)bufB, sizeof(int)*2);
		memcpy((double *)(B + sizeof(int)*2), (double *)(bufB + sizeof(int)*2), bufB_size-sizeof(int)*2);
		// clear bufB
		memset(bufB, 0, bufB_size);
	}
	/*
	for(int k = 0; k < rootp*rootp; k++){
		if(myrank == k){
			printf("shuffle B core %d, ((int*)B)[0] = %d, ((int*)B)[1] = %d\n", k, ((int*)A)[0], ((int*)A)[1]);
			for(i = 0; i < ((int*)B)[0]; i++){
				for(j = 0; j < ((int*)B)[1]; j++){
					printf("%f ", ((double *)(B + sizeof(int)*2))[i*maxcols_b+j]);
				}
				printf("\n");
			}
		}
	}
	*/
}

void cannon(double *A, double *bufA, int bufA_size, double *B, double *bufB, int bufB_size, 
	double *C, int bufC_size, int n1, int n2, int n3, int rootp, int myrank){
	MPI_Status status;
	int i, j;
	int maxrows_a = (n1 + rootp - 1) / rootp;
	int maxcols_a = (n2 + rootp - 1) / rootp;
	int maxcols_b = (n3 + rootp - 1) / rootp;
	memset(C, 0, bufC_size);
	int cur_row = myrank / rootp;
	int cur_col = myrank % rootp;
	// calculate block of matrix C in every proc
	for(i = 0; i < rootp; i++){
		// matrix multiply for current part
		matmul((double*)(A+sizeof(int)*2), (double*)(B+sizeof(int)*2), (double*)(C+sizeof(int)*2), 
			((int*)A)[0], ((int*)A)[1], ((int*)B)[1], maxcols_a, maxcols_b);
		// matrix block of A shift left for one time and get the block from the last proc
		MPI_Sendrecv((int *)A, 2, MPI_INT, next_proc(cur_row, cur_col-1, rootp), 100, 
			(int *)bufA, 2, MPI_INT, next_proc(cur_row, cur_col+1, rootp), 100, MPI_COMM_WORLD, &status);
		MPI_Sendrecv((double *)(A+sizeof(int)*2), maxrows_a*maxcols_a, MPI_DOUBLE, next_proc(cur_row, cur_col-1, rootp), 108, 
			(double *)(bufA+sizeof(int)*2), maxrows_a*maxcols_a, MPI_DOUBLE, next_proc(cur_row, cur_col+1, rootp), 108, MPI_COMM_WORLD, &status);
		// matrix block of B shift up for one time and get the block from the last proc
		MPI_Sendrecv((int *)B, 2, MPI_INT, next_proc(cur_row-1, cur_col, rootp), 101, 
			(int *)bufB, 2, MPI_INT, next_proc(cur_row+1, cur_col, rootp), 101, MPI_COMM_WORLD, &status);
		MPI_Sendrecv((double *)(B+sizeof(int)*2), maxcols_a*maxcols_b, MPI_DOUBLE, next_proc(cur_row-1, cur_col, rootp), 109, 
			(double *)(bufB+sizeof(int)*2), maxcols_a*maxcols_b, MPI_DOUBLE, next_proc(cur_row+1, cur_col, rootp), 109, MPI_COMM_WORLD, &status);
		// move bufA to A
		// memcpy(A, bufA, bufA_size);
		memcpy((int *)A, (int *)bufA, sizeof(int)*2);
		memcpy((double *)(A + sizeof(int)*2), (double *)(bufA + sizeof(int)*2), bufA_size-sizeof(int)*2);
		// move bufB to B
		// memcpy(B, bufB, bufB_size);
		memcpy((int *)B, (int *)bufB, sizeof(int)*2);
		memcpy((double *)(B + sizeof(int)*2), (double *)(bufB + sizeof(int)*2), bufB_size-sizeof(int)*2);
	}
}

int main(int argc, char** argv){
	// Suppose A: n1*n2, B: n2*n3. n1~n3 are read from input files
	int n1, n2, n3;
	// Buffers for matrix A, B, C. Because A, B will be shifted, so they each have two buffers.
	double *A, *B, *C, *bufA, *bufB;
	// On proc 0, buffers to cache matrix files of A, B and C
	char *fstreama, *fstreamb, *fstreamc;
	// rank of current proc, total number of procs and root of p
	int myrank, p, rootp;
	// MPI status
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	rootp = sqrt(p);
	if(p != rootp*rootp){
		// fail("Processor number must be a square!\n");
		printf("Processor number must be a square!\n");
	}

	// On proc 0, preprocess the command line, read in files for A, B and put their sizes in dim[]
	int dim[3];
	if(myrank == 0){
		if(setup(argc, argv, fstreama, fstreamb, dim)){
			MPI_Finalize(); // Something error during preprocessing
			exit(-1);
		}
	}

	MPI_Bcast(dim, 3, MPI_INT, 0, MPI_COMM_WORLD);
	n1 = dim[0];
	n2 = dim[1];
	n3 = dim[2];

	// Allocate memories for A, B, C, bufA and bufB.
	// Suppose an m*n matrix is 2D block-distributed on a rootp*rootp processor grid.
	// If rootp doesn't divide m or n, then submatrixes won't have the same size.
	// Because we will shift A, B, we allocate memories according to the max.
	// rows and cols of A and B.
	int maxrows_a = (n1 + rootp - 1) / rootp;
	int maxcols_a = (n2 + rootp - 1) / rootp;
	int maxrows_b = maxcols_a;
	int maxcols_b = (n3 + rootp - 1) / rootp;
	int bufA_size = sizeof(int)*2 + sizeof(double)*maxrows_a*maxcols_a;
	int bufB_size = sizeof(int)*2 + sizeof(double)*maxrows_b*maxcols_b;
	int bufC_size = sizeof(int)*2 + sizeof(double)*maxrows_a*maxcols_b;
	int pad_size = sizeof(int)*128;

	char *buf;
	buf=(char *)malloc(bufA_size*2 + bufB_size*2 + bufC_size + pad_size*5);
	if(!buf){
		// fail("Memory allocation failed\n");
		printf("Memory allocation failed\n");
	}
	A = (double*)buf;
	bufA = (double*)(buf + bufA_size + pad_size);
	B = (double*)(buf + bufA_size*2 + pad_size*2);
	bufB = (double*)(buf + bufA_size*2 + bufB_size + pad_size*3);
	C = (double*)(buf + bufA_size*2 + bufB_size*2 + pad_size*4);

	if(myrank == 0){
		printf("setup finish!\n");
	}

	// Proc 0 scatters A, B to other procs in a 2D block distribution fashion
	/*
	if(myrank == 0){
		printf("scatter_matrix A, n1 = %d, n2 = %d\n", n1, n2);
	}
	*/
	scatter_matrix((double*)(fstreama + sizeof(int)*2), n1, n2, A, rootp, 99);
	MPI_Barrier(MPI_COMM_WORLD);
	/*
	if(myrank == 0){
		printf("scatter_matrix B, n2 = %d, n3 = %d\n", n2, n3);
	}
	*/
	scatter_matrix((double*)(fstreamb + sizeof(int)*2), n2, n3, B, rootp, 105);
	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = MPI_Wtime();

	// pre set blocks of matrix A and B in every proc
	shuffle(A, bufA, bufA_size, B, bufB, bufB_size, n1, n2, n3, rootp, myrank);
	MPI_Barrier(MPI_COMM_WORLD);
	// Compute C = A * B by Cannon Algorithm
	cannon(A, bufA, bufA_size, B, bufB, bufB_size, C, bufC_size, n1, n2, n3, rootp, myrank);

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime() - elapsed_time;
	/*
	if(myrank == 0){
		printf("cannon finish!\n");
	}
	*/
	// Proc 0 gathers C from other procs and write it out
	FILE *fhc;
	int fsizec = sizeof(int)*2 + sizeof(double)*n1*n3;
	if(myrank == 0){
		if(!(fhc = fopen(argv[3], "w"))){
			printf("Can't open file %s, Errno=%d\n", argv[3], errno);
			MPI_Finalize();
		}
		fstreamc = (char *)malloc(fsizec);
		((int *)fstreamc)[0] = n1;
		((int *)fstreamc)[1] = n3;
		printf("open filec\n");
	}
	gather_matrix((double *)(fstreamc + sizeof(int)*2), n1, n3, C, rootp);
	MPI_Barrier(MPI_COMM_WORLD); // Make sure proc 0 read all it needs
	/*
	if(myrank == 0){
		printf("gather_matrix finish!\n");
	}
	*/
	if(myrank == 0){
		printf("Cannon Algorithm: multiply a %d*%d with a %d*%d, use %.2f(s)\n", n1, n2, n2, n3, elapsed_time);
		fwrite(fstreamc, sizeof(char), fsizec, fhc);
		fclose(fhc);
		free(fstreama);
		free(fstreamb);
		free(fstreamc);
	}
	if(buf){
		free(buf);
	}
	MPI_Finalize();
	return 0;
}