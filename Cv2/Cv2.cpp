#include "pch.h"

int main(int argc, char **argv) {
	/*
	int nprocs, myrank;
	int isend[5], irecv[15];
	int sendcount, displs[5] = {3, 10, 0, 1, 6};
	int recvcounts[5] = { 3, 5, 1, 2, 4 };

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0) sendcount = 3;
	else if (myrank == 1) sendcount = 5;
	else if (myrank == 2) sendcount = 1;
	else if (myrank == 3) sendcount = 2;
	else if (myrank == 4) sendcount = 4;

	for (int i = 0; i < sendcount; i++) {
		isend[i] = sendcount;
	}
	
	MPI_Gatherv(isend, sendcount, MPI_INT, irecv, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	// MPI_Gather(&isend, 1, MPI_INT, irecv, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		for (int i = 0; i < 15; i++) {
			printf("%d\t", irecv[i]);
		}
		printf("\n");
	}

	MPI_Finalize();*/

	/*
	const int N = 100;
	int a[N];

	// MPI_Bcast(a, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// =======================================
	// ====== PARALLEL STATS =================

	int istart, iend, nprocs, myrank;
	int max = -999999, min = 999999, sum = 0, loc_min, loc_max;
	int isend_max[2], irecv_max[2], isend_min[2], irecv_min[2], sumrecv;

	// printf("max = %d\n\n", max);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int nlocal = N / nprocs;
	int nlast = N - N % nprocs - 1;
	// ======================================
    // ==== RNG =============================
	if (myrank == 0) {
		srand((unsigned)time(NULL));

		for (int i = 0; i < N; i++) {
			a[i] = rand() % 20001 - 10000;
			// printf("p%d: a[%d] = %d\n", myrank, i, a[i]);
		}

		printf("\n ::::: p%d : N = %d, nlocal = %d, nlast = %d ::::\n\n", myrank, N, nlocal, nlast);
	}

	MPI_Bcast(a, N, MPI_INT, 0, MPI_COMM_WORLD);

	istart = myrank * nlocal;
	iend = (myrank < nprocs - 1 ? istart + nlocal - 1 : N - 1);
	// printf("p%d : istart = %d, iend = %d\n", myrank, istart, iend, nlocal);

	for (int i = istart; i <= iend; i++) {
		if (a[i] > max) {
			max = a[i];
			loc_max = i;
		}
		if (a[i] < min) {
			min = a[i];
			loc_min = i;
		}
		sum += a[i];
	}

	isend_max[0] = max;
	isend_max[1] = loc_max;

	isend_min[0] = min;
	isend_min[1] = loc_min;

	MPI_Reduce(isend_max, irecv_max, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
	MPI_Reduce(isend_min, irecv_min, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sum, &sumrecv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		printf("p%d : celkove max: %d na pozicii: %d\n", myrank, irecv_max[0], irecv_max[1]);
		printf("p%d : celkove min: %d na pozicii: %d\n", myrank, irecv_min[0], irecv_min[1]);
		printf("p%d : celkovy priemer: %.4lf \n", myrank, (double)(sumrecv / N));
	}

	MPI_Finalize();*/

	const int N = 10;
	double A[N][N], b[N], c[N];

	int istart, iend, nprocs, myrank;
	double sum = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int nlocal = N / nprocs;
	int nlast = N - N % nprocs - 1;
	// ======================================
	// ==== RNG =============================
	if (myrank == 0) {
		srand((unsigned)time(NULL));

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				A[i][j] = rand() % 20 - 10;
			}
			b[i] = rand() % 100;
		}

		// printf("\n ::::: p%d : N = %d, nlocal = %d, nlast = %d ::::\n\n", myrank, N, nlocal, nlast);
	}

	MPI_Bcast(&A, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	istart = myrank * nlocal;
	iend = (myrank < nprocs - 1 ? istart + nlocal - 1 : N - 1);

	double *tempv;

	tempv = new double [iend - istart];
	// recvv = new double[iend - istart];

	for (int i = istart; i <= iend; i++) {
		tempv[i - istart] = 0;
		for (int k = 0; k < N; k++) {
			tempv[i - istart] += A[i][k] * b[k];
		}
		// printf("p%d: c[%d] = %lf\n", myrank, i, tempv[i - istart]);
	}

	MPI_Gather(tempv, iend - istart, MPI_DOUBLE, c, iend - istart, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		for (int i = 0; i < N; i++) {
			printf("%lf\n", c[i]);
		}
	}

	MPI_Finalize();
}