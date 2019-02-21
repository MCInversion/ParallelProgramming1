#include "pch.h"

int main(int argc, char **argv) {
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

	MPI_Finalize();
}