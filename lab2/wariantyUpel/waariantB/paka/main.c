#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"

#define N 1000 // rozmiar macierzy M: NxN


void print_matrix(float ** A);
// void fprint_matrix(float ** A, FILE * fout);
// void invert_matrix(float ** M, int * indx, float * col, float ** invM);
// float matrix_norm(float ** M);
// void scalar_multiply_matrix(float ** A, float ** B, float ** mulM);

int main(void)
{
	float **A, *b;
	int m = 5;
	float x[3][3] = {{1,1,1}, {2,2,2}, {3,3,3}};
	float (*px)[][];

	//	alokacja macierzy
	A = matrix(1, N, 1, N);
	b = vector(1, N);
	// x = vector(1, N);

	// 	wypelnienie macierzy
	for (int i = 0; i < N; ++i)
	{
		b[i+1] = i+1;
		// x[i+1] = i+1;
		for (int j = 0; j < N; ++j)
			A[i+1][j+1] = (abs(i-j) > m) ? 0 : 1.0 / (1.0 + abs(i-j));
	}
	


	// // rozklad LU
	// ludcmp(A, N, indxA, &dA);
	// ludcmp(B, N, indxB, &dB);

	// // wynik
	// printf("LU decomposition of A:\n");
	// print_matrix(A);
	// printf("LU decomposition of B:\n");
	// print_matrix(B);

	// // odwrocenie macierzy
	// float **invA, **invB, *col;
	// invA = matrix(1, N, 1, N);
	// invB = matrix(1, N, 1, N);
	// col = vector(1, N);

	// invert_matrix(A, indxA, col, invA);
	// invert_matrix(B, indxB, col, invB);

	// // wynik
	// printf("Inverted A:\n");
	// print_matrix(invA);
	// printf("Inverted B:\n");
	// print_matrix(invB);

	// // wskazniki uwarunkowania macierzy
	// float condA, condB;
	// condA = matrix_norm(A_copy) * matrix_norm(invA);
	// condB = matrix_norm(B_copy) * matrix_norm(invB);

	// // wyniki
	// printf("Wskaznik uwarunkowania macierzy A: %g\n", condA);
	// printf("Wskaznik uwarunkowania macierzy B: %g\n\n", condB);

	// //iloczyny macierzy
	// float **mulA, **mulB;
	// mulA = matrix(1, N, 1, N);
	// mulB = matrix(1, N, 1, N);

	// scalar_multiply_matrix(A_copy, invA, mulA);
	// scalar_multiply_matrix(B_copy, invB, mulB);

	// // wyniki
	// printf("Multiplied A*invA:\n");
	// print_matrix(mulA);
	// printf("Multiplied B*invB:\n");
	// print_matrix(mulB);

	// // wypisanie do pliku tekstowego
	// FILE * fout = fopen("wynik.txt", "w+");

	// fprintf(fout, "Wskaznik uwarunkowania macierzy A: %g\n", condA);
	// fprintf(fout, "Wskaznik uwarunkowania macierzy B: %g\n\n", condB);

	// fprintf(fout, "Iloczyn A*A^-1:\n");
	// fprint_matrix(mulA, fout);
	// fprintf(fout, "Iloczyn B*B^-1:\n");
	// fprint_matrix(mulB, fout);

	// fclose(fout);

	// //	Zwolnienie pamieci
	// free_matrix(A, 1, N, 1, N);
	// free_matrix(B, 1, N, 1, 1);
	// free_matrix(A_copy, 1, N, 1, N);
	// free_matrix(B_copy, 1, N, 1, 1);
	// free_matrix(invA, 1, N, 1, N);
	// free_matrix(invB, 1, N, 1, N);
	// free_matrix(mulA, 1, N, 1, N);
	// free_matrix(mulB, 1, N, 1, 1);
	// free_ivector(indxA, 1, N);
	// free_ivector(indxB, 1, N);
	// free_vector(col, 1, N);

	// return 0;
}

void print_matrix(float ** A)
{
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
			printf("%4.2g\t", A[i][j]);
		
		printf("\n");
	}
	printf("\n");
}

// void fprint_matrix(float ** A, FILE * fout)
// {
// 	for (int i = 1; i <= N; ++i)
// 	{
// 		for (int j = 1; j <= N; ++j)
// 			fprintf(fout, "%5.4g\t\t\t", A[i][j]);
		
// 		fprintf(fout, "\n");
// 	}
// 	fprintf(fout, "\n");
// }

// void invert_matrix(float ** M, int * indx, float * col, float ** invM)
// {
// 	for (int j = 1; j <= N; j++)
// 	{
// 		for (int i = 1; i <= N; i++)
// 			col[i] = 0.0;
// 		col[j] = 1.0;

// 		lubksb(M, N, indx, col);

// 		for (int i = 1; i <= N; i++)
// 			invM[i][j] = col[i];
// 	}
// }

 float matrix_norm(float ** M)
 {
 	float biggest = 0.0;
 	for (int i = 1; i <= N; ++i)
 		for (int j = 1; j <= N; ++j)
 			if (biggest < fabs(M[i][j]))
 				biggest = fabs(M[i][j]);		
 	return biggest;
 }

 void scalar_multiply_matrix(float ** A, float ** B, float ** mulM)
 {
 	for (int i = 1; i <= N; ++i)
 	{
 		for (int j = 1; j <= N; ++j)
 		{
 			mulM[i][j] = 0.0;
 			for	(int k = 1; k <= N; k++)
 				mulM[i][j] += A[i][k] * B[k][j];
 		}
 	}
 }