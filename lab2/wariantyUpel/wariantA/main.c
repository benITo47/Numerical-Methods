#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/ludcmp.c"
#include "/home/NR/numerical_recipes.c/lubksb.c"

#define N 3 // rozmiar macierzy M: NxN



void print_matrix(float ** A);
void fprint_matrix(float **A, FILE *fout);
void invert_matrix(float **Matrix, int *indx, float* col, float** invertedMatrix);
 float matrix_norm(float ** M);

 void matrix_multiply(float** M1, float** M2, float** resM);

int main(void)
{
    float **A = NULL;
    float **B = NULL;

    A = matrix(1, N, 1, N);
    B = matrix(1,N,1,N);

    float value = 1; 
    int* indA = ivector(1,N);
    int* indB = ivector(1,N);
    float dA; 
    float dB; 

    float* a = vector(1,N);
    float* b = vector(1,N);

    for(int i=0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A[i+1][j+1] = value;
            B[i+1][j+1] = value;
            value++;
        }
    }
    B[1][1] = 1.1; 
    //rozklad LU dla A i B; 
    ludcmp(A, N, indA, &dA);
    ludcmp(B, N, indB, &dB);

    FILE* fout = fopen("wynikAzad1.txt", "w+");
    fprintf(fout, "Rozklad LU macierzy A: \n");
    fprint_matrix(A,fout);
    fclose(fout);


    fout = fopen("wynikBzad1.txt", "w+");
    fprintf(fout, "Rozklad LU macierzy B: \n");
    fprint_matrix(B,fout);
    fclose(fout);

    print_matrix(A);
    printf("\n\n");
    print_matrix(B);


    //zad 2 macierz odwrotna: 

    float **invA, **invB, *colA, *colB;
	invA = matrix(1, N, 1, N);
	invB = matrix(1, N, 1, N);
	colA = vector(1, N);
    colB = vector(1, N);

    invert_matrix(A, indA, colA, invA);
    invert_matrix(B, indB, colB, invB);
    
    fout = fopen("wynikAinverted_zad2.txt", "w+");
    fprintf(fout, "Macierz odwrotna A: \n");
    fprint_matrix(invA,fout);
    fclose(fout);
    

    fout = fopen("wynikBinverted_zad2.txt", "w+");
    fprintf(fout, "Macierz odwrotna B: \n");
    fprint_matrix(invB,fout);
    fclose(fout);
    
    
    print_matrix(invA);
    printf("\n\n");
    print_matrix(invB);

    //zad 3

    float** A_copy = matrix(1,N,1,N);
    float** B_copy = matrix(1,N,1,N);

    value = 1; 
    for(int i=0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A_copy[i+1][j+1] = value;
            B_copy[i+1][j+1] = value;
            value++;
        }
    }
    B_copy[1][1] = 1.1; 
    
    float wskUwarA = matrix_norm(A_copy) * matrix_norm(invA);
    float wskUwarB = matrix_norm(B_copy) * matrix_norm(invB);

    printf("Wskaznik uwarunkowania macierzy A: %g , norma A: %g, norma invertedA: %g\n", wskUwarA,matrix_norm(A_copy), matrix_norm(invA) );
    printf("Wskaznik uwarunkowania macierzy B: %g , norma B: %g, norma invertedB: %g\n", wskUwarB,matrix_norm(B_copy), matrix_norm(invB) );


    float **AA = matrix(1,N,1,N);
    float **BB = matrix(1,N,1,N);

    matrix_multiply(A_copy, invA, AA);
    matrix_multiply(B_copy, invB, BB);

    fout = fopen("iloczyn_AA.txt", "w+");
    fprintf(fout, "iloczyn macierzy A invA: \n");
    fprint_matrix(AA,fout);
    fclose(fout);
    

    fout = fopen("iloczynBB.txt", "w+");
    fprintf(fout, "iloczyn macierzy B invB: \n");
    fprint_matrix(BB,fout);
    fclose(fout);

    print_matrix(AA);
    printf("\n\n");
    print_matrix(BB);

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

void fprint_matrix(float **A, FILE *fout)
{
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
			fprintf(fout, "%5.4g\t\t\t", A[i][j]);
		fprintf(fout, "\n");
	}
	fprintf(fout, "\n");
}

void invert_matrix(float **Matrix, int *indx, float* col, float** invertedMatrix)
{
    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            col[i+1] = 0.0; 
        }
        col[j+1] = 1.0; 
        lubksb(Matrix, N, indx, col);
        for(int i = 0; i <N; i++)
        {
            invertedMatrix[i+1][j+1] = col[i+1];
        }
    }
}

 float matrix_norm(float ** M)
 {
 	float biggest = 0.0;
 	for (int i = 1; i <= N; ++i)
 		for (int j = 1; j <= N; ++j)
 			if (biggest < fabs(M[i][j]))
 				biggest = fabs(M[i][j]);		
 	return biggest;
 }

 void matrix_multiply(float** M1, float** M2, float** resM)
 {
 	for (int i = 1; i <= N; ++i)
 	{
 		for (int j = 1; j <= N; ++j)
 		{
 			resM[i][j] = 0.0;
 			for	(int k = 1; k <= N; k++)
 				resM[i][j] += M1[i][k] * M2[k][j];
 		}
 	}
 }