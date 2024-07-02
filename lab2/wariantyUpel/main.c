#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//Sciezka taurus; 

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/ludcmp.c"
#include "/home/NR/numerical_recipes.c/lubksb.c"


#define N 3 

void matrix_fprint(float **A, FILE *fout);
void matrix_print(float ** A);
void calculateInvertedMatrix(float** Matrix, int* indx, float* wyrazywolne, float** invertedMatrix);
float matrix_norm(float** Matrix);
void matrix_multiply(float** M1, float** M2, float** resultM);
float matrix_det(float** M, float d );


int main(void)
{
    //macierze uzywane w zadaniu
    float** A = matrix(1,N,1,N);
    float** B = matrix(1,N,1,N);
    float** invA = matrix(1,N,1,N);
    float** invB = matrix(1,N,1,N);
    float** cpyA = matrix(1,N,1,N);
    float** cpyB = matrix(1,N,1,N);

    int* indA = ivector(1,N);
    int* indB = ivector(1,N);
    
    float* wyrazyWolneA = vector(1,N);
    float* wyrazyWolneB = vector(1,N);

    float dA, dB;

    float value = 1;

    for(int i =1; i<=N; i++)
    {
        for(int j=1; j <=N; j++)
        {
            A[i][j] = value;
            B[i][j] = value;
            cpyA[i][j] = value;
            cpyB[i][j] = value;
            value++;
        }   
    }

    B[1][1] = 1.1; 
    cpyB[1][1] = 1.1; 

    FILE* plik = fopen("Wyniki Lab 2.txt", "w+"); 
    

    //-------Wypis macierzy A i B---------------
    printf("Macierz A: \n");
        fprintf(plik,"Macierz A: \n");
        matrix_fprint(A,plik);
    matrix_print(A);
        fprintf(plik,"Macierz B: \n");
        matrix_fprint(B,plik);
    printf("Macierz B: \n");
    matrix_print(B);

    //
    //-------------Zad1 Rozklad LU---------------------
    //

    ludcmp(A,N,indA,&dA);
    ludcmp(B,N,indB,&dB);


    //---------Wypis rozkladu LU macierzy A i B-------------- 

    printf("Rozklad LU macierzy A: \n");
    matrix_print(A);

        fprintf(plik,"\nRozklad LU macierzy A: \n");
        matrix_fprint(A, plik);

    printf("Rozklad LU macierzy B: \n");
    matrix_print(B);

        fprintf(plik,"\nRozklad LU macierzy B: \n");
        matrix_fprint(B,plik);


    //Obliczanie wyznacznikow Macierzy

    float detA = matrix_det(A,dA);
    float detB = matrix_det(B,dB);

    printf("\nWyznacznik A: %e\n", detA);
    printf("\nWyznacznik B: %f\n", detB);

    fprintf(plik,"\nWyznacznik A: %e\n", detA);
    fprintf(plik,"\nWyznacznik B: %f\n", detB);

    //
    //-------Zad2 Macierze odwrotne A B------------
    //


    calculateInvertedMatrix(A, indA, wyrazyWolneA, invA);
    calculateInvertedMatrix(B, indB, wyrazyWolneB, invB);

    //----Wypis odwotnosci A i B---


    printf("Macierz invA: \n");
    matrix_print(invA);

    fprintf(plik,"\nMacierz invA: \n");
    matrix_fprint(invA, plik);
    

    printf("Macierz invB: \n");
    matrix_print(invB);

    fprintf(plik,"\nMacierz invB: \n");
    matrix_fprint(invB, plik);


    //
    // -----Zad3 Wskazniki uwarunkowania Macierzy A i B ------
    //

    float NormaA = matrix_norm(cpyA);
    float NormaInvA = matrix_norm(invA);
    float wskUwarA = NormaA * NormaInvA;

    float NormaB = matrix_norm(cpyB);
    float NormaInvB = matrix_norm(invB);
    float wskUwarB = NormaB* NormaInvB;

    printf("\nNorma A: %f\n, Norma invertedA: %e\n, Wskaznik Uwarunkowania: %e\n", NormaA, NormaInvA, wskUwarA);
    printf("\nNorma B: %f\n, Norma invertedB: %e\n, Wskaznik Uwarunkowania: %e\n", NormaB, NormaInvB, wskUwarB);

    fprintf(plik, "\nNorma A: %f", NormaA);
    fprintf(plik, "\nNorma invertedA: %e", NormaInvA);
    fprintf(plik, "\nWskaznik UwarunkowaniaA: %e", wskUwarA);
    fprintf(plik, "\nNorma B: %f", NormaB);
    fprintf(plik, "\nNorma invertedB: %e", NormaInvB);
    fprintf(plik, "\nWskaznik Uwarunkowania B: %e", wskUwarB);
    

    //
    //------Zad 4 iloczyn Macierzy A i inv A oraz B i invB-------
    //

    float** AA = matrix(1,N,1,N);
    float** BB = matrix(1,N,1,N);

    matrix_multiply(cpyA, invA, AA);
    matrix_multiply(cpyB, invB, BB);

    printf("Iloczyn macierzy A i invA:\n");
    matrix_print(AA);

        fprintf(plik,"\nIloczyn macierzy A i invA:\n");
        matrix_fprint(AA, plik);
    
    
    printf("Iloczyn macierzy B i invB:\n");
    matrix_print(BB);

        fprintf(plik,"Iloczyn macierzy B i invB:\n");
    matrix_fprint(BB,plik);

    fclose(plik);

}

void matrix_fprint(float **A, FILE *fout)
{
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
			fprintf(fout, "%5.4g\t\t\t", A[i][j]);
		fprintf(fout, "\n");
	}
	fprintf(fout, "\n");
}

void matrix_print(float ** A)
{
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
			printf("%4.2g\t", A[i][j]);
		
		printf("\n");
	}
	printf("\n");
}

void calculateInvertedMatrix(float** Matrix, int* indx, float* wyrazywolne, float** invertedMatrix)
{
    for(int j = 1; j <= N ; j++)
    {
        for(int i = 1; i <= N; i++)
        {
            wyrazywolne[i] = 0.0;
        }
        wyrazywolne[j] = 1.0; 
        lubksb(Matrix, N, indx, wyrazywolne);
        for(int i = 1; i <= N; i++)
        {
            invertedMatrix[i][j] = wyrazywolne[i];
        }

    }
}

float matrix_norm(float** Matrix)
{
    float biggest = 0.0; 
    for(int i = 1; i<= N; i ++)
    {
        for(int j =1; j<= N ; j++)
        {
            if(biggest < fabs(Matrix[i][j]))
                biggest = fabs(Matrix[i][j]);
        }
    }
    return biggest;
}

void matrix_multiply(float** M1, float** M2, float** resultM)
{
    for (int i = 1; i <= N; ++i)
 	{
 		for (int j = 1; j <= N; ++j)
 		{
 			resultM[i][j] = 0.0;
 			for	(int k = 1; k <= N; k++)
 				resultM[i][j] += M1[i][k] * M2[k][j];
 		}
 	}
}

float matrix_det(float** M, float d )
{   
   
    for(int i = 1; i <= N ; i++)
    {
        d = d * M[i][i];
    }
    return d;
}