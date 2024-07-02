#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/ludcmp.c"
#include "/home/NR/numerical_recipes.c/lubksb.c"

#define N 4 // rozmiar macierzy M: NxN

void print_matrix(float **A)
{
    for (int i = 1; i <= N; ++i)
    {
        for (int j = 1; j <= N; ++j)
            printf("%4.2f\t", A[i][j]);

        printf("\n");
    }
    printf("\n");
}

int main(void)
{
    float **A = NULL;
    float **B = NULL;

    A = matrix(1, N, 1, N);
    
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            A[i][j] = (1.0 / (i + j));
        }
    }

    print_matrix(A);

    int *indA = ivector(1, N);
    float d;

   ludcmp(A, N, indA, &d);

    print_matrix(A);

    FILE* plik; 
    plik = fopen("Diagonala U", "w+");
    for(int i = 1; i <=N; i++)
    {
        fprintf(plik, "A[%d]", i);
        fprintf(plik, "[%d]", i);
        fprintf(plik, " = %f \n", A[i][i]);      
    }
    fclose(plik);

    
}