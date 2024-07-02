#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Dyrektywy zakladajace, ze te trzy pliki sa skopiowane do aktualnego katalogu. */

// #include "nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
// #include "gaussj.c" // To tez mozna usunac, jesli plik jest dodany w poleceniu kompilacji.

/* Dyrektywy dla Taurusa (nie wymagaja kopiowania plikow, ale Taurus musi dzialac...) */
#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/gaussj.c"

#define N 300 // rozmiar macierzy M: NxN

int main(void)
{
    float **M = NULL;
    float **b = NULL;
    //	Alokacja macierzy
    M = matrix(1, N, 1, N);
    b = matrix(1, N, 1, 1);

    float h = 0.1;
    float V0 = 0;
    float A = 1;
    float Omega = 1;

    for (int i = 1; i <= N; i++)
    {
        b[i][1] = 0.0;
        for (int j = 1; j <= N; j++)
            M[i][j] = 0.0;
    }

    b[1][1] = A;
    b[2][1] = V0 * h;

    for (int i = 1; i <= N; i++)
    {
        M[i][i] = 1;
        if (i >= 3)
        {
            M[i][i - 1] = (Omega * h * h - 2);
            M[i][i - 2] = 1;
        }
    }

    M[2][1] = -1;

    /*
    for(int i = 1; i <= N; i++)          //Print Macierzy
     {
         for(int j = 1; j <= N; j++)
         {
             printf("%2.2f , ", M[i][j]);
         }
         printf("\n\n");
     }
   */
    gaussj(M, N, b, 1);

    //	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
    for (int i = 1; i <= N; ++i)
        printf("%g    %g  \n", i * h, b[i][1]);

    FILE *fptr;
    fptr = fopen("Dane.txt", "w");
    for (int i = 1; i <= N; ++i)
    {
        fprintf(fptr, "%f", i * h);
        fprintf(fptr, " %g \n", b[i][1]);
    }
    fclose(fptr);
    

    //	Zwolnienie pamieci
    free_matrix(M, 1, N, 1, N);
    free_matrix(b, 1, N, 1, 1);

    return 0;
}