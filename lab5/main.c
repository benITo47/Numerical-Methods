#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/tred2.c"
#include "/home/NR/numerical_recipes.c/tqli.c"
#include "/home/NR/numerical_recipes.c/pythag.c"

// Funkcja do mnożenia macierzy przez wektor
void pomnozMacierzWektor(float *v_to, float **macierz, float *v_from, int n);
// Funkcja do obliczania iloczynu skalarnego dwóch wektorów
float iloczynSkalarny(float *v_poziomo, float *v_pionowo, int n);
// Funkcja do dzielenia wszystkich elementów wektora przez daną wartość
void podzielWektor(float *vec, float val, int n);
// Funkcja do kopiowania zawartości jednego wektora do drugiego
void skopiujWektor(float *v_from, float *v_to, int n);
// Funkcja do redukcji macierzy przy użyciu metody Hotellinga
void redukcjaHotelling(float **macierz, float val, float *v,float* v2, int n);

int main(void)
{
    int n = 7;

    float **A = matrix(1, n, 1, n);
    float **W = matrix(1, n, 1, n); // Kopia macierzy A - do drugiej części zadania
    float *d = vector(1, n);        // Wektor przechowujący wartości własne macierzy
    float *e = vector(1, n);        // Wektor przechowujący wartości pomocnicze


    // 1. Wypełnienie macierzy A wartościami zgodnie z przepisem
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            W[i][j] = A[i][j] = sqrt(i + j);
            // Kopia macierzy A
        }
    }

    // Inicjalizacja wektorów d i e wartościami początkowymi
    for (int i = 1; i <= n; i++)
        d[i] = e[i] = 1.0f;

    // 2. Redukcja macierzy A do postaci trójdiagonalnej
    tred2(A, n, d, e);

    // 3. Wyznaczenie wartości własnych macierzy trójdiagonalnej
    tqli(d, e, n, A);

    // 4. Metoda iteracyjna wyznaczania wartości własnych
    float *x_stary = vector(1, n); // Wektor poprzedniego przybliżenia
    float *x_nowy = vector(1, n);  // Wektor aktualnego przybliżenia
    float *lambda = vector(1, n);  // Tablica przechowująca wartości własne

    FILE *fp1; // plik z zapisem danych
    fp1 = fopen("dane.txt", "w");
    FILE *fp2 = fopen("dane_wszystkie.txt", "w");
    fprintf(fp1, "i\t\tlambda - iteracyjna\t\t\td-nrecpies\n");
    fprintf(fp2, "k\ti\tlambda[i]\n\n");
    printf("k\t\ti\t\tlambda[i]\n\n");
    for (int k = 1; k <= n; k++)
    {
        for (int j = 1; j <= n; j++){
            x_stary[j] = 1;
        }

        for (int i = 1; i <= 8; i++)
        {
            pomnozMacierzWektor(x_nowy, W, x_stary, n); // x_{i+1} = W_k * x_i

            lambda[i] = iloczynSkalarny(x_nowy, x_stary, n) / iloczynSkalarny(x_stary, x_stary, n); // lambda_i = (x^T_{i+1} * x_i)/(x^T_i * x_i)

            //podzielWektor(x_nowy, sqrt(iloczynSkalarny(x_nowy, x_nowy, n)), n); // x_{i+1} = x_{i+1}/||x_{i+1}||_2

            skopiujWektor(x_stary,x_nowy, n); // x_i = x_{i+1}
        }
        redukcjaHotelling(W, lambda[8], x_stary, x_stary, n); // W_{k+1} = W_k - Lambda_k * x_k * x^T_k

        fprintf(fp1, "%d\t\t\t%g\t\t\t%g\n", k, lambda[8], d[n + 1 - k]);

        for (int i = 1; i <= 8; i++)
        {
            printf("%d\t%d\t%g\n", k, i, lambda[i]);
            fprintf(fp2, "%d\t%d\t%g\n", k, i, lambda[i]);
        }

        printf("\n");
        fprintf(fp2, "\n");
    }

    fclose(fp1);

    // Zwolnienie zaalokowanej pamięci
    free_matrix(A, 1, n, 1, n);
    free_matrix(W, 1, n, 1, n);
    free_vector(d, 1, n);
    free_vector(e, 1, n);
    free_vector(x_stary, 1, n);
    free_vector(x_nowy, 1, n);
    free_vector(lambda, 1, n);
}

// Funkcja mnoży macierz przez wektor
void pomnozMacierzWektor(float *v_to, float **macierz, float *v_from, int n)
{
    for (int i = 1; i <= n; i++)
    {
        float suma = 0;

        for (int j = 1; j <= n; j++)
            suma += macierz[i][j] * v_from[j];

        v_to[i] = suma;
    }
}

// Funkcja obliczająca iloczyn skalarny dwóch wektorów
float iloczynSkalarny(float *v_poziomo, float *v_pionowo, int n)
{
    float wynik = 0;

    for (int i = 1; i <= n; i++)
        wynik += v_poziomo[i] * v_pionowo[i];

    return wynik;
}

// Funkcja dzieląca wszystkie elementy wektora przez daną wartość
void podzielWektor(float *vec, float val, int n)
{
    for (int i = 1; i <= n; i++)
        vec[i] /= val;
}

// Funkcja kopiująca zawartość jednego wektora do drugiego
void skopiujWektor(float *v_to, float *v_from, int n)
{

    for (int i = 1; i <= n; i++)
        v_to[i] = v_from[i];
}
// Funkcja redukująca macierz
void redukcjaHotelling(float **macierz, float lambda, float *v,float* v2, int n)
{
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            macierz[i][j] -= lambda * v[i] * v2[j];
}
