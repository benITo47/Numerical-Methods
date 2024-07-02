#include <stdio.h>
#include <stdlib.h>
#include <math.h>  // PRZY KOMPILACJI NALEZY DODAC FLAGE -lm



#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/tred2.c"
#include "/home/NR/numerical_recipes.c/tqli.c"
#include "/home/NR/numerical_recipes.c/pythag.c"

int main(void)
{
    // 1 wartosci stalych
    int nx = 20;
    int ny = 20;
    int n = nx * ny;
    int m = 10;
    float t = -0.021;


    // 2 Tworzymy macierze oraz wektory
    float** H = matrix(1, n, 1, n);
    float** Y = matrix(1, n, 1, n);
    float** X = matrix(1, n, 1, n);


    float* d = vector(1, n);
    float* e = vector(1, n);


    int* indx = ivector(1, n);

    // 3 Wypleniamy macierz H 
    int l;
    for(int i=1; i<=nx; i++)
    {
        for(int j=1; j<=ny; j++)
        {
            l = j + (i-1) * ny;
            for (int k=1; k<=n; k++) H[l][k] = 0.;

            if (i>1) H[l][l-ny] = t; //dla i=1 nie ma sasiada z lewej strony
            if (i<nx) H[l][l+ny] = t; //dla i=nx nie ma sasiada z prawej strony
           
            H[l][l] = -4*t;
            
            if (j>1) H[l][l-1] = t; //dla j=1 nie ma sasiada ponizej siatki
            if (j<ny) H[l][l+1] = t; //dla j=ny nie ma sasiada powyzej siatki
        }
    }

    //Wypelnienie diagonali macierzy Y, reszta zerami,
    //Wypelnienie macierzy X zerami 

    for (int i = 1; i <= n; i++)
    {
        for(int j = 1; j <= n; j++)
        {
            Y[i][j] = (i == j) ? 1. : 0.0f;
            X[i][j] = 0.0f;
        }
    }

    //wypelnienie wektorw e i d 

    for(int i = 1; i <= n; i++)
        d[i] = e[i] = 0.;

    // 4 Przeksztacenie do postaci trjdiagonalnej
    tred2(H, n, d, e);


    // 5 Diagonalizacja macierzy T 
    tqli(d, e, n, Y);

    // 6 Wektory wlasne 
    for(int i=1; i<=n; i++)
    {
        for(int j=1; j<=n; j++)
        {
            for (int k=1; k <= n; k++)
            
            //   x_k = P  * y_k
                X[i][j] += H[i][k] * Y[k][j];
        }
    }

    // 7 Sortowanie wektorow

    for(int l=1;l<=n;l++) 
        indx[l]=l; // inicjalizacja

    float e1, e2;
    int l1, l2;
    for(int l=1;l<=n-1;l++)
    {
        for(int k=n;k>=l+1;k--)
        {
            e1=d[k-1];
            e2=d[k];
            l1=indx[k-1];
            l2=indx[k];
            if(e2<e1)
            { //wymieniamy energie i indeksy wektorów miejscami
                d[k]=e1;
                d[k-1]=e2;
                indx[k]=l1;
                indx[k-1]=l2;
            }
        }
    }

    // 8 Zapis wektorow do pliku 
    FILE *fp;
    fp=fopen("dane.dat","w");
    for(int i=1;i<=nx;i++)
    {
        for(int j=1;j<=ny;j++)
        {
            l=j+(i-1)*ny;
            fprintf(fp,"%6d %6d ",i,j);
            for(int k=1;k<=m;k++)
                fprintf(fp," %12.6f ", X[l][indx[k]]);
            fprintf(fp,"\n");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    //9 wykresy w gnuplot: 


    // 10  zapis wartosci wlasnych; 
    FILE *fp2;

    fp2=fopen("dane_m.txt","w");

    for (int i = 1; i <= 10; i++)
        fprintf(fp2, "%d\t%g\t%d\n", i, d[i], indx[i]);

    fclose(fp2);

    free_matrix(H,1,n,1,n);
    free_matrix(Y,1, n, 1, n);

    free_matrix(X,1, n, 1, n);


    free_vector(d,1, n);
    free_vector(e,1, n);


    free_ivector(indx,1, n);

}