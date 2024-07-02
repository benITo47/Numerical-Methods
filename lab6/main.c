#include <stdio.h>
#include <math.h>

double licz_r(double * a, double * b, int stopien, double x0);

int main(void)
{
    int n = 5;
    double a[6] = {240, -196, -92, 33, 14, 1};  //wielomian do policzena
    double b[6] = {0};                          //tablica pomocnicza
    double c[6] = {0};                          //tablica pomocnicza
    double x0 = 0, x1 = 0;;
    int ITmax = 30, it = 0;
    int stopien = 0;
    double Rj = 0, Rj1 = 0;

    double miejscaZerowe[6] = {0};

    FILE * fout = fopen("out.txt", "w");
    fprintf(fout,"---------------------------Wartosc srartowa x0 = %d ----------------------------\n",x0);
    fprintf(fout,"-------------------------Wspolczynniki wielomianu a[n] -------------------------\n");
    for(int i = 0 ; i<6; i++)
    {
        fprintf(fout,"a[%d] = %g\n",i, a[i]);
    }

    for (int i = 1; i <= n; i++)
    {
        stopien = n-i+1;

        x0 = 0;

        fprintf(fout,"-------------------aktualny stopien wielomianu n = N-i+1 = %d-------------------\n",n-i+1);
        for (it = 1; it <= ITmax; it++)
        {
            Rj = licz_r(a, b, stopien, x0);
            Rj1 = licz_r(b, c, stopien-1, x0);

            x1 = x0 - Rj / Rj1;
            
            if (fabs(x1 - x0) < 1.0e-7)     //waurnek iteracyjny
                break;

            x0 = x1;
            fprintf(fout, "i=%d\titeration=%d\tx0=%.4g\tRj=%.4g\tRj1=%.4g\n", i, it, x0, Rj, Rj1);
        }

        fprintf(fout,"\n");
        miejscaZerowe[i] = x1; 


        for (int j = 0; j <= n-1; j++)
            a[j] = b[j];
    }

    fprintf(fout,"-------------------------Mijesca zerowe wielomianu a[n] ------------------------\n");
    for(int i = 1 ; i<6; i++)
    {
        fprintf(fout,"x0 = %g\n",i, miejscaZerowe[i]);
    }
    fclose(fout);
}

double licz_r(double*a, double*b, int stopien, double x0)
{
    b[stopien] = 0; //zerujemy element o najwiekszej potedze 

    for(int k = stopien - 1; k >= 0; k--)
    {
        b[k] = a[k+1] + x0 * b[k+1];    //wielomian stopnia mniejszego zapsiujemy w wektorze b; 
    }

    return a[0] + x0 * b[0];
}