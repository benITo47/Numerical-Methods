#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void wykonajDlaN(int n, double xmin, double xmax, FILE* fout1, FILE* fout2);
double ZeraCzebyszewa(double xmin, double xmax, int n, int i);
double interpolacjaLagrange(double x, int n, double * tabX, double * tabY);

double func(double x) 
{
    return exp(-x*x); // f(x) = e^(-x^2) 
}

int main()
{
    double x1 = -5.0; 
    double x2 = 5.0;

    //---------UNIFORM-------------
    FILE* plik1 = fopen("plik1.dat", "w");
    //---------CZEBYSZEW-----------
    FILE* plik2 = fopen("plik2.dat", "w");

    for(int n = 5; n < 25; n+=5)
    {
        wykonajDlaN(n,x1,x2,plik1,plik2);
    }


    //---------UNIFORM-------------
    fclose(plik1);
    //---------CZEBYSZEW-----------
    fclose(plik2);
}




double ZeraCzebyszewa(double xmin, double xmax, int n, int i)
{
    return ((xmax - xmin) * cos(3.141 * (double) (2 * i + 1) / (double) (2 * n + 2)) + (xmin+ xmax)) / 2.0;
}

double interpolacjaLagrange(double x, int n, double * tabX, double * tabY)
{
    double suma = 0.0, iloczyn = 1.0;

    for (int j = 0; j <= n; j++)
    {
        iloczyn = 1.0;
        
        for (int k = 0; k <= n; k++){
            if (j != k){
                iloczyn *= (x - tabX[k]) / (tabX[j] - tabX[k]);  // (x-x_k)/(x_j - x_k) 
            }
        }
        suma += tabY[j] * iloczyn;
    }

    return suma;
}


void wykonajDlaN(int n, double xmin, double xmax, FILE* fout1, FILE* fout2)
{
    double h = (xmax-xmin)/n;
    

    //
    //----------- tworzymy kontenery na X i Y -----------
    //

    //---------UNIFORM------------- 
    double tabX[n+1];
    double tabY[n+1];

    //---------CZEBYSZEW-----------
    double tabX_czybyszew[n+1];
    double tabY_czybyszew[n+1];


    //
    //-----------Wyeplniamy Kontenery wartosciami-----------
    //
    for (int i = 0; i <= n; i++)
    {
        //---------UNIFORM-------------
        tabX[i] = xmin + i*h;
        tabY[i] = func(tabX[i]);

        //---------CZEBYSZEW-----------
        tabX_czybyszew[i] = ZeraCzebyszewa(xmin, xmax, n, i);
        tabY_czybyszew[i] = func(tabX_czybyszew[i]);
    }


    //
    // ----------- Obliczamy zinterpolowana wspolrzedna Y i zapisujemy do PLIKU-----------
    //

    for (double x = xmin; x <= xmax; x += 0.01)
    {
        //---------UNIFORM-------------
        fprintf(fout1, "%g\t%g\n", x, interpolacjaLagrange(x, n, tabX, tabY));
        //---------CZEBYSZEW-----------
        fprintf(fout2, "%g\t%g\n", x, interpolacjaLagrange(x, n, tabX_czybyszew, tabY_czybyszew));
    }

    //---------UNIFORM-------------
    fprintf(fout1, "\n\n");
    
    //---------CZEBYSZEW-----------
    fprintf(fout2, "\n\n");

}