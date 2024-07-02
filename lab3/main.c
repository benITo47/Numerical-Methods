#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"

#define N 2000

void policz_i_zapisz(float* B, float* D0, float* D1, float*D2, float* Xs, float* Xn,float Beta, float F0, float Omega_, const char* fileName)
{
    //Stale dla kazdego przypadku
    double V0 = 0, x0 = 1, omega = 1.0, h = 0.02;

    //Inicjalizacja wektorow;
   
    long double suma_n = 0.0;
    long double suma_s = 0.0;
    //Obliczanie wartosci na podstawie otrzymanych argumentow
    float F0hh = F0*h*h;
    float a2 = (omega*omega)*(h*h) - 2 - Beta*h;
    float a3 = 1 + Beta*h;
    float a1 = 1.0f;

    //Wypelnianie dwoch pierwszych elementow wektorow
    B[0] = 1.0f;
    B[1] = 0.0f;

    D0[0] = 1.0f;
    D0[1] = 1.0f;
    
    D1[0] = 0.0f;
    D1[1] = -1.0f;
    
    D2[0] = 0.0f;
    D2[1] = 0.0f;

    Xs[0] = 1.0f;
    Xs[1] = 1.0f;

    //Wypelnianie reszty wektorow

    for (int i = 2; i <= N; i++)
    {
        B[i] = F0hh*sin(Omega_*h*i);
        D0[i] = a3;
        D1[i] = a2;
        D2[i] = a1;
        Xs[i] = 1.0;
    }


    int iterator = 0;
    
    while (iterator++ < 10000)
    {
        
        Xn[0] = B[0]/D0[0];
        Xn[1] = (B[1] - D1[1]*Xs[0]) / D0[1];

        for (int i = 2; i <= N; i++)
            Xn[i] = (B[i] - D1[i]*Xs[i-1] - D2[i]*Xs[i-2]) / D0[i];

        suma_n = 0.0;
        suma_s = 0.0;

        for(int i=0; i <=N; i++)
        {
            suma_n += Xn[i] * Xn[i];
            suma_s += Xs[i] * Xs[i];
        }

        if(fabs(suma_s-suma_n) < 1e-6)
        {
            break;
        }
        
        for (int i = 0; i <= N; i++)
          Xs[i] = Xn[i];
        

    }   
printf("%lf, %lf, \t iteracja: %d \n", suma_n, suma_s, iterator);

    //Zapis do pliku

    FILE * fout = fopen(fileName, "w");

    fprintf(fout, "%d\n\n", iterator);

    for (int i = 0; i <= N; i++)
        fprintf(fout, "%g\t%g\n", i*h, Xn[i]);

    fclose(fout);

}

int main(void)
{   
    float* B = vector(0, N); 
    float* D0 = vector(0, N);
    float* D1 = vector(0, N); 
    float* D2 = vector( 0, N); 
    float* Xs = vector( 0, N);
    float* Xn = vector( 0, N);
    srand(time(NULL));
    // Przypadek 1:
    float Beta = 0.0;
    float F0 = 0.0;
    float Omega_ = 0.8;
    char* fileName = "output1.txt";
    policz_i_zapisz(B,D0,D1,D2,Xs,Xn,Beta,F0,Omega_,fileName);
    //Przypadek 2:
    Beta = 0.4; 
    F0 = 0.0;
    Omega_ = 0.8; 
    fileName = "output2.txt";
   policz_i_zapisz(B,D0,D1,D2,Xs,Xn,Beta,F0,Omega_,fileName);
    //Przypadek 3:
    Beta = 0.4; 
    F0 = 0.1;
    Omega_ = 0.8; 
    fileName = "output3.txt";
    policz_i_zapisz(B,D0,D1,D2,Xs,Xn,Beta,F0,Omega_,fileName);

	return 0;
}


