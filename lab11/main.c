#include <stdio.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "/home/NR/numerical_recipes.c/nrutil.h"
#include "/home/NR/numerical_recipes.c/nrutil.c"
#include "/home/NR/numerical_recipes.c/four1.c"



void FFT_Odszumianie(int n, FILE * fout); //Funkcja przeprowadza faktycznie całe zadanie - generuje tablice, zaszumienie, oblicza transformaty, zapisuje wyniki do pliku 
float f0(float t, float T); // Sygnał wejsciowy niezaburzony
float f_rand ( const float min , const float max ); // delta - generuj delte z zakresu [min,max] 
float g(float t, float sigma); // Funkcja wagowa

int main()
{
    FILE * fout8 = fopen("k8.dat", "w");
    FILE * fout10 = fopen("k10.dat", "w");
    FILE * fout12 = fopen("k12.dat", "w");

    FFT_Odszumianie(pow(2, 8), fout8);
    FFT_Odszumianie(pow(2, 10), fout10);
    FFT_Odszumianie(pow(2, 12), fout12);

    

    fclose(fout8);
    fclose(fout10);
    fclose(fout12);
}

void FFT_Odszumianie(int n, FILE * fout)
{
    float T = 1.0;
    float t_max = 3.0*T;
    float sigma = T/20.0;
    float dt = t_max/(float)n;

    float* g_k = vector(1,2*n);
    float* f = vector(1, 2*n); // 2* n bo część rzeczywista i urojona 
    float* g1 = vector(1, 2*n); // 
    float* g2 = vector(1, 2*n); //  Urojona jest zerami dla kadego z trzech wektorow, ale potrzebujemy tych miejsc dla bibliotek 

    

    // 3.1 - Wypełnianie tablicy wartosciami 

    //Zerowanie vectorow
    for (int i = 1; i <= 2*n; i++)
        f[i] = g1[i] = g2[i] = 0.;
    
    for (int i = 1; i <= n; i++)
    {
        float t_i = dt * (i-1);
        f[2*i-1] = f0(t_i, T) + f_rand(-0.5,0.5); // TO jest sygnał początkowy wraz z zaszumieniem
        g1[2*i-1] = g2[2*i-1] = g(t_i, sigma);
        
        fprintf(fout, "%g\t%g\n", t_i, f[2*i-1]); // zapis sygnału zaburzonego 
    }

    fprintf(fout, "\n\n");

    // 3.2 Obliczanie f_k, g_1, g_2 
    four1(f, n, 1);
    four1(g1, n, 1);
    four1(g2, n, -1); // t < 0 

    // 3 Obliczanie f_k * (g_1(k) + g_2(k))


 //g[i] = g_1[i] + g_2[i] splot to suma transformat dla dwoch g (z t<0 i dla g gdzie t > 0) 
    for(int i =1 ; i <=n; i++)
    {
        g_k[2*i-1] = g1[2*i-1] + g2[2*i-1];
        g_k[2*i] = g1[2*i] + g2[2*i];
    }

//Obliczanie transformaty splotu 
    float a1 = 0., b1 = 0., a2 = 0., b2 = 0.;
    for (int i = 1; i <= n; i++)
    {
        a1 = f[2*i-1]; // Re{f(k_i)}
        b1 = f[2*i]; // Im {f(k_i)}

       

        a2 = g_k[2*i-1]; // Re{g(k_i)}
        b2 = g_k[2*i]; // Im{g(k_i)}

        f[2*i-1] = a1*a2 - b1*b2;
        f[2*i] = a1*b2 + a2*b1;
    }
    
    // 4 Transforamata odwrotna - wyliczamy splot f(t) * g(t)
    four1(f, n, -1);

    // 5 Szukanie maksimum w tablicu

    float f_max = abs(f[1]);
    for (int i = 1; i <= n; i++)
    {
        if (f_max < abs(f[2*i-1]))
            {
                f_max = abs(f[2*i-1]);
            }
    }


    // 6 Zapis do pliku 
    for (int i = 1; i <= n; i++)
        fprintf(fout, "%g\t%g\n", (i-1)*dt, f[2*i-1]*2.5/f_max); // zapis znormalizowanego splotu 

    free_vector(f, 1, 2*n);
    free_vector(g1, 1, 2*n);
    free_vector(g2, 1, 2*n);
    free_vector(g_k, 1, 2*n);
}




float f0(float t, float T)
{
    float omega = 2*M_PI/T; 
    return sin(omega * t) + sin(2* omega * t) + sin(3*omega * t);
}


float g(float t, float sigma)
{
    return exp(-t*t/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}

float f_rand ( const float min , const float max ){
 float r = ( float ) rand () / RAND_MAX ; // Przedzial [0 , 1]
 r = r * ( max - min ) + min ; // Przeskalowanie do [min , max]
 return r;
}

