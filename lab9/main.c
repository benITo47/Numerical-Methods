#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define frand() ((double)rand())/(RAND_MAX+1.0)


float funkcja(float x, float x_min, float x_max, float x_0, float o);
float f_szum(float y);
void approx(FILE * fout, int m, float * c, float * s, int n, float phi[51][201], float * x_k);

int main()
{
    // deklaracja zmiennych
    float x_min = -4., x_max = 4., x_0 = 2.0;
    float o = (x_max - x_min) / 16.;
    int n = 201;
    int m = 50;

    float x_k[n], y_k[n];
    float phi[m+1][n];


    float h = (x_max - x_min) / (n-1);




    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < m+1; j++)
            phi[j][i] = 0.;
        
        x_k[i] = x_min + h * i;
        y_k[i] = f_szum(funkcja(x_k[i], x_min, x_max, x_0, o));
    }

// 1
    float denominator, nominator, alpha = 0., beta = 0.;

    // rzad -1 i 0
    float phi1 = 0.;
    for(int i = 0; i < n; i++)
        phi[0][i] = 1.;

    // rzad 1
    for(int i = 0; i < n; i++)
    {     
        // alpha
        denominator = nominator = 0.0;
        for (int k = 0; k < n; k++)
        {
            nominator += x_k[k] * phi[0][k] * phi[0][k];
            denominator += phi[0][k] * phi[0][k];
        }
        alpha = nominator/denominator;
        
        // result
        phi[1][i] = (x_k[i] - alpha)*phi[0][i];
    }

    // reszta rzedow
    for (int j = 1; j < m; j++)
    {
        for (int i = 0; i < n; i++)
        {
            // alpha
            denominator = nominator = 0.;
            for (int k = 0; k < n; k++)
            {
                nominator += x_k[k] * phi[j][k] * phi[j][k];
                denominator += phi[j][k] * phi[j][k];
            }
            alpha = nominator/denominator;

            // beta
            denominator = nominator = 0.;
            for (int k = 0; k < n; k++)
            {
                nominator += x_k[k] * phi[j-1][k] * phi[j][k];
                denominator += phi[j-1][k] * phi[j-1][k];
            }
            beta = nominator/denominator;

            // result
            phi[j+1][i] = (x_k[i] - alpha)*phi[j][i] - beta*phi[j-1][i];
        }
    }

    // wyniki dla wielomianow
    FILE * fout1 = fopen("Gram.dat", "w");
    for (int k = 0; k < n; k++)
    {
        fprintf(fout1, "%g ", x_k[k]);
        for (int j = 0; j < 7; j++)
        {
            fprintf(fout1, "%g ", phi[j][k] / phi[j][0]);
        }
        fprintf(fout1, "\n");
    }

    fclose(fout1);

// 4
    FILE * fout2 = fopen("pkt.dat", "w");

    for (int i = 0; i < n; i++)
        fprintf(fout2, "%g\t%g\n", x_k[i], y_k[i]);
    
    fclose(fout2);
    
    FILE * fout3 = fopen("approx.dat", "w");

    float c[m+1], s[m+1];

    for (int i = 0; i < m+1; i++)
        c[i] = s[i] = 0.;

    for (int j = 0; j < 51; j++)
    {
        c[j] = s[j] = 0.;
    
        for (int i = 0; i < n; i++)
        {
            c[j] += y_k[i] * phi[j][i];
            s[j] += phi[j][i] * phi[j][i];
        }
    }

    approx(fout3, 10, c, s, n, phi, x_k);
    approx(fout3, 30, c, s, n, phi, x_k);
    approx(fout3, 50, c, s, n, phi, x_k);

    fclose(fout3);
}

void approx(FILE * fout, int m, float * c, float * s, int n, float phi[51][201], float * x_k)
{
    float F[n];

    for (int i = 0; i < m+1; i++)
        F[i] = 0.;

    for (int k = 0; k < n; k++)
    {
        F[k] = 0.;
        for (int j = 0; j < m+1; j++)
        {
            F[k] += c[j] / s[j] * phi[j][k];
        }
    }
    
    for (int i = 0; i < n; i++)
        fprintf(fout, "%g\t%g\n", x_k[i], F[i]);
    
    fprintf(fout, "\n\n");
}

float funkcja(float x, float x_min, float x_max, float x_0, float o)
{
    return sin(14*3.14*x / (x_max-x_min)) * (exp(-pow(x-x_0, 2) / (2*o*o)) + exp(-pow(x+x_0, 2) / (2*o*o)));
}

float f_szum(float y)
{
    return y + (frand() - 0.5) / 5.;
}